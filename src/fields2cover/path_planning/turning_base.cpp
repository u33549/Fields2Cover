//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>
#include "fields2cover/path_planning/turning_base.h"

namespace f2c::pp {

namespace {

// Perpendicular distance from a point to the swath that runs from pos in
// direction ang, for the first strip_length of it only.
//
// The strip is not endless. What excuses a turn for clipping it is that the
// vehicle is settling onto that swath either side of the turn, and settling
// takes about a turning radius; a path that runs down the swath line far
// beyond that is driving over standing crop, however well aligned it is.
// Measured before this bound: Reeds-Shepp put a whole u-turn under the
// boundary and half of it was excused.
double distanceToSwath(const F2CPoint& p, const F2CPoint& pos, double ang,
    double strip_length) {
  const double dx = p.getX() - pos.getX();
  const double dy = p.getY() - pos.getY();
  const double ahead = dx * std::cos(ang) + dy * std::sin(ang);
  if (ahead > strip_length) {
    return std::numeric_limits<double>::infinity();
  }
  if (ahead < 0.0) {
    return std::hypot(dx, dy);  // behind the swath end: measure to the end
  }
  return std::fabs(-dx * std::sin(ang) + dy * std::cos(ang));
}

// The point test runs on every sample of every candidate turn -- on the order
// of a hundred million times for one field -- so it goes through the prepared
// area. The depth below only runs on samples that already fell outside, which
// is rare, so that one still asks the geometry.
void measureOutside(const F2CPath& path, const F2CCells& free_space,
    const f2c::PreparedArea& free_area,
    double step, double half_swath, double strip_length,
    const F2CPoint& in_pos, double in_ang,
    const F2CPoint& out_pos, double out_ang,
    double* length_out, double* deepest_out, double* in_swath_out,
    double abort_over = std::numeric_limits<double>::infinity()) {
  // abort_over lets a caller that only asks "does it fit?" stop at the first
  // metre that says no. Measuring a point against the ground is the expensive
  // step, and a candidate that leaves it early leaves it whatever the rest of
  // the path does.
  *length_out = 0.0;
  *deepest_out = 0.0;
  *in_swath_out = 0.0;
  const auto& states = path.getStates();
  for (size_t i = 0; i + 1 < states.size(); ++i) {
    const F2CPoint a = states[i].point;
    const F2CPoint b = states[i + 1].point;
    const double len = a.distance(b);
    if (len <= 0.0) {
      continue;
    }
    const int n = std::max(1, static_cast<int>(std::ceil(len / step)));
    for (int j = 0; j < n; ++j) {
      const double t = (j + 0.5) / n;
      const F2CPoint p {a.getX() + (b.getX() - a.getX()) * t,
                        a.getY() + (b.getY() - a.getY()) * t};
      if (free_area.holds(p.getX(), p.getY())) {
        continue;
      }
      if (half_swath > 0.0 &&
          std::min(distanceToSwath(p, in_pos, in_ang, strip_length),
                   distanceToSwath(p, out_pos, out_ang, strip_length))
              <= half_swath) {
        *in_swath_out += len / n;
        continue;
      }
      *length_out += len / n;
      if (*length_out > abort_over) {
        return;
      }
      *deepest_out = std::max(*deepest_out, p.distance(free_space));
    }
  }
}

// Smallest radius the path actually turns at. A turn planner keeps its own
// limit, but a candidate is only worth keeping if the joined path does too.
double minRadius(const F2CPath& path) {
  double r = std::numeric_limits<double>::infinity();
  const auto& states = path.getStates();
  for (size_t i = 0; i + 1 < states.size(); ++i) {
    const double ds = states[i].point.distance(states[i + 1].point);
    const double da = F2CPoint::getAngleDiffAbs(
        states[i].angle, states[i + 1].angle);
    if (da > 1e-6 && ds > 1e-9) {
      r = std::min(r, ds / da);
    }
  }
  return r;
}

}  // namespace


std::vector<double> TurningBase::transformToNormalTurn(
    const F2CPoint& start_pos, double start_angle,
    const F2CPoint& end_pos, double end_angle) {
  double dist_start_end = start_pos.distance(end_pos);
  auto dir = end_pos - start_pos;
  auto angle = F2CPoint::mod_2pi(dir.getAngleFromPoint());
  bool inverted {false};
  start_angle = F2CPoint::mod_2pi(start_angle - angle);
  end_angle = F2CPoint::mod_2pi(end_angle - angle);
  if (start_angle > boost::math::constants::pi<double>()) {
    start_angle = F2CPoint::mod_2pi(-start_angle);
    end_angle = F2CPoint::mod_2pi(-end_angle);
    inverted = true;
  }
  return {dist_start_end, angle, start_angle, end_angle,
    static_cast<double>(inverted)};
}


F2CPath TurningBase::plainTurn(const F2CRobot& robot,
    const F2CPoint& start_pos, double start_angle,
    const F2CPoint& end_pos, double end_angle) {
  auto turn_values =
    transformToNormalTurn(start_pos, start_angle, end_pos, end_angle);
  double dist_start_end = turn_values[0];
  double rot_angle = turn_values[1];
  double start_angle_t = turn_values[2];
  double end_angle_t = turn_values[3];
  double inverted = turn_values[4];

  F2CPath path;
  if (using_cache) {
    path = createTurnIfNotCached(robot,
        dist_start_end, start_angle_t, end_angle_t);
  } else {
    path = createSimpleTurn(robot,
        dist_start_end, start_angle_t, end_angle_t);
  }
  if (path.size() <= 1) {return F2CPath();}

  if (inverted) {
    std::for_each(path.begin(), path.end(), [] (auto& s) {
        s.point.setY(-s.point.getY());
        s.angle = F2CPoint::mod_2pi(-s.angle);});
  }
  for (auto&& s : path) {
    s.point = F2CPoint(.0, .0).rotateFromPoint(rot_angle, s.point) + start_pos;
    s.angle = F2CPoint::mod_2pi(s.angle + rot_angle);
  }

  correctPath(path, start_pos, end_pos);
  return path;
}

void TurningBase::correctPath(F2CPath& path, const F2CPoint& start_pos,
    const F2CPoint& end_pos, float max_error_dist) {
  if (path.size() < 2) {return;}

  auto is_near = [max_error_dist](const F2CPoint& a, const F2CPoint& b) {
      return (a.distance(b) < max_error_dist);
  };
  if (is_near(path[0].point, start_pos)) {
    path[0].point = start_pos;
  }
  if (is_near(path.back().point, end_pos)) {
    path.back().point = end_pos;
  }
}



std::vector<std::pair<F2CPoint, double>> TurningBase::concaveCorners(
    const F2CCells& region, double eps) {
  // Two point tests per corner, and a field's ground has hundreds of them.
  const f2c::PreparedArea area(region);
  std::vector<std::pair<F2CPoint, double>> corners;
  for (auto&& cell : region) {
    for (auto&& ring : cell) {
      std::vector<F2CPoint> v;
      for (auto&& p : ring) {
        v.emplace_back(p);
      }
      if (v.size() > 1 && v.front().distance(v.back()) < 1e-9) {
        v.pop_back();
      }
      const size_t n = v.size();
      if (n < 3) {
        continue;
      }
      for (size_t i = 0; i < n; ++i) {
        const F2CPoint& prev = v[(i + n - 1) % n];
        const F2CPoint& here = v[i];
        const F2CPoint& next = v[(i + 1) % n];
        double ax = prev.getX() - here.getX();
        double ay = prev.getY() - here.getY();
        double bx = next.getX() - here.getX();
        double by = next.getY() - here.getY();
        const double la = std::hypot(ax, ay);
        const double lb = std::hypot(bx, by);
        if (la < 1e-9 || lb < 1e-9) {
          continue;
        }
        ax /= la; ay /= la; bx /= lb; by /= lb;
        double mx = ax + bx;
        double my = ay + by;
        const double lm = std::hypot(mx, my);
        if (lm < 1e-6) {
          continue;  // straight vertex: no bisector to speak of
        }
        mx /= lm; my /= lm;
        // The bisector of the *narrow* side is m. Step both ways: a convex
        // corner has the region on the narrow side, a concave one on the wide
        // side. Asking the region rather than the winding keeps this immune to
        // whichever orientation the ring came back with.
        const bool narrow_in = region.isPointIn(
            F2CPoint(here.getX() + eps * mx, here.getY() + eps * my));
        const bool wide_in = region.isPointIn(
            F2CPoint(here.getX() - eps * mx, here.getY() - eps * my));
        if (narrow_in || !wide_in) {
          continue;
        }
        corners.emplace_back(here, std::atan2(-my, -mx));
      }
    }
  }
  return corners;
}


std::vector<F2CPath> TurningBase::alternativeTurns(const F2CRobot&,
    const F2CPoint&, double, const F2CPoint&, double) {
  return {};
}

F2CPath TurningBase::createTurn(const F2CRobot& robot,
    const F2CPoint& start_pos, double start_angle,
    const F2CPoint& end_pos, double end_angle, TurnReport* report) {
  const F2CPath plain =
      plainTurn(robot, start_pos, start_angle, end_pos, end_angle);
  if (this->free_space_.isEmpty()) {
    if (report != nullptr) { *report = TurnReport(); }
    return plain;
  }

  const double half_swath = 0.5 * this->swath_width_;
  const double strip_length = robot.getMinTurningRadius();
  const double in_ang = start_angle + M_PI;
  // Does it fit? -- the cheap question: stop at the first metre off the ground.
  auto fits = [&](const F2CPath& p) {
    double out = 0.0, deep = 0.0, in_swath = 0.0;
    measureOutside(p, this->free_space_, this->free_area_,
          this->discretization * 10.0,
        half_swath, strip_length, start_pos, in_ang, end_pos, end_angle,
        &out, &deep, &in_swath, 0.05);
    return out <= 0.05;
  };
  // The second, softer question: does it also keep to the ground the route
  // travels through? Ground outside it is still clear of the crop.
  auto measurePreferred = [&](const F2CPath& p, TurnReport* r) {
    if (this->preferred_space_.isEmpty()) {
      r->in_preferred = true;
      r->length_off_preferred = 0.0;
      return;
    }
    double deep = 0.0, in_swath = 0.0;
    measureOutside(p, this->preferred_space_, this->preferred_area_,
          this->discretization * 10.0,
        half_swath, strip_length, start_pos, in_ang, end_pos, end_angle,
        &r->length_off_preferred, &deep, &in_swath);
    r->in_preferred = r->length_off_preferred <= 0.05;
  };
  // The cheap form of the same question, for ranking candidates.
  auto fitsPreferred = [&](const F2CPath& p) {
    if (this->preferred_space_.isEmpty()) { return true; }
    double out = 0.0, deep = 0.0, in_swath = 0.0;
    measureOutside(p, this->preferred_space_, this->preferred_area_,
          this->discretization * 10.0,
        half_swath, strip_length, start_pos, in_ang, end_pos, end_angle,
        &out, &deep, &in_swath, 0.05);
    return out <= 0.05;
  };
  auto measure = [&](const F2CPath& p, TurnReport* r) {
    measureOutside(p, this->free_space_, this->free_area_,
          this->discretization * 10.0,
        half_swath, strip_length, start_pos, in_ang, end_pos, end_angle,
        &r->length_outside, &r->deepest_outside, &r->length_in_swath);
    r->inside = r->length_outside <= 0.05;
    measurePreferred(p, r);
  };

  // The shortest this planner can drive is the one it just gave; if that
  // stays on the ground, nothing else it offers can be both shorter and
  // inside, so do not ask.
  TurnReport best_rep;
  F2CPath best;
  // Two answers are tracked: the shortest that stays off the crop, and the
  // shortest that also keeps to the ground the route travels through.
  TurnReport pref_rep;
  F2CPath pref;
  const double radius = robot.getMinTurningRadius();
  // What the preference is allowed to cost. Staying in the corridor is worth
  // a detour, but not one longer than driving a full circle -- past that the
  // shorter turn that is merely clear of the crop is the better answer.
  const double preferred_budget = 2.0 * M_PI * radius;
  measure(plain, &best_rep);
  if (best_rep.inside) {
    best = plain;
    if (best_rep.in_preferred) {
      if (report != nullptr) { *report = best_rep; }
      return plain;
    }
  }

  // It does not. Whatever else this planner can drive between the two poses
  // is longer, but length is not what is wrong with the answer.
  TurnReport shallow_rep = best_rep;
  F2CPath shallow = plain;
  for (const F2CPath& p :
       alternativeTurns(robot, start_pos, start_angle, end_pos, end_angle)) {
    if (p.size() < 2) { continue; }
    TurnReport rep;
    measure(p, &rep);
    if (rep.inside) {
      if (best.size() == 0 || p.length() < best.length()) {
        best = p; best_rep = rep;
      }
      if (rep.in_preferred && (pref.size() == 0 || p.length() < pref.length())) {
        pref = p; pref_rep = rep;
      }
    } else if (rep.deepest_outside < shallow_rep.deepest_outside) {
      shallow = p; shallow_rep = rep;
    }
  }
  if (pref.size() > 1) {
    if (report != nullptr) { *report = pref_rep; }
    return pref;
  }

  // Nothing this planner drives stays on the ground it was given. The turn
  // leaves it where the ground turns inward, so retry through a waypoint set
  // a little inside each concave corner and headed along it. How far inside
  // decides how long the detour is and the best distance differs per corner,
  // so try a ladder and keep the shortest that fits.
  // Only corners the turn could plausibly go through. Going start -> way ->
  // end costs at least the two straight legs, so a corner outside the ellipse
  // whose foci are the poses and whose extra length is one full circle cannot
  // give a turn worth driving. Two circles, not one: measured on the bench,
  // one circle prunes waypoints that were being used and leaves 18 turns on
  // the crop instead of 10; two matches the unbounded answer exactly.
  // Without any bound the ladder walks every
  // concave corner of the free space -- on a real headland ring that is tens
  // to hundreds of corners, five offsets and two sides each, every candidate
  // measured against the whole field.
  const double reach = start_pos.distance(end_pos) + 4.0 * M_PI * radius;
  // Plain numbers, not paths: F2CPath is expensive to keep by the hundred,
  // and the winner is cheap to build again once it is known which it is.
  struct Candidate {
    double x, y, angle, length;
  };
  std::vector<Candidate> candidates;
  for (const auto& corner : this->free_corners_) {
    if (start_pos.distance(corner.first) +
        corner.first.distance(end_pos) > reach) {
      continue;
    }
    for (const double frac : {1.0, 0.5, 0.25, 0.1, 0.05}) {
      const double offset = this->waypoint_offset_ * frac * radius;
      const F2CPoint way {
          corner.first.getX() + offset * std::cos(corner.second),
          corner.first.getY() + offset * std::sin(corner.second)};
      if (!this->free_space_.isPointIn(way)) { continue; }
      for (int side = 0; side < 2; ++side) {
        const double way_angle = corner.second + (side ? -M_PI_2 : M_PI_2);
        const F2CPath first =
            plainTurn(robot, start_pos, start_angle, way, way_angle);
        if (first.size() < 2) { continue; }
        const F2CPath second =
            plainTurn(robot, way, way_angle, end_pos, end_angle);
        if (second.size() < 2) { continue; }
        F2CPath joined = first;
        joined += second;
        if (minRadius(joined) < 0.95 * radius) { continue; }
        candidates.push_back(
            {way.getX(), way.getY(), way_angle, joined.length()});
      }
    }
  }
  // Shortest first, and stop at the first that fits: that one is the shortest
  // that fits. Measuring every candidate to then pick the shortest asks the
  // expensive question of candidates already known to be worse.
  std::sort(candidates.begin(), candidates.end(),
      [](const Candidate& a, const Candidate& b) {
        return a.length < b.length;
      });
  for (const Candidate& c : candidates) {
    // Sorted shortest first, so once a candidate costs more than the answer
    // in hand plus what the preference is worth, neither it nor anything
    // after it can be chosen.
    if (best.size() > 1 && c.length > best.length() + preferred_budget) {
      break;
    }
    const F2CPoint way {c.x, c.y};
    F2CPath joined = plainTurn(robot, start_pos, start_angle, way, c.angle);
    joined += plainTurn(robot, way, c.angle, end_pos, end_angle);
    if (!fits(joined)) { continue; }
    // Measure in full only what is actually taken: the cheap questions above
    // settle the ranking, and most candidates lose it.
    const bool keeps_to_route = fitsPreferred(joined);
    if (best.size() == 0 || joined.length() < best.length()) {
      best = joined;
      measure(joined, &best_rep);
      best_rep.used_waypoint = true;
      best_rep.waypoint = way;
    }
    if (keeps_to_route) {
      // Sorted shortest first, so this is the shortest that keeps to the
      // route's ground; nothing further along can beat it.
      pref = joined;
      measure(joined, &pref_rep);
      pref_rep.used_waypoint = true;
      pref_rep.waypoint = way;
      break;
    }
  }

  // The answer, in order: off the crop and on the route's ground, unless
  // that costs more than driving a full circle; then merely off the crop;
  // then, if nothing fits, the one that goes least deep into it.
  if (pref.size() > 1 &&
      (best.size() < 2 || pref.length() <= best.length() + preferred_budget)) {
    if (report != nullptr) { *report = pref_rep; }
    return pref;
  }
  if (best.size() > 1) {
    if (report != nullptr) { *report = best_rep; }
    return best;
  }
  if (report != nullptr) { *report = shallow_rep; }
  return shallow.size() > 1 ? shallow : plain;
}

const F2CCells& TurningBase::getFreeSpace() const {
  return this->free_space_;
}

const F2CCells& TurningBase::getPreferredSpace() const {
  return this->preferred_space_;
}

void TurningBase::setPreferredSpace(const F2CCells& preferred) {
  this->preferred_space_ = preferred;
  this->preferred_area_ = f2c::PreparedArea(preferred);
}

void TurningBase::setFreeSpace(const F2CCells& free_space) {
  this->free_space_ = free_space;
  // Read once here rather than on every sample of every candidate turn.
  this->free_area_ = f2c::PreparedArea(free_space);
  // The corners do not move while the ground does not, and createTurn asked
  // for them again on every single turn.
  this->free_corners_ = concaveCorners(free_space);
}

double TurningBase::getSwathWidth() const {
  return this->swath_width_;
}

void TurningBase::setSwathWidth(double width) {
  this->swath_width_ = std::fabs(width);
}

double TurningBase::getWaypointOffset() const {
  return this->waypoint_offset_;
}

void TurningBase::setWaypointOffset(double offset) {
  this->waypoint_offset_ = std::fabs(offset);
}

F2CPath TurningBase::createTurnIfNotCached(const F2CRobot& robot,
    double dist_start_end, double start_angle, double end_angle) {
  std::vector<int> v_turn {
        static_cast<int>(1e3 * robot.getMaxCurv()),
        static_cast<int>(1e3 * robot.getMaxDiffCurv()),
        static_cast<int>(1e3 * dist_start_end),
        static_cast<int>(1e3 * start_angle),
        static_cast<int>(1e3 * end_angle)
  };
  auto it = path_cache_.find(v_turn);
  if (it != path_cache_.end()) {
    return it->second;
  }
  auto path = createSimpleTurn(robot, dist_start_end, start_angle, end_angle);
  path_cache_.insert({v_turn, path});
  return path;
}


bool TurningBase::isTurnValid(const types::Path& path,
    double dist_start_end, double end_angle,
    double max_dist, double max_rot_error) {
  for (auto&& s : path) {
    if (s.point.getY() < -max_dist) {
      return false;
    }
  }
  F2CPoint p_end = path.atEnd();
  return (cos(path.back().angle - end_angle) >= 1 - max_rot_error) &&
    (fabs(p_end.getX() - dist_start_end) < max_dist) &&
    (fabs(p_end.getY()) < max_dist);
}

bool TurningBase::hasContinuousCurvature() const {
  return false;
}

double TurningBase::getDiscretization() const {
  return this->discretization;
}

void TurningBase::setDiscretization(double d) {
  this->discretization = d;
}

bool TurningBase::getUsingCache() const {
  return this->using_cache;
}

void TurningBase::setUsingCache(bool c) {
  this->using_cache = c;
}


}  // namespace f2c::pp

