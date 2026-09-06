//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <algorithm>
#include <limits>
#include <cmath>
#include <utility>
#include <vector>
#include "fields2cover/path_planning/constrained_turn.h"

namespace f2c::pp {

namespace {

// Length of path that runs outside the region, and how far outside it gets.
// Sampled rather than intersected: the answer feeds a comparison between
// candidates, and every candidate is measured the same way.
void measureOutside(const F2CPath& path, const F2CCells& free_space,
    double step, double* length_out, double* deepest_out) {
  *length_out = 0.0;
  *deepest_out = 0.0;
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
      if (free_space.isPointIn(p)) {
        continue;
      }
      *length_out += len / n;
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

std::vector<std::pair<F2CPoint, double>> ConstrainedTurn::concaveCorners(
    const F2CCells& free_space, double eps) {
  std::vector<std::pair<F2CPoint, double>> corners;
  for (auto&& cell : free_space) {
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
        const bool narrow_in = free_space.isPointIn(
            F2CPoint(here.getX() + eps * mx, here.getY() + eps * my));
        const bool wide_in = free_space.isPointIn(
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

ConstrainedTurn::Result ConstrainedTurn::createTurn(const F2CRobot& robot,
    const F2CCells& free_space, const F2CPoint& start_pos, double start_angle,
    const F2CPoint& end_pos, double end_angle, TurningBase& turn) const {
  Result res;
  res.path = turn.createTurn(
      robot, start_pos, start_angle, end_pos, end_angle);
  measureOutside(res.path, free_space, this->sample_step_,
      &res.length_outside, &res.deepest_outside);
  if (res.length_outside <= this->tolerance_) {
    res.inside = true;
    return res;
  }

  // The plain turn leaves the free space. Retry through one waypoint per
  // concave corner, offset inward so the turning circle has room, and headed
  // along the corner rather than into it.
  const double radius = robot.getMinTurningRadius();
  const auto corners = concaveCorners(free_space);
  F2CPath best;
  double best_len = 0.0;
  for (const auto& corner : corners) {
    const F2CPoint way {
        corner.first.getX() + this->waypoint_offset_ * radius *
            std::cos(corner.second),
        corner.first.getY() + this->waypoint_offset_ * radius *
            std::sin(corner.second)};
    if (!free_space.isPointIn(way)) {
      continue;
    }
    for (int side = 0; side < 2; ++side) {
      const double way_angle = corner.second +
          (side ? -M_PI_2 : M_PI_2);
      F2CPath first = turn.createTurn(
          robot, start_pos, start_angle, way, way_angle);
      if (first.size() < 2) {
        continue;
      }
      const F2CPath second = turn.createTurn(
          robot, way, way_angle, end_pos, end_angle);
      if (second.size() < 2) {
        continue;
      }
      F2CPath joined = first;
      joined += second;
      double out_len = 0.0;
      double out_deep = 0.0;
      measureOutside(joined, free_space, this->sample_step_,
          &out_len, &out_deep);
      if (out_len > this->tolerance_) {
        continue;
      }
      if (minRadius(joined) < 0.95 * radius) {
        continue;
      }
      const double len = joined.length();
      if (best.size() == 0 || len < best_len) {
        best = joined;
        best_len = len;
      }
    }
  }
  if (best.size() > 0) {
    res.path = best;
    res.inside = true;
    res.used_waypoint = true;
    measureOutside(res.path, free_space, this->sample_step_,
        &res.length_outside, &res.deepest_outside);
  }
  return res;
}

double ConstrainedTurn::getWaypointOffset() const {
  return this->waypoint_offset_;
}

void ConstrainedTurn::setWaypointOffset(double offset) {
  this->waypoint_offset_ = std::fabs(offset);
}

double ConstrainedTurn::getSampleStep() const {
  return this->sample_step_;
}

void ConstrainedTurn::setSampleStep(double step) {
  this->sample_step_ = std::fabs(step);
}

double ConstrainedTurn::getTolerance() const {
  return this->tolerance_;
}

void ConstrainedTurn::setTolerance(double tol) {
  this->tolerance_ = std::fabs(tol);
}

}  // namespace f2c::pp
