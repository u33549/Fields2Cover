//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <set>
#include <thread>
#include <utility>
#include <vector>
#include "fields2cover/route_planning/free_space_route_planner.h"

namespace f2c::rp {

namespace {

// A ring flattened into plain arrays, with its box. Asking a polygon whether it
// holds a point is the inner loop of a quadratic graph, so it is worth keeping
// out of the geometry library.
struct Ring {
  std::vector<double> x, y;
  double x0 {1e18}, x1 {-1e18}, y0 {1e18}, y1 {-1e18};

  void add(double px, double py) {
    x.push_back(px);
    y.push_back(py);
    x0 = std::min(x0, px); x1 = std::max(x1, px);
    y0 = std::min(y0, py); y1 = std::max(y1, py);
  }

  bool holds(double px, double py) const {
    if (px < x0 || px > x1 || py < y0 || py > y1) {
      return false;
    }
    bool in = false;
    const size_t n = x.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
      if (((y[i] > py) != (y[j] > py)) &&
          (px < (x[j] - x[i]) * (py - y[i]) / (y[j] - y[i]) + x[i])) {
        in = !in;
      }
    }
    return in;
  }
};

struct Poly {
  Ring outer;
  std::vector<Ring> holes;

  bool holds(double px, double py) const {
    if (!outer.holds(px, py)) {
      return false;
    }
    for (const auto& h : holes) {
      if (h.holds(px, py)) {
        return false;
      }
    }
    return true;
  }
};

class Area {
 public:
  explicit Area(const F2CCells& cs) {
    for (size_t i = 0; i < cs.size(); ++i) {
      const F2CCell c = cs.getGeometry(i);
      Poly p;
      readRing(c.getExteriorRing(), &p.outer);
      for (size_t r = 0; r + 1 < c.size(); ++r) {
        Ring h;
        readRing(c.getInteriorRing(r), &h);
        p.holes.push_back(std::move(h));
      }
      polys_.push_back(std::move(p));
    }
  }

  bool holds(double px, double py) const {
    for (const auto& p : polys_) {
      if (p.holds(px, py)) {
        return true;
      }
    }
    return false;
  }

  bool isEmpty() const { return polys_.empty(); }

  // Every ring edge, for the crossing test.
  void edges(std::vector<std::array<double, 4>>* out) const {
    for (const auto& p : polys_) {
      pushRing(p.outer, out);
      for (const auto& h : p.holes) {
        pushRing(h, out);
      }
    }
  }

 private:
  static void readRing(const F2CLinearRing& r, Ring* out) {
    for (size_t k = 0; k < r.size(); ++k) {
      const F2CPoint v = r.getGeometry(k);
      out->add(v.getX(), v.getY());
    }
  }
  static void pushRing(const Ring& r, std::vector<std::array<double, 4>>* out) {
    for (size_t i = 0; i + 1 < r.x.size(); ++i) {
      out->push_back({r.x[i], r.y[i], r.x[i + 1], r.y[i + 1]});
    }
  }
  std::vector<Poly> polys_;
};

double side(double ox, double oy, double ax, double ay, double bx, double by) {
  return (ax - ox) * (by - oy) - (ay - oy) * (bx - ox);
}

int sign(double v, double eps) { return v > eps ? 1 : (v < -eps ? -1 : 0); }

// Whether two segments cross properly. Meeting at an end does not count: the
// nodes sit on the border, so every edge touches it.
bool crosses(const std::array<double, 4>& s, double px, double py,
    double qx, double qy, double eps) {
  const int d1 = sign(side(s[0], s[1], s[2], s[3], px, py), eps);
  const int d2 = sign(side(s[0], s[1], s[2], s[3], qx, qy), eps);
  const int d3 = sign(side(px, py, qx, qy, s[0], s[1]), eps);
  const int d4 = sign(side(px, py, qx, qy, s[2], s[3]), eps);
  return d1 * d2 < 0 && d3 * d4 < 0;
}

void ringCorners(const F2CCells& cs, std::vector<F2CPoint>* out) {
  for (size_t i = 0; i < cs.size(); ++i) {
    const F2CCell c = cs.getGeometry(i);
    for (size_t r = 0; r < c.size(); ++r) {
      const F2CLinearRing ring =
          (r == 0) ? c.getExteriorRing() : c.getInteriorRing(r - 1);
      for (size_t k = 0; k + 1 < ring.size(); ++k) {
        out->push_back(ring.getGeometry(k));
      }
    }
  }
}

// Union-find, to count the pieces the graph falls into.
class Pieces {
 public:
  explicit Pieces(size_t n) : parent_(n) {
    std::iota(parent_.begin(), parent_.end(), 0);
  }
  size_t root(size_t a) {
    return parent_[a] == a ? a : parent_[a] = root(parent_[a]);
  }
  void join(size_t a, size_t b) { parent_[root(a)] = root(b); }
  size_t count() {
    std::set<size_t> roots;
    for (size_t i = 0; i < parent_.size(); ++i) {
      roots.insert(root(i));
    }
    return roots.size();
  }

 private:
  std::vector<size_t> parent_;
};

}  // namespace


F2CGraph2D FreeSpaceRoutePlanner::createShortestGraph(
    const F2CCells& cells, const F2CSwathsByCells& swaths_by_cells,
    double d_tol) const {
  F2CGraph2D g;
  this->ground_ = cells;   // what the connections built on this graph align to
  if (cells.size() == 0) {
    n_components_ = 0;
    return g;
  }

  const Area ground(cells);
  std::vector<std::array<double, 4>> border;
  ground.edges(&border);

  // The band along the edge the route should stay out of where it can. Priced,
  // not forbidden: a swath end sits on the border and has no other way out.
  const bool has_band = (clearance_ > 0.0 && clearance_cost_ > 0.0);
  const F2CCells inner =
      has_band ? cells.buffer(-clearance_) : F2CCells();
  const Area room(inner);

  // Nodes: every swath end, and every corner of the ground.
  std::vector<F2CPoint> nodes;
  const F2CCells outline =
      (corner_tol_ > 0.0) ? cells.simplify(corner_tol_) : cells;
  ringCorners(outline, &nodes);
  const size_t n_corners = nodes.size();
  // A swath end sits on the border of the ground, and the machine standing
  // there is still on the swath's own line. Joined straight to whatever corner
  // is nearest, the first metres of the connection cut across that line at an
  // angle no machine leaving the swath can hold. Walk in along the swath's axis
  // first and let the route start from there.
  const double entry = this->getTurnRoom();
  std::vector<std::pair<size_t, F2CPoint>> entries;   // swath end -> its entry
  for (auto&& swaths : swaths_by_cells) {
    for (auto&& s : swaths) {
      const F2CPoint ends[2] = {s.startPoint(), s.endPoint()};
      const double aways[2] = {s.getInAngle() + M_PI, s.getOutAngle()};
      for (int e = 0; e < 2; ++e) {
        nodes.push_back(ends[e]);
        if (entry <= 0.0) {
          continue;
        }
        // A swath end stands off the ground by half the strip its own swath
        // covers: that is the band the ground was eroded by, and the machine
        // crossing it is still on what it has just cut. An end further off
        // than that is not standing at the edge of its swath -- it is inside
        // ground the machine cannot enter, where the ground has pinched shut,
        // and its axis finds ground again only on the far side of that.
        // Standing exactly that far off is the defining case, not a borderline
        // one, so the comparison is made to the tolerance the caller gave.
        if (ends[e].distance(cells) >
            0.5 * s.getWidth() + std::max(d_tol, 1e-9)) {
          continue;
        }
        // Along the axis until the ground holds it, and then the room asks
        // for on top of that: the end itself may stand outside the ground --
        // half the machine's width from the crop is where a swath ends -- so
        // the room has to be measured from where the ground starts, not from
        // the end. The leg the turn planner is given starts at that point.
        // How far away the ground may be is a separate limit from the room
        // measured once it is found: a swath whose end stands deep inside the
        // crop -- where the ground pinches shut -- has ground along its axis
        // only on the far side of that crop, and an entry set there is only
        // reached by driving through it.
        const double step = (sample_step_ > 0.0) ? sample_step_ : 0.5;
        const double reach = 2.0 * entry;
        const double dx = std::cos(aways[e]), dy = std::sin(aways[e]);
        const Area area(cells);
        F2CPoint found = ends[e];
        bool on = false;
        double reached = 0.0;
        for (double t = step; t <= reach + entry + step; t += step) {
          const F2CPoint q(ends[e].getX() + t * dx, ends[e].getY() + t * dy);
          if (!area.holds(q.getX(), q.getY())) {
            if (on || t > reach) { break; }
            continue;
          }
          if (!on) { reached = t; }
          found = q;
          on = true;
          if (t - reached >= entry) { break; }
        }
        if (on && found.distance(ends[e]) > 1e-9) {
          entries.push_back({nodes.size() - 1, found});
        }
      }
    }
  }
  const size_t n_before_entries = nodes.size();
  for (const auto& e : entries) {
    nodes.push_back(e.second);
  }
  // An end that has an entry leaves along it and nowhere else: a chord straight
  // off the end is a sideways start no machine can hold. An end without one --
  // no ground along its axis -- keeps its edges, so no swath is cut off.
  std::vector<bool> leaves_by_entry(nodes.size(), false);
  for (const auto& e : entries) {
    leaves_by_entry[e.first] = true;
  }
  if (nodes.size() < 2) {
    n_components_ = nodes.size();
    return g;
  }

  const double step = (sample_step_ > 0.0) ? sample_step_ : 0.5;
  const double eps = std::max(d_tol, 1e-9);

  // Whether the segment between two nodes stays on the ground, and how much of
  // it runs inside the clearance band.
  auto reaches = [&](const F2CPoint& a, const F2CPoint& b,
      bool a_end, bool b_end,
      double* in_band) -> bool {
    const double ax = a.getX(), ay = a.getY();
    const double bx = b.getX(), by = b.getY();
    const double dx = bx - ax, dy = by - ay;
    const double len2 = dx * dx + dy * dy;
    if (len2 < 1e-12) {
      return false;
    }
    for (const auto& s : border) {
      if (crosses(s, ax, ay, bx, by, eps)) {
        return false;
      }
    }
    // Crossing alone is not enough: a chord between two corners can clip a spur
    // by entering and leaving at corners, which is not a proper crossing.
    const double len = std::sqrt(len2);
    const int n = std::max(2, static_cast<int>(std::ceil(len / step)));
    const double dl = len / n;
    const double free_of_band = 1.5 * clearance_;
    double band = 0.0;
    for (int t = 0; t < n; ++t) {
      const double u = (t + 0.5) / n;
      const double mx = ax + dx * u, my = ay + dy * u;
      if (!ground.holds(mx, my)) {
        return false;
      }
      if (has_band && !room.isEmpty() && !room.holds(mx, my)) {
        const double from_a = a_end ? std::hypot(mx - ax, my - ay) : 1e18;
        const double from_b = b_end ? std::hypot(mx - bx, my - by) : 1e18;
        if (std::min(from_a, from_b) > free_of_band) {
          band += dl;
        }
      }
    }
    if (in_band != nullptr) {
      *in_band = band;
    }
    return true;
  };

  // The pairs are independent, so share them out; F2CGraph2D is filled after.
  struct Edge { size_t i, k; double band; };
  unsigned n_threads = std::thread::hardware_concurrency();
  n_threads = std::max(1u, std::min(n_threads, 16u));
  std::vector<std::vector<Edge>> found(n_threads);
  std::vector<std::thread> workers;
  for (unsigned t = 0; t < n_threads; ++t) {
    workers.emplace_back([&, t] {
      for (size_t i = t; i < nodes.size(); i += n_threads) {
        for (size_t k = i + 1; k < nodes.size(); ++k) {
          if (leaves_by_entry[i] || leaves_by_entry[k]) {
            continue;
          }
          double band = 0.0;
          if (reaches(nodes[i], nodes[k], i >= n_corners, k >= n_corners,
                &band)) {
            found[t].push_back({i, k, band});
          }
        }
      }
    });
  }
  for (auto& w : workers) {
    w.join();
  }

  Pieces pieces(nodes.size());
  for (const auto& v : found) {
    for (const auto& e : v) {
      const double len = nodes[e.i].distance(nodes[e.k]);
      const double cost = len + clearance_cost_ * e.band;
      g.addEdge(nodes[e.i], nodes[e.k], static_cast<int64_t>(1e3 * cost));
      pieces.join(e.i, e.k);
    }
  }
  // The step in from a swath end is the machine's own line, not a shortcut the
  // graph has to justify -- and the end stands on the border, where a sampled
  // segment is not held. Join the two unconditionally, or the entry is a node
  // nothing reaches and the end keeps leaving sideways.
  for (size_t e = 0; e < entries.size(); ++e) {
    const size_t a = entries[e].first, b = n_before_entries + e;
    const double len = nodes[a].distance(nodes[b]);
    g.addEdge(nodes[a], nodes[b], static_cast<int64_t>(1e3 * len));
    pieces.join(a, b);
  }
  n_components_ = pieces.count();
  return g;
}

void FreeSpaceRoutePlanner::setClearance(double clearance) {
  this->clearance_ = std::fabs(clearance);
}

double FreeSpaceRoutePlanner::getClearance() const {
  return this->clearance_;
}

void FreeSpaceRoutePlanner::setClearanceCost(double cost) {
  this->clearance_cost_ = std::fabs(cost);
}

double FreeSpaceRoutePlanner::getClearanceCost() const {
  return this->clearance_cost_;
}

void FreeSpaceRoutePlanner::setCornerTolerance(double tol) {
  this->corner_tol_ = std::fabs(tol);
}

double FreeSpaceRoutePlanner::getCornerTolerance() const {
  return this->corner_tol_;
}

void FreeSpaceRoutePlanner::setSampleStep(double step) {
  this->sample_step_ = std::fabs(step);
}

double FreeSpaceRoutePlanner::getSampleStep() const {
  return this->sample_step_;
}

void FreeSpaceRoutePlanner::setTurnRoom(double room) {
  this->turn_room_ = std::fabs(room);
}

double FreeSpaceRoutePlanner::getTurnRoom() const {
  return (this->turn_room_ >= 0.0) ? this->turn_room_ : this->clearance_;
}

size_t FreeSpaceRoutePlanner::getComponentCount() const {
  return this->n_components_;
}

namespace {

// A leg, as the line it runs on.
struct Line {
  F2CPoint p;
  double ux, uy;
};

bool meet(const Line& a, const Line& b, F2CPoint* out) {
  const double c = a.ux * b.uy - a.uy * b.ux;
  if (std::fabs(c) < 1e-3) {
    return false;              // parallel
  }
  const double dx = b.p.getX() - a.p.getX(), dy = b.p.getY() - a.p.getY();
  const double t = (dx * b.uy - dy * b.ux) / c;
  *out = F2CPoint(a.p.getX() + t * a.ux, a.p.getY() + t * a.uy);
  return true;
}

// How far the ground reaches from a point on its border, along `n`.
double roomFrom(const Area& ground, const F2CPoint& foot, double nx, double ny,
    double reach) {
  auto in = [&](double t) {
    return ground.holds(foot.getX() + t * nx, foot.getY() + t * ny);
  };
  double inside = 0.05, outside = -1.0;
  for (double t = 0.05; t < reach; t += 0.25) {
    if (!in(t)) {
      outside = t;
      break;
    }
    inside = t;
  }
  if (outside < 0.0) {
    return inside;
  }
  for (int k = 0; k < 14; ++k) {
    const double m = 0.5 * (inside + outside);
    (in(m) ? inside : outside) = m;
  }
  return inside;
}

bool onGround(const Area& ground, const F2CPoint& a, const F2CPoint& b,
    double step) {
  const double len = a.distance(b);
  if (len < 1e-9) {
    return true;
  }
  const int n = std::max(2, static_cast<int>(std::ceil(len / step)));
  for (int t = 0; t < n; ++t) {
    const double u = (t + 0.5) / n;
    if (!ground.holds(a.getX() + (b.getX() - a.getX()) * u,
                      a.getY() + (b.getY() - a.getY()) * u)) {
      return false;
    }
  }
  return true;
}

// The border a leg runs along: nearly parallel, nearest, and facing it.
struct Border {
  bool found {false};
  F2CPoint foot;
  double ux {0}, uy {0};     // along the border, pointing the way the leg goes
  double nx {0}, ny {0};     // into the ground
  double room {0}, at {0};   // ground left across it, and where the leg is now
};

Border borderOf(const std::vector<std::array<double, 4>>& edges,
    const Area& ground, const F2CPoint& a, const F2CPoint& b) {
  Border w;
  const double len = a.distance(b);
  if (len < 1e-9) {
    return w;
  }
  const double ux = (b.getX() - a.getX()) / len, uy = (b.getY() - a.getY()) / len;
  const F2CPoint mid(0.5 * (a.getX() + b.getX()), 0.5 * (a.getY() + b.getY()));
  const double parallel = std::sin(10.0 * M_PI / 180.0);
  double best = 1e18;
  size_t pick = edges.size();
  for (size_t e = 0; e < edges.size(); ++e) {
    const double ex = edges[e][2] - edges[e][0], ey = edges[e][3] - edges[e][1];
    const double el = std::hypot(ex, ey);
    if (el < 1e-6) {
      continue;
    }
    const double dx = ex / el, dy = ey / el;
    if (std::fabs(ux * dy - uy * dx) > parallel) {
      continue;
    }
    const double along =
        (mid.getX() - edges[e][0]) * dx + (mid.getY() - edges[e][1]) * dy;
    if (along < -1.0 || along > el + 1.0) {
      continue;
    }
    const double off = std::fabs(
        (mid.getX() - edges[e][0]) * dy - (mid.getY() - edges[e][1]) * dx);
    if (off < best) {
      best = off;
      pick = e;
    }
  }
  if (pick == edges.size()) {
    return w;
  }
  const double ex = edges[pick][2] - edges[pick][0];
  const double ey = edges[pick][3] - edges[pick][1];
  const double el = std::hypot(ex, ey);
  double dx = ex / el, dy = ey / el;
  if (dx * ux + dy * uy < 0.0) {
    dx = -dx; dy = -dy;
  }
  const double along =
      (mid.getX() - edges[pick][0]) * dx + (mid.getY() - edges[pick][1]) * dy;
  const F2CPoint foot(edges[pick][0] + along * dx, edges[pick][1] + along * dy);
  double nx = -dy, ny = dx;
  if (!ground.holds(foot.getX() + 0.05 * nx, foot.getY() + 0.05 * ny)) {
    nx = -nx; ny = -ny;
    if (!ground.holds(foot.getX() + 0.05 * nx, foot.getY() + 0.05 * ny)) {
      return w;
    }
  }
  w.found = true;
  w.foot = foot;
  w.ux = dx; w.uy = dy;
  w.nx = nx; w.ny = ny;
  w.at = best;
  w.room = roomFrom(ground, foot, nx, ny, 4.0 * best + 200.0);
  return w;
}

// `pinned_next`: the point next to each end is the step in along the swath's
// own axis. It is where the machine actually leaves the swath, so it is as
// fixed as the end itself -- aligning it to a border would put the start of
// the connection back across the swath's line.
F2CMultiPoint alignConnection(const F2CMultiPoint& mp, const F2CCells& cells,
    double margin, double step, bool pinned_next) {
  const size_t n = mp.size();
  if (margin <= 0.0 || n < 3 || cells.isEmpty()) {
    return mp;
  }
  std::vector<F2CPoint> p;
  for (size_t i = 0; i < n; ++i) {
    p.push_back(mp.getGeometry(i));
  }
  const Area ground(cells);
  std::vector<std::array<double, 4>> edges;
  ground.edges(&edges);

  // The ends are swath ends and stay put, so the legs that reach them keep
  // their line; only what runs between borders is laid along one.
  std::vector<Line> legs(n - 1);
  std::vector<Border> borders(n - 1);
  for (size_t i = 0; i + 1 < n; ++i) {
    const double len = p[i].distance(p[i + 1]);
    if (len < 1e-9) {
      return mp;
    }
    legs[i] = Line{p[i], (p[i + 1].getX() - p[i].getX()) / len,
                         (p[i + 1].getY() - p[i].getY()) / len};
    if (i == 0 || i + 2 == n) {
      continue;                // reaches a swath end: its line cannot move
    }
    if (pinned_next && (i == 1 || i + 3 == n)) {
      continue;                // reaches the step in from one, same thing
    }
    borders[i] = borderOf(edges, ground, p[i], p[i + 1]);
  }

  // What each corner asks for, measured on the legs as they will run.
  auto heading = [&](size_t i) {
    return borders[i].found ? std::atan2(borders[i].uy, borders[i].ux)
                            : std::atan2(legs[i].uy, legs[i].ux);
  };
  auto asks = [&](size_t i) {
    if (i == 0 || i + 1 >= n) {
      return 0.0;
    }
    const double turn = F2CPoint::getAngleDiffAbs(heading(i - 1), heading(i));
    return margin * (1.0 - std::cos(0.5 * turn));
  };
  // Whether a leg's border is the one its corner turns towards: the arc only
  // cuts the inside of the turn, so that is the side to keep away from.
  auto inside = [&](size_t leg, size_t corner) {
    if (corner == 0 || corner + 1 >= n) {
      return false;
    }
    const F2CPoint a = p[corner - 1] - p[corner], b = p[corner + 1] - p[corner];
    const double la = std::hypot(a.getX(), a.getY());
    const double lb = std::hypot(b.getX(), b.getY());
    if (la < 1e-9 || lb < 1e-9) {
      return false;
    }
    const double bx = a.getX() / la + b.getX() / lb;
    const double by = a.getY() / la + b.getY() / lb;
    return bx * -borders[leg].nx + by * -borders[leg].ny > 0.0;
  };

  for (size_t i = 0; i + 1 < n; ++i) {
    const Border& w = borders[i];
    if (!w.found) {
      continue;
    }
    double want = 0.0;
    if (inside(i, i)) {
      want = std::max(want, asks(i));
    }
    if (inside(i, i + 1)) {
      want = std::max(want, asks(i + 1));
    }
    // Never towards the border: the leg is only laid parallel to it, and moved
    // out when what it carries needs more room than it has.
    const double off = std::min(std::max(w.at, want), w.room);
    legs[i] = Line{F2CPoint(w.foot.getX() + off * w.nx, w.foot.getY() + off * w.ny),
                   w.ux, w.uy};
  }

  std::vector<F2CPoint> out = p;
  for (size_t i = 1; i + 1 < n; ++i) {
    if (pinned_next && (i == 1 || i + 2 == n)) {
      continue;                // the step in from a swath end stays put
    }
    F2CPoint x;
    if (meet(legs[i - 1], legs[i], &x)) {
      out[i] = x;
    }
  }

  // A corner both of whose legs reach a swath end -- a connection that is one
  // corner, the common case -- cannot be given room by laying a leg along a
  // border, since neither leg may move. Push the corner itself away from the
  // side the turn cuts towards, by what the arc cuts off there.
  for (size_t i = 1; i + 1 < n; ++i) {
    if (pinned_next && (i == 1 || i + 2 == n)) {
      continue;
    }
    if (out[i].distance(p[i]) > 1e-9) {
      continue;
    }
    const double turn = F2CPoint::getAngleDiffAbs(heading(i - 1), heading(i));
    const double half = std::cos(0.5 * turn);
    if (turn < 0.05 || half < 1.0 / 3.0) {
      continue;                // straight, or too sharp to round at all
    }
    const F2CPoint a = p[i - 1] - p[i], b = p[i + 1] - p[i];
    const double la = std::hypot(a.getX(), a.getY());
    const double lb = std::hypot(b.getX(), b.getY());
    if (la < 1e-9 || lb < 1e-9) {
      continue;
    }
    const double bx = a.getX() / la + b.getX() / lb;
    const double by = a.getY() / la + b.getY() / lb;
    const double bl = std::hypot(bx, by);
    if (bl < 1e-6) {
      continue;
    }
    const double want = margin * (1.0 / half - 1.0);
    for (const double s : {1.0, 0.75, 0.5, 0.25}) {
      const F2CPoint c(p[i].getX() - s * want * bx / bl,
                       p[i].getY() - s * want * by / bl);
      if (ground.holds(c.getX(), c.getY()) &&
          onGround(ground, out[i - 1], c, step) &&
          onGround(ground, c, out[i + 1], step)) {
        out[i] = c;
        break;
      }
    }
  }
  for (size_t i = 1; i + 1 < n; ++i) {
    if (out[i].distance(p[i]) > 4.0 * margin || !ground.holds(out[i].getX(), out[i].getY())) {
      return mp;
    }
  }
  for (size_t i = 0; i + 1 < n; ++i) {
    if (!onGround(ground, out[i], out[i + 1], step)) {
      return mp;
    }
  }
  F2CMultiPoint res;
  for (const auto& q : out) {
    res.addPoint(q);
  }
  return res;
}

}  // namespace

F2CRoute FreeSpaceRoutePlanner::transformSolutionToRoute(
    const std::vector<long long int>& route_ids,
    const F2CSwathsByCells& swaths_by_cells,
    const F2CGraph2D& coverage_graph,
    F2CGraph2D& shortest_graph) const {
  F2CRoute route = RoutePlannerBase::transformSolutionToRoute(
      route_ids, swaths_by_cells, coverage_graph, shortest_graph);
  const double room = this->getTurnRoom();
  if (room <= 0.0 || this->ground_.isEmpty()) {
    return route;
  }
  const double step = (this->sample_step_ > 0.0) ? this->sample_step_ : 0.5;
  for (size_t i = 0; i < route.sizeConnections(); ++i) {
    route.getConnection(i) = alignConnection(
        route.getConnection(i), this->ground_, room, step, room > 0.0);
  }
  return route;
}

}  // namespace f2c::rp
