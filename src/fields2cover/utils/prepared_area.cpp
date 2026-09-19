//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <algorithm>
#include "fields2cover/utils/prepared_area.h"

namespace f2c {

void PreparedArea::Ring::add(double px, double py) {
  x.push_back(px);
  y.push_back(py);
  x0 = std::min(x0, px); x1 = std::max(x1, px);
  y0 = std::min(y0, py); y1 = std::max(y1, py);
}

bool PreparedArea::Ring::holds(double px, double py) const {
  // The box answers most points without touching the ring at all.
  if (px < x0 || px > x1 || py < y0 || py > y1) {
    return false;
  }
  // Ray cast: count the edges a ray to the right of the point crosses.
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

bool PreparedArea::Poly::holds(double px, double py) const {
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

PreparedArea::PreparedArea(const F2CCells& cells) {
  for (size_t i = 0; i < cells.size(); ++i) {
    const F2CCell c = cells.getGeometry(i);
    Poly p;
    readRing(c.getExteriorRing(), &p.outer);
    // Cell::size() counts the exterior ring as well as the holes.
    for (size_t r = 0; r + 1 < c.size(); ++r) {
      Ring h;
      readRing(c.getInteriorRing(r), &h);
      p.holes.push_back(std::move(h));
    }
    polys_.push_back(std::move(p));
  }
}

bool PreparedArea::holds(double x, double y) const {
  for (const auto& p : polys_) {
    if (p.holds(x, y)) {
      return true;
    }
  }
  return false;
}

bool PreparedArea::isEmpty() const { return polys_.empty(); }

void PreparedArea::edges(std::vector<std::array<double, 4>>* out) const {
  const auto push = [out](const Ring& r) {
    for (size_t i = 0; i + 1 < r.x.size(); ++i) {
      out->push_back({r.x[i], r.y[i], r.x[i + 1], r.y[i + 1]});
    }
  };
  for (const auto& p : polys_) {
    push(p.outer);
    for (const auto& h : p.holes) {
      push(h);
    }
  }
}

void PreparedArea::readRing(const F2CLinearRing& r, Ring* out) {
  for (size_t k = 0; k < r.size(); ++k) {
    const F2CPoint v = r.getGeometry(k);
    out->add(v.getX(), v.getY());
  }
}

}  // namespace f2c
