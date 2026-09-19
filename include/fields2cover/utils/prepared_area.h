//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#pragma once
#ifndef FIELDS2COVER_UTILS_PREPARED_AREA_H_
#define FIELDS2COVER_UTILS_PREPARED_AREA_H_

#include <array>
#include <vector>
#include "fields2cover/types.h"

namespace f2c {

/// Cells read once so that asking "is this point on them" is cheap.
///
/// `Cells::isPointIn` goes through OGR to GEOS, which rebuilds both geometries
/// for every call: fine now and then, ruinous in a loop that samples a path.
/// Reading the rings into plain arrays once and testing with a bounding box
/// and a ray cast answers the same question about 60 times faster on a
/// 360-corner field.
///
/// Two differences from `Cells::isPointIn` matter. A point exactly on a ring
/// is **undefined** here, where GEOS reports it as outside; callers that care
/// about the edge have to say so themselves. And the test is plain floating
/// point, not GEOS's robust predicates, so a point within rounding distance of
/// an edge can go either way.
class PreparedArea {
 public:
  PreparedArea() = default;

  /// Read the rings of \a cells. Holes are read too and count as not held.
  explicit PreparedArea(const F2CCells& cells);

  /// Whether the point is on the cells: inside an outer ring and not inside
  /// one of its holes. A point on a ring itself is undefined.
  bool holds(double x, double y) const;

  /// Whether no cell was read.
  bool isEmpty() const;

  /// Append every ring edge as {x1, y1, x2, y2}, for segment tests.
  void edges(std::vector<std::array<double, 4>>* out) const;

 private:
  struct Ring {
    std::vector<double> x, y;
    double x0 {1e18}, x1 {-1e18}, y0 {1e18}, y1 {-1e18};
    void add(double px, double py);
    bool holds(double px, double py) const;
  };
  struct Poly {
    Ring outer;
    std::vector<Ring> holes;
    bool holds(double px, double py) const;
  };
  static void readRing(const F2CLinearRing& ring, Ring* out);
  std::vector<Poly> polys_;
};

}  // namespace f2c

#endif  // FIELDS2COVER_UTILS_PREPARED_AREA_H_
