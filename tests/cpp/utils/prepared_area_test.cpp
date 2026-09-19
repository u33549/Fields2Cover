//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <gtest/gtest.h>
#include <array>
#include <vector>
#include "fields2cover/types.h"
#include "fields2cover/utils/prepared_area.h"

namespace {
F2CCell square(double x0, double y0, double x1, double y1) {
  return F2CCell(F2CLinearRing({
      F2CPoint(x0, y0), F2CPoint(x1, y0), F2CPoint(x1, y1),
      F2CPoint(x0, y1), F2CPoint(x0, y0)}));
}
}  // namespace

TEST(fields2cover_utils_prepared_area, holdsWhatTheCellsHold) {
  const F2CCells cells {square(0, 0, 10, 10)};
  const f2c::PreparedArea area(cells);

  EXPECT_FALSE(area.isEmpty());
  EXPECT_TRUE(area.holds(5, 5));
  EXPECT_TRUE(area.holds(0.01, 0.01));
  EXPECT_FALSE(area.holds(-1, 5));
  EXPECT_FALSE(area.holds(15, 5));
  EXPECT_FALSE(area.holds(5, -1));

  // The same answers the cells themselves give, away from the ring.
  for (const auto& p : {F2CPoint(5, 5), F2CPoint(-1, 5), F2CPoint(15, 15)}) {
    EXPECT_EQ(area.holds(p.getX(), p.getY()), cells.isPointIn(p))
        << "at " << p.getX() << "," << p.getY();
  }
}

TEST(fields2cover_utils_prepared_area, doesNotHoldAHole) {
  F2CCell donut = square(0, 0, 10, 10);
  donut.addRing(F2CLinearRing({
      F2CPoint(4, 4), F2CPoint(6, 4), F2CPoint(6, 6),
      F2CPoint(4, 6), F2CPoint(4, 4)}));
  const f2c::PreparedArea area(F2CCells{donut});

  EXPECT_TRUE(area.holds(1, 1));
  EXPECT_FALSE(area.holds(5, 5)) << "the hole is not ground";
}

TEST(fields2cover_utils_prepared_area, holdsEveryPiece) {
  F2CCells two {square(0, 0, 10, 10)};
  two.addGeometry(square(20, 0, 30, 10));
  const f2c::PreparedArea area(two);

  EXPECT_TRUE(area.holds(5, 5));
  EXPECT_TRUE(area.holds(25, 5));
  EXPECT_FALSE(area.holds(15, 5)) << "the gap between them is not ground";
}

TEST(fields2cover_utils_prepared_area, emptyHoldsNothing) {
  const f2c::PreparedArea area;
  EXPECT_TRUE(area.isEmpty());
  EXPECT_FALSE(area.holds(0, 0));
}

TEST(fields2cover_utils_prepared_area, handsBackEveryRingEdge) {
  F2CCell donut = square(0, 0, 10, 10);
  donut.addRing(F2CLinearRing({
      F2CPoint(4, 4), F2CPoint(6, 4), F2CPoint(6, 6),
      F2CPoint(4, 6), F2CPoint(4, 4)}));
  const f2c::PreparedArea area(F2CCells{donut});

  std::vector<std::array<double, 4>> edges;
  area.edges(&edges);
  // Four sides of the square and four of the hole.
  EXPECT_EQ(edges.size(), 8);
}
