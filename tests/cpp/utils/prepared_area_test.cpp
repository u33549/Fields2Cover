//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <gtest/gtest.h>
#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "fields2cover/utils/transformation.h"
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

TEST(fields2cover_utils_prepared_area, agreesWithTheCellsAwayFromTheirRings) {
  // The point of the type is to answer what Cells::isPointIn answers, only
  // without going to GEOS. On a real, heavily digitised border the two have to
  // agree everywhere that is not on a ring, where the answer is undefined
  // here and GEOS calls it outside.
  std::ifstream f(std::string(DATA_PATH) + "ee_field_130.wkt");
  ASSERT_TRUE(f.is_open()) << "ee_field_130.wkt not found";
  std::stringstream ss;
  ss << f.rdbuf();
  F2CCell raw;
  raw.importFromWkt(ss.str());
  ASSERT_GT(raw.area(), 0);
  F2CField field(F2CCells(raw), "ee_field_130");
  field.setCRS("EPSG:4326");
  f2c::Transform::transformToUTM(field);
  const F2CCells cells = field.getField();

  F2CMultiLineString rings;
  for (size_t i = 0; i < cells.size(); ++i) {
    const F2CCell c = cells.getGeometry(i);
    rings.addGeometry(F2CLineString(c.getExteriorRing()));
    for (size_t r = 1; r < c.size(); ++r) {
      rings.addGeometry(F2CLineString(c.getInteriorRing(r - 1)));
    }
  }

  const f2c::PreparedArea area(cells);
  const F2CPoint lo = cells.getDimMinX() < 1e17
      ? F2CPoint(cells.getDimMinX(), cells.getDimMinY())
      : F2CPoint(0, 0);
  const double w = cells.getDimMaxX() - cells.getDimMinX();
  const double h = cells.getDimMaxY() - cells.getDimMinY();

  int checked = 0, disagreed = 0, inside = 0;
  for (int i = 0; i <= 80; ++i) {
    for (int j = 0; j <= 80; ++j) {
      const F2CPoint p(lo.getX() + w * i / 80.0, lo.getY() + h * j / 80.0);
      if (p.distance(rings) < 1e-3) {
        continue;      // on a ring: undefined here, outside for GEOS
      }
      const bool ours = area.holds(p.getX(), p.getY());
      if (ours != cells.isPointIn(p)) {
        ++disagreed;
      }
      if (ours) {
        ++inside;
      }
      ++checked;
    }
  }
  EXPECT_GT(checked, 5000) << "the grid has to actually cover the field";
  EXPECT_GT(inside, 500) << "and land inside it often enough to mean something";
  EXPECT_EQ(disagreed, 0);
}
