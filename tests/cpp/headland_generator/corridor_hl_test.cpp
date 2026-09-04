//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "fields2cover/types.h"
#include "fields2cover/decomposition/trapezoidal_decomp.h"
#include "fields2cover/headland_generator/constant_headland.h"
#include "fields2cover/headland_generator/corridor_headland.h"
#include "fields2cover/path_planning/dubins_curves.h"
#include "fields2cover/path_planning/reeds_shepp_curves.h"

TEST(fields2cover_hl_corridor_gen, onlyBetweenTouchingCells) {
  f2c::hg::CorridorHL corridor;
  // Two 50x100 cells sharing the x=50 border: a 2 m corridor over a 100 m
  // border removes 200 m2, and both cells stay in one piece.
  F2CCells cells;
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0,0), F2CPoint(50,0), F2CPoint(50,100), F2CPoint(0,100), F2CPoint(0,0)})));
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(50,0), F2CPoint(100,0), F2CPoint(100,100),
      F2CPoint(50,100), F2CPoint(50,0)})));

  // Both cells are the same size, so they share the corridor evenly.
  F2CCells carved = corridor.generateHeadlands(cells, 2.0);
  EXPECT_EQ(carved.size(), 2);
  EXPECT_NEAR(cells.area() - carved.area(), 200, 1e-3);
  EXPECT_NEAR(carved.getGeometry(0).area(), 49 * 100, 1e-3);
  EXPECT_NEAR(carved.getGeometry(1).area(), 49 * 100, 1e-3);
}

TEST(fields2cover_hl_corridor_gen, leavesLoneCellsAlone) {
  f2c::hg::CorridorHL corridor;
  F2CCell left {F2CLinearRing({
      F2CPoint(0,0), F2CPoint(50,0), F2CPoint(50,100), F2CPoint(0,100), F2CPoint(0,0)})};

  // Cells across a gap wider than the corridor share no border.
  F2CCells apart;
  apart.addGeometry(left);
  apart.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(60,0), F2CPoint(110,0), F2CPoint(110,100),
      F2CPoint(60,100), F2CPoint(60,0)})));
  EXPECT_NEAR(corridor.generateHeadlands(apart, 2.0).area(), apart.area(), 1e-3);

  // A single cell has no neighbour at all.
  F2CCells one {left};
  EXPECT_NEAR(corridor.generateHeadlands(one, 2.0).area(), one.area(), 1e-3);
}

TEST(fields2cover_hl_corridor_gen, ignoresCornerContact) {
  f2c::hg::CorridorHL corridor;
  // Cells that meet at a single corner share no border, so nothing is carved.
  F2CCells cells;
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0,0), F2CPoint(10,0), F2CPoint(10,10), F2CPoint(0,10), F2CPoint(0,0)})));
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(10,10), F2CPoint(20,10), F2CPoint(20,20),
      F2CPoint(10,20), F2CPoint(10,10)})));

  EXPECT_NEAR(corridor.generateHeadlands(cells, 2.0).area(), cells.area(), 1e-6);
}

TEST(fields2cover_hl_corridor_gen, onPartiallyTouchingEdge) {
  f2c::hg::CorridorHL corridor;
  // Only the right half of the lower cell's top edge is shared. An edge test
  // that reads a whole edge at a time misses this and carves nothing.
  F2CCells cells;
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0,0), F2CPoint(20,0), F2CPoint(20,10), F2CPoint(0,10), F2CPoint(0,0)})));
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(10,10), F2CPoint(20,10), F2CPoint(20,20),
      F2CPoint(10,20), F2CPoint(10,10)})));

  // Only the shared part is cut: 10 m of border x 2 m = 20 m2, nothing more.
  // The smaller cell gives it up, so the larger one is untouched.
  F2CCells carved = corridor.generateHeadlands(cells, 2.0);
  EXPECT_EQ(carved.size(), 2);
  EXPECT_NEAR(cells.area() - carved.area(), 20, 1e-2);
  EXPECT_NEAR(carved.getGeometry(0).area(), 200, 1e-2);
  EXPECT_NEAR(carved.getGeometry(1).area(), 80, 1e-2);
}

TEST(fields2cover_hl_corridor_gen, keepsBothHalvesOfASplitCell) {
  f2c::hg::CorridorHL corridor;
  // The corridor cuts through the bridge of the U, leaving two arms. Both must
  // survive: dropping the extra piece would silently lose half the cell.
  F2CCells cells;
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0,0), F2CPoint(20,0), F2CPoint(20,10), F2CPoint(13,10),
      F2CPoint(13,2), F2CPoint(7,2), F2CPoint(7,10), F2CPoint(0,10), F2CPoint(0,0)})));
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(7,2), F2CPoint(13,2), F2CPoint(13,100), F2CPoint(7,100), F2CPoint(7,2)})));

  // The U has the shorter border, so it gives up the corridor and is split.
  F2CCells carved = corridor.generateHeadlands(cells, 6.0);
  EXPECT_EQ(carved.size(), 3);
  EXPECT_NEAR(carved.getGeometry(0).area(), 22, 1e-2);
  EXPECT_NEAR(carved.getGeometry(1).area(), 22, 1e-2);
  EXPECT_NEAR(carved.getGeometry(2).area(), 588, 1e-2);
}

TEST(fields2cover_hl_corridor_gen, leavesNoHairlineSpur) {
  f2c::hg::CorridorHL corridor;
  // A boolean difference can leave a zero-width spur on the outline. It keeps
  // the area right but derails whatever reads the shape next -- the swath
  // generator picks a skewed angle -- so carving must not hand one back:
  // opening the result by a hair must not change its outline.
  F2CCell field(F2CLinearRing({
      F2CPoint(0,0), F2CPoint(120,0), F2CPoint(120,90), F2CPoint(0,90), F2CPoint(0,0)}));
  field.addRing(F2CLinearRing({
      F2CPoint(45,35), F2CPoint(75,35), F2CPoint(75,60), F2CPoint(45,60), F2CPoint(45,35)}));

  f2c::decomp::TrapezoidalDecomp decomp;
  decomp.setSplitAngle(0.0);
  F2CCells carved = corridor.generateHeadlands(decomp.decompose(
      f2c::hg::ConstHL().generateHeadlands(F2CCells(field), 3.75)), 2.5);

  for (size_t i = 0; i < carved.size(); ++i) {
    const F2CCell cell = carved.getGeometry(i);
    F2CCells opened = F2CCells(cell).buffer(-1e-9).buffer(1e-9);
    ASSERT_EQ(opened.size(), 1);
    EXPECT_EQ(opened.getGeometry(0).getGeometry(0).size(),
              cell.getGeometry(0).size()) << "cell " << i << " has a spur";
  }
}

TEST(fields2cover_hl_corridor_gen, treatsEquallySizedCellsAlike) {
  f2c::hg::CorridorHL corridor;
  // Five equal slices meeting at the centre. Their perimeters are equal in
  // theory but differ in the last few bits, so a strict comparison lets that
  // noise pick which cell gives up the corridor -- and some borders end up
  // with no corridor at all while others get a full one.
  const double r = 50.0;
  F2CCells pie;
  for (int k = 0; k < 5; ++k) {
    const double a0 = 2 * M_PI * k / 5, a1 = 2 * M_PI * (k + 1) / 5;
    pie.addGeometry(F2CCell(F2CLinearRing({
        F2CPoint(0, 0), F2CPoint(r * cos(a0), r * sin(a0)),
        F2CPoint(r * cos((a0 + a1) / 2) * 1.02, r * sin((a0 + a1) / 2) * 1.02),
        F2CPoint(r * cos(a1), r * sin(a1)), F2CPoint(0, 0)})));
  }

  F2CCells carved = corridor.generateHeadlands(pie, 4.0);
  ASSERT_EQ(carved.size(), 5);
  for (size_t i = 1; i < carved.size(); ++i) {
    EXPECT_NEAR(carved.getGeometry(i).area(), carved.getGeometry(0).area(), 1e-3)
        << "slice " << i << " lost a different amount than slice 0";
  }
}

TEST(fields2cover_hl_corridor_gen, sharesTellWhoGivesTheCorridor) {
  f2c::hg::CorridorHL corridor;
  // The lower cell is the larger one, and only 10 m of its top edge is shared.
  F2CCells cells;
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0,0), F2CPoint(20,0), F2CPoint(20,10), F2CPoint(0,10), F2CPoint(0,0)})));
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(10,10), F2CPoint(20,10), F2CPoint(20,20),
      F2CPoint(10,20), F2CPoint(10,10)})));

  const std::vector<f2c::hg::CorridorShare> shares = corridor.corridorShares(cells);
  // One row per side of the border, and both know the border is 10 m long.
  ASSERT_EQ(shares.size(), 2);
  for (const auto& s : shares) {
    EXPECT_FALSE(s.same_size);
    EXPECT_NEAR(s.shared_length, 10, 1e-2);
    EXPECT_NEAR(s.shared_border.length(), 10, 1e-2);
    EXPECT_NEAR(s.perimeter_i, cells.getCellBorder(s.cell_i).length(), 1e-9);
    EXPECT_NEAR(s.perimeter_k, cells.getCellBorder(s.cell_k).length(), 1e-9);
  }
  // The smaller cell gives the whole corridor, the larger one gives nothing.
  const auto& from_big = shares[0].cell_i == 0 ? shares[0] : shares[1];
  const auto& from_small = shares[0].cell_i == 0 ? shares[1] : shares[0];
  EXPECT_EQ(from_big.cell_i, 0);
  EXPECT_EQ(from_big.cell_k, 1);
  EXPECT_NEAR(from_big.share, 0.0, 1e-9);
  EXPECT_EQ(from_small.cell_i, 1);
  EXPECT_EQ(from_small.cell_k, 0);
  EXPECT_NEAR(from_small.share, 1.0, 1e-9);
}

TEST(fields2cover_hl_corridor_gen, sharesSplitEvenlyBetweenEqualCells) {
  f2c::hg::CorridorHL corridor;
  F2CCells cells;
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0,0), F2CPoint(50,0), F2CPoint(50,100), F2CPoint(0,100), F2CPoint(0,0)})));
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(50,0), F2CPoint(100,0), F2CPoint(100,100),
      F2CPoint(50,100), F2CPoint(50,0)})));

  const std::vector<f2c::hg::CorridorShare> shares = corridor.corridorShares(cells);
  ASSERT_EQ(shares.size(), 2);
  for (const auto& s : shares) {
    EXPECT_TRUE(s.same_size);
    EXPECT_NEAR(s.share, 0.5, 1e-9);
    EXPECT_NEAR(s.shared_length, 100, 1e-2);
  }
}

TEST(fields2cover_hl_corridor_gen, sharesSkipCellsThatDoNotShareABorder) {
  f2c::hg::CorridorHL corridor;
  // Cells apart, and cells meeting at a single corner: neither shares a border.
  F2CCells apart;
  apart.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0,0), F2CPoint(50,0), F2CPoint(50,100), F2CPoint(0,100), F2CPoint(0,0)})));
  apart.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(60,0), F2CPoint(110,0), F2CPoint(110,100),
      F2CPoint(60,100), F2CPoint(60,0)})));
  EXPECT_TRUE(corridor.corridorShares(apart).empty());

  F2CCells corner;
  corner.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0,0), F2CPoint(10,0), F2CPoint(10,10), F2CPoint(0,10), F2CPoint(0,0)})));
  corner.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(10,10), F2CPoint(20,10), F2CPoint(20,20),
      F2CPoint(10,20), F2CPoint(10,10)})));
  EXPECT_TRUE(corridor.corridorShares(corner).empty());

  F2CCells one {apart.getGeometry(0)};
  EXPECT_TRUE(corridor.corridorShares(one).empty());
}

TEST(fields2cover_hl_corridor_gen, sharesAccountForEveryCarvedCorridor) {
  f2c::hg::CorridorHL corridor;
  // Three cells in a row, the middle one the widest: the outer cells each give
  // a full corridor, and the shares add up to the area the carving removed.
  F2CCells cells;
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0,0), F2CPoint(10,0), F2CPoint(10,100), F2CPoint(0,100), F2CPoint(0,0)})));
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(10,0), F2CPoint(40,0), F2CPoint(40,100),
      F2CPoint(10,100), F2CPoint(10,0)})));
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(40,0), F2CPoint(50,0), F2CPoint(50,100),
      F2CPoint(40,100), F2CPoint(40,0)})));

  const std::vector<f2c::hg::CorridorShare> shares = corridor.corridorShares(cells);
  ASSERT_EQ(shares.size(), 4);
  double carved_by_shares {0.0};
  for (const auto& s : shares) {
    carved_by_shares += s.share * s.shared_length * 2.0;
  }
  const double carved = cells.area() - corridor.generateHeadlands(cells, 2.0).area();
  EXPECT_NEAR(carved_by_shares, carved, 1e-2);
  EXPECT_NEAR(carved, 2 * 100 * 2.0, 1e-2);
}

TEST(fields2cover_hl_corridor_gen, sharesSymmetricModeAlwaysSplitsEvenly) {
  f2c::hg::CorridorHL corridor;
  // Same field as sharesTellWhoGivesTheCorridor: sizes differ, so the default
  // rule would give the whole corridor to the smaller cell.
  F2CCells cells;
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0,0), F2CPoint(20,0), F2CPoint(20,10), F2CPoint(0,10), F2CPoint(0,0)})));
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(10,10), F2CPoint(20,10), F2CPoint(20,20),
      F2CPoint(10,20), F2CPoint(10,10)})));

  const std::vector<f2c::hg::CorridorShare> shares = corridor.corridorShares(
      cells, f2c::hg::CorridorShareMode::SYMMETRIC);
  ASSERT_EQ(shares.size(), 2);
  for (const auto& s : shares) {
    EXPECT_FALSE(s.same_size);
    EXPECT_NEAR(s.share, 0.5, 1e-9);
  }
}

TEST(fields2cover_hl_corridor_gen, defaultShareModeIsAsymmetric) {
  f2c::hg::CorridorHL corridor;
  EXPECT_EQ(corridor.getShareMode(), f2c::hg::CorridorShareMode::ASYMMETRIC);
}

TEST(fields2cover_hl_corridor_gen, generateHeadlandsUsesConfiguredShareMode) {
  f2c::hg::CorridorHL corridor;
  // Two cells of the same size, so a symmetric split reduces both by the same
  // amount instead of taking the whole corridor out of just one.
  F2CCells cells;
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0,0), F2CPoint(50,0), F2CPoint(50,100), F2CPoint(0,100), F2CPoint(0,0)})));
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(50,0), F2CPoint(100,0), F2CPoint(100,100),
      F2CPoint(50,100), F2CPoint(50,0)})));

  corridor.setShareMode(f2c::hg::CorridorShareMode::SYMMETRIC);
  EXPECT_EQ(corridor.getShareMode(), f2c::hg::CorridorShareMode::SYMMETRIC);
  F2CCells carved = corridor.generateHeadlands(cells, 2.0);
  EXPECT_NEAR(carved.getGeometry(0).area(), 49 * 100, 1e-3);
  EXPECT_NEAR(carved.getGeometry(1).area(), 49 * 100, 1e-3);
}

TEST(fields2cover_hl_corridor_gen, sharesTreatACornerAsNoBorder) {
  f2c::hg::CorridorHL corridor;
  // Five equal slices meeting at the centre. Each slice borders two others;
  // the other two only reach the same point. Buffering a neighbour makes a
  // corner look like a few millimetres of shared border, and a table built on
  // that would call all four neighbours.
  const double r = 50.0;
  F2CCells pie;
  for (int k = 0; k < 5; ++k) {
    const double a0 = 2 * M_PI * k / 5, a1 = 2 * M_PI * (k + 1) / 5;
    pie.addGeometry(F2CCell(F2CLinearRing({
        F2CPoint(0, 0), F2CPoint(r * cos(a0), r * sin(a0)),
        F2CPoint(r * cos((a0 + a1) / 2) * 1.02, r * sin((a0 + a1) / 2) * 1.02),
        F2CPoint(r * cos(a1), r * sin(a1)), F2CPoint(0, 0)})));
  }

  const std::vector<f2c::hg::CorridorShare> shares = corridor.corridorShares(pie);
  // Five borders, seen from both sides.
  EXPECT_EQ(shares.size(), 10);
  for (const auto& s : shares) {
    EXPECT_TRUE(s.same_size);
    EXPECT_NEAR(s.share, 0.5, 1e-9);
    EXPECT_NEAR(s.shared_length, r, 1e-2) << "cells " << s.cell_i << "-" << s.cell_k;
  }
}

// Two 100x50 cells sharing the whole y=50 border. Same perimeter, so each
// gives half of whatever the corridor turns out to need.
namespace {
F2CCells twoCellsSharingAHorizontalBorder() {
  F2CCells cells;
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0,0), F2CPoint(100,0), F2CPoint(100,50),
      F2CPoint(0,50), F2CPoint(0,0)})));
  cells.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0,50), F2CPoint(100,50), F2CPoint(100,100),
      F2CPoint(0,100), F2CPoint(0,50)})));
  return cells;
}

// A planner that cannot turn at all, to stand in for one that fails on a
// turn it is handed.
class NoTurnPlanner : public f2c::pp::TurningBase {
 public:
  F2CPath createSimpleTurn(const F2CRobot&, double, double, double) override {
    return F2CPath();
  }
};
}  // namespace

TEST(fields2cover_hl_corridor_gen, turnExtentIsOneRadiusWhenTheSwathsLeaveRoom) {
  f2c::hg::CorridorHL corridor;
  f2c::pp::DubinsCurves dubins;
  F2CRobot robot(2.0, 10.0);
  robot.setMinTurningRadius(2.0);

  // Swaths 10 m apart with a 2 m turning radius leave room to spare, so the
  // turn is two quarter circles joined by a straight and only reaches one
  // radius past the swath ends -- half of what the classic bound assumes.
  EXPECT_NEAR(corridor.turnExtent(robot, dubins),
      robot.getMinTurningRadius(), 1e-2);
  EXPECT_LT(corridor.turnExtent(robot, dubins),
      2.0 * robot.getMinTurningRadius());
}

TEST(fields2cover_hl_corridor_gen, turnExtentPassesTheBoundWhereTheTurnHasToLoop) {
  f2c::hg::CorridorHL corridor;
  f2c::pp::DubinsCurves dubins;
  F2CRobot robot(2.0, 6.0);
  robot.setMinTurningRadius(5.0);
  const double r = robot.getMinTurningRadius();
  const double d = robot.getCovWidth();

  // Swaths closer together than twice the radius leave the turn no room to
  // join them directly: it loops, and reaches r + sqrt(4r^2 - (d/2 + r)^2)
  // past their ends -- further than the classic bound of twice the radius,
  // which is the case that bound gets wrong.
  const double loop = r + std::sqrt(4.0 * r * r - std::pow(0.5 * d + r, 2));
  EXPECT_GT(loop, 2.0 * r);
  EXPECT_NEAR(corridor.turnExtent(robot, dubins), loop, 1e-2);
}

TEST(fields2cover_hl_corridor_gen, turnExtentIsZeroForATurnThatBacksUp) {
  f2c::hg::CorridorHL corridor;
  f2c::pp::ReedsSheppCurves reeds_shepp;
  F2CRobot robot(2.0, 2.0);
  robot.setMinTurningRadius(5.0);

  // Reeds-Shepp may back up rather than drive round, and then the turn never
  // reaches past the end of the swaths at all. That is an answer, not a
  // failure: no room is needed there.
  EXPECT_NEAR(corridor.turnExtent(robot, reeds_shepp), 0.0, 1e-2);
}

TEST(fields2cover_hl_corridor_gen, turnExtentFallsBackOnTheBoundWithoutATurn) {
  f2c::hg::CorridorHL corridor;
  NoTurnPlanner no_turn;
  F2CRobot robot(2.0, 6.0);
  robot.setMinTurningRadius(5.0);

  // Nothing came back to measure, so the classic bound stands.
  EXPECT_NEAR(corridor.turnExtent(robot, no_turn),
      2.0 * robot.getMinTurningRadius(), 1e-9);
}

TEST(fields2cover_hl_corridor_gen, corridorIsAsDeepAsThePlannedTurn) {
  f2c::hg::CorridorHL corridor;
  f2c::pp::DubinsCurves dubins;
  F2CCells cells = twoCellsSharingAHorizontalBorder();
  F2CRobot robot(2.0, 10.0);
  robot.setMinTurningRadius(2.0);

  // The corridor over the 100 m border is the turn's reach plus half the
  // robot, and the two cells are the same size so each gives half of it.
  const double width =
      corridor.turnExtent(robot, dubins) + 0.5 * robot.getWidth();
  F2CCells carved = corridor.generateHeadlands(cells, robot, dubins);
  EXPECT_EQ(carved.size(), 2);
  EXPECT_NEAR(cells.area() - carved.area(), 100 * width, 1e-2);
  EXPECT_NEAR(carved.getGeometry(0).area(), 100 * (50 - 0.5 * width), 1e-2);
}

TEST(fields2cover_hl_corridor_gen, dropsAPieceNarrowerThanTheImplement) {
  f2c::hg::CorridorHL corridor;
  f2c::pp::DubinsCurves dubins;
  F2CRobot robot(1.0, 4.0);   // body 1 m, implement 4 m
  robot.setMinTurningRadius(1.0);

  // Two groups of three stacked cells, far enough apart to share no border.
  // The middle band of each group is the smallest of its group, so it gives
  // up the whole corridor on both of its borders.
  //   band of  5 m -> 2 m left: narrower than the implement, and wider than
  //                    the body, so only the implement's width can decide it
  //   band of  9 m -> 6 m left: wide enough to cover, and must survive
  auto band = [](double x0, double y0, double y1) {
    return F2CCell(F2CLinearRing({
        F2CPoint(x0, y0), F2CPoint(x0 + 100, y0), F2CPoint(x0 + 100, y1),
        F2CPoint(x0, y1), F2CPoint(x0, y0)}));
  };
  F2CCells cells;
  cells.addGeometry(band(0, 0, 50));
  cells.addGeometry(band(0, 50, 55));     // 5 m
  cells.addGeometry(band(0, 55, 105));
  cells.addGeometry(band(200, 0, 50));
  cells.addGeometry(band(200, 50, 59));   // 9 m
  cells.addGeometry(band(200, 59, 109));

  const double width =
      corridor.turnExtent(robot, dubins) + 0.5 * robot.getWidth();
  ASSERT_NEAR(width, 1.5, 1e-2);

  F2CCells carved = corridor.generateHeadlands(cells, robot, dubins);
  // The 2 m band is gone; the 6 m one is not.
  EXPECT_EQ(carved.size(), 5);
  EXPECT_NEAR(carved.area(), 4 * 100 * 50 + 100 * 6, 1e-2);
  double thinnest = 1e9;
  for (size_t i = 0; i < carved.size(); ++i) {
    thinnest = std::min(thinnest, carved.getGeometry(i).area() / 100.0);
    EXPECT_GT(F2CCells::buffer(
        carved.getGeometry(i), -0.5 * robot.getCovWidth()).area(), 0.0);
  }
  EXPECT_NEAR(thinnest, 6.0, 1e-2);
}

TEST(fields2cover_hl_corridor_gen,
    angsNarrowTheCorridorWhenSwathsRunAlongTheBorder) {
  f2c::hg::CorridorHL corridor;
  f2c::pp::DubinsCurves dubins;
  F2CCells cells = twoCellsSharingAHorizontalBorder();
  F2CRobot robot(2.0, 10.0);
  robot.setMinTurningRadius(2.0);

  // Both cells' swaths run along the border (track angle 0, same as the
  // border itself): no turn happens there and nothing drives there either,
  // so no corridor is opened and the two cells come back as one piece.
  const std::vector<double> angs {0.0, 0.0};
  F2CCells carved = corridor.generateHeadlands(cells, robot, dubins, angs);
  EXPECT_NEAR(carved.area(), cells.area(), 1e-2);
  EXPECT_EQ(carved.size(), 1);
}

TEST(fields2cover_hl_corridor_gen,
    angsMatchTheUniformCorridorWhenSwathsMeetTheBorderHeadOn) {
  f2c::hg::CorridorHL corridor;
  f2c::pp::DubinsCurves dubins;
  F2CCells cells = twoCellsSharingAHorizontalBorder();
  F2CRobot robot(2.0, 10.0);
  robot.setMinTurningRadius(2.0);

  // Both cells' swaths meet the border head-on (track angle pi/2, square to
  // the border): the worst case turnExtent() already covers, so the
  // angle-aware corridor should come out exactly as deep as the uniform one.
  const std::vector<double> angs {M_PI_2, M_PI_2};
  F2CCells angled = corridor.generateHeadlands(cells, robot, dubins, angs);
  F2CCells uniform = corridor.generateHeadlands(cells, robot, dubins);
  EXPECT_NEAR(angled.area(), uniform.area(), 1e-2);
}

TEST(fields2cover_hl_corridor_gen, angsWrongSizeThrows) {
  f2c::hg::CorridorHL corridor;
  f2c::pp::DubinsCurves dubins;
  F2CCells cells = twoCellsSharingAHorizontalBorder();
  F2CRobot robot(2.0, 10.0);
  robot.setMinTurningRadius(2.0);

  const std::vector<double> too_short {0.0};
  EXPECT_THROW(
      corridor.generateHeadlands(cells, robot, dubins, too_short),
      std::invalid_argument);
}

TEST(fields2cover_hl_corridor_gen, turnExtentAtRightAngleIsTheSquareCase) {
  f2c::hg::CorridorHL corridor;
  f2c::pp::DubinsCurves dubins;
  F2CRobot robot(2.0, 3.0);
  robot.setMinTurningRadius(2.0);

  EXPECT_NEAR(corridor.turnExtent(robot, dubins, M_PI_2),
      corridor.turnExtent(robot, dubins), 1e-6);
}

TEST(fields2cover_hl_corridor_gen, turnExtentStaysAboveTheTurningRadius) {
  f2c::hg::CorridorHL corridor;
  f2c::pp::DubinsCurves dubins;
  F2CRobot robot(2.0, 3.0);
  robot.setMinTurningRadius(2.0);

  // Swaths further and further from square to the border cross it further
  // apart, but the turn between them still has to swing round: it never
  // reaches less than the radius, so it never scales to nothing the way
  // sin(angle) does.
  for (int deg = 5; deg <= 90; deg += 5) {
    const double reach = corridor.turnExtent(robot, dubins, deg * M_PI / 180.0);
    EXPECT_GE(reach, robot.getMinTurningRadius() - 1e-6) << "at " << deg;
  }
  const double shallow = corridor.turnExtent(robot, dubins, 5.0 * M_PI / 180.0);
  EXPECT_GT(shallow,
      std::sin(5.0 * M_PI / 180.0) * corridor.turnExtent(robot, dubins));
}
