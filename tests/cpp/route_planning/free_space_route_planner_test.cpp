//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <gtest/gtest.h>
#include <vector>
#include "fields2cover/types.h"
#include "fields2cover/route_planning/free_space_route_planner.h"
#include "fields2cover.h"

namespace {

// A square field with a square of crop in the middle: the ground the robot may
// drive on is the ring between them, and the only way from one side to the
// other is around.
F2CCells ringGround(double outer, double inner) {
  F2CCell c(F2CLinearRing({
      F2CPoint(0, 0), F2CPoint(outer, 0), F2CPoint(outer, outer),
      F2CPoint(0, outer), F2CPoint(0, 0)}));
  const double lo = 0.5 * (outer - inner), hi = lo + inner;
  c.addRing(F2CLinearRing({
      F2CPoint(lo, lo), F2CPoint(lo, hi), F2CPoint(hi, hi),
      F2CPoint(hi, lo), F2CPoint(lo, lo)}));
  return F2CCells(c);
}

double pathLength(const std::vector<F2CPoint>& p) {
  double len = 0.0;
  for (size_t i = 0; i + 1 < p.size(); ++i) {
    len += p[i].distance(p[i + 1]);
  }
  return len;
}

F2CCells block(double x0, double y0, double x1, double y1) {
  return F2CCells(F2CCell(F2CLinearRing({
      F2CPoint(x0, y0), F2CPoint(x1, y0), F2CPoint(x1, y1),
      F2CPoint(x0, y1), F2CPoint(x0, y0)})));
}

}  // namespace

TEST(fields2cover_rp_free_space, every_edge_stays_on_the_ground) {
  const F2CCells ground = ringGround(100.0, 60.0);
  F2CSwathsByCells swaths;
  f2c::rp::FreeSpaceRoutePlanner planner;
  const F2CGraph2D g = planner.createShortestGraph(ground, swaths, 1e-4);

  EXPECT_GT(g.numNodes(), 0u);
  const auto edges = g.getEdges();
  size_t checked = 0;
  for (const auto& from : edges) {
    for (const auto& e : from.second) {
      const F2CPoint a = g.indexToNode(from.first);
      const F2CPoint b = g.indexToNode(e.first);
      // Sample the middle of the edge; the crop is a hole, so an edge that
      // crossed it would land inside.
      for (int k = 1; k < 8; ++k) {
        const double u = k / 8.0;
        const F2CPoint m(a.getX() + (b.getX() - a.getX()) * u,
                         a.getY() + (b.getY() - a.getY()) * u);
        EXPECT_TRUE(ground.isPointIn(m) || ground.isPointInBorder(m))
            << "edge leaves the drivable ground at " << m.getX()
            << ", " << m.getY();
      }
      ++checked;
    }
  }
  EXPECT_GT(checked, 0u);
}

TEST(fields2cover_rp_free_space, goes_around_the_crop_not_through_it) {
  const F2CCells ground = ringGround(100.0, 60.0);
  F2CSwathsByCells swaths;
  f2c::rp::FreeSpaceRoutePlanner planner;
  F2CGraph2D g = planner.createShortestGraph(ground, swaths, 1e-4);

  // Two corners of the field on opposite sides of the crop. The straight line
  // between them runs through the crop, so the path has to be longer than it.
  const F2CPoint a(0.0, 0.0), b(100.0, 100.0);
  ASSERT_TRUE(g.hasNode(a));
  ASSERT_TRUE(g.hasNode(b));
  const std::vector<F2CPoint> path = g.shortestPath(a, b);
  ASSERT_GE(path.size(), 2u);
  EXPECT_EQ(path.front(), a);
  EXPECT_EQ(path.back(), b);
  EXPECT_GT(pathLength(path), a.distance(b));
  EXPECT_EQ(planner.getComponentCount(), 1u);
}

TEST(fields2cover_rp_free_space, counts_the_pieces_the_ground_falls_into) {
  // Two squares that do not touch: nothing can be routed between them, and the
  // planner says so instead of leaving the caller to find out from the route.
  F2CCells ground;
  ground.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(0, 0), F2CPoint(10, 0), F2CPoint(10, 10),
      F2CPoint(0, 10), F2CPoint(0, 0)})));
  ground.addGeometry(F2CCell(F2CLinearRing({
      F2CPoint(50, 0), F2CPoint(60, 0), F2CPoint(60, 10),
      F2CPoint(50, 10), F2CPoint(50, 0)})));

  F2CSwathsByCells swaths;
  f2c::rp::FreeSpaceRoutePlanner planner;
  const F2CGraph2D g = planner.createShortestGraph(ground, swaths, 1e-4);
  EXPECT_EQ(planner.getComponentCount(), 2u);
  EXPECT_GT(g.numNodes(), 0u);
}

TEST(fields2cover_rp_free_space, swath_ends_are_reachable_nodes) {
  const F2CCells ground = ringGround(100.0, 60.0);
  // A swath lying in the ring, its ends on the ground.
  F2CSwaths sw;
  sw.emplace_back(F2CLineString({F2CPoint(5, 20), F2CPoint(5, 80)}), 2.0);
  F2CSwathsByCells swaths;
  swaths.emplace_back(sw);

  f2c::rp::FreeSpaceRoutePlanner planner;
  const F2CGraph2D g = planner.createShortestGraph(ground, swaths, 1e-4);
  EXPECT_TRUE(g.hasNode(F2CPoint(5, 20)));
  EXPECT_TRUE(g.hasNode(F2CPoint(5, 80)));
}

TEST(fields2cover_rp_free_space, clearance_is_priced_not_forbidden) {
  const F2CCells ground = ringGround(100.0, 60.0);
  F2CSwathsByCells swaths;

  f2c::rp::FreeSpaceRoutePlanner plain;
  const F2CGraph2D g0 = plain.createShortestGraph(ground, swaths, 1e-4);

  f2c::rp::FreeSpaceRoutePlanner kept_clear;
  kept_clear.setClearance(3.0);
  kept_clear.setClearanceCost(10.0);
  const F2CGraph2D g1 = kept_clear.createShortestGraph(ground, swaths, 1e-4);

  // The band changes what an edge costs, never whether it exists: the ground
  // stays one piece and no node is cut off.
  EXPECT_EQ(g0.numNodes(), g1.numNodes());
  EXPECT_EQ(kept_clear.getComponentCount(), 1u);
  EXPECT_DOUBLE_EQ(kept_clear.getClearance(), 3.0);
  EXPECT_DOUBLE_EQ(kept_clear.getClearanceCost(), 10.0);
}

TEST(fields2cover_rp_free_space, a_connection_leaves_room_for_its_turn_at_a_corner) {
  // The ground is an L 6 m wide around the corner of the crop, with a swath in
  // each arm, both far enough along it that the connection between them runs
  // nearly parallel to the borders and turns a right angle at the corner. A
  // geodesic hugs whatever it goes around, so that connection runs through the
  // corner itself -- and a turn of any radius has to cut into the crop to
  // follow it. A quarter turn of radius 6 asks for 6 (1 - cos 45) = 1.76 m.
  const F2CCells crop = block(0.0, 0.0, 194.0, 194.0);
  const F2CCells ground = block(0.0, 0.0, 200.0, 200.0).difference(crop);

  F2CSwaths sw;
  sw.emplace_back(F2CLineString({F2CPoint(197, 5), F2CPoint(197, 60)}), 2.0);
  sw.emplace_back(F2CLineString({F2CPoint(60, 197), F2CPoint(5, 197)}), 2.0);
  F2CSwathsByCells swaths;
  swaths.emplace_back(sw);

  f2c::rp::FreeSpaceRoutePlanner planner;
  planner.setClearance(6.0);   // the room a turn of that radius needs
  const F2CRoute route = planner.genRoute(ground, swaths, false, 1e-4);

  ASSERT_GT(route.sizeConnections(), 0u);
  size_t checked = 0;
  for (const F2CMultiPoint& mp : route.getConnections()) {
    for (size_t i = 0; i < mp.size(); ++i) {
      const F2CPoint p = mp.getGeometry(i);
      EXPECT_GE(p.distance(crop), 1.5)
          << "connection point " << p.getX() << ", " << p.getY()
          << " hugs the crop, leaving its turn nothing to cut";
      ++checked;
    }
  }
  EXPECT_GT(checked, 0u);
}

TEST(fields2cover_rp_free_space, empty_ground_is_not_a_crash) {
  F2CCells ground;
  F2CSwathsByCells swaths;
  f2c::rp::FreeSpaceRoutePlanner planner;
  const F2CGraph2D g = planner.createShortestGraph(ground, swaths, 1e-4);
  EXPECT_EQ(g.numNodes(), 0u);
  EXPECT_EQ(planner.getComponentCount(), 0u);
}
