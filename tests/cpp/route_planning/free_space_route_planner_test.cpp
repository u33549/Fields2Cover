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

TEST(fields2cover_rp_free_space, turn_room_says_how_much_and_whether) {
  // The same corner, asked three ways: unset, the clearance answers for it;
  // zero leaves the connection where the graph put it; and a value of its own
  // is what the corner is held clear by.
  const F2CCells crop = block(0.0, 0.0, 194.0, 194.0);
  const F2CCells ground = block(0.0, 0.0, 200.0, 200.0).difference(crop);
  F2CSwaths sw;
  sw.emplace_back(F2CLineString({F2CPoint(197, 5), F2CPoint(197, 60)}), 2.0);
  sw.emplace_back(F2CLineString({F2CPoint(60, 197), F2CPoint(5, 197)}), 2.0);
  F2CSwathsByCells swaths;
  swaths.emplace_back(sw);

  auto cornerRoom = [&](bool set_room, double room) {
    f2c::rp::FreeSpaceRoutePlanner planner;
    planner.setClearance(6.0);
    if (set_room) {
      planner.setTurnRoom(room);
    }
    const F2CRoute route = planner.genRoute(ground, swaths, false, 1e-4);
    double worst = 1e18;
    for (const F2CMultiPoint& mp : route.getConnections()) {
      for (size_t i = 1; i + 1 < mp.size(); ++i) {
        worst = std::min(worst, mp.getGeometry(i).distance(crop));
      }
    }
    return worst;
  };

  f2c::rp::FreeSpaceRoutePlanner unset;
  unset.setClearance(6.0);
  EXPECT_DOUBLE_EQ(unset.getTurnRoom(), 6.0);
  f2c::rp::FreeSpaceRoutePlanner own;
  own.setClearance(6.0);
  own.setTurnRoom(3.0);
  EXPECT_DOUBLE_EQ(own.getTurnRoom(), 3.0);

  EXPECT_GE(cornerRoom(false, 0.0), 1.5);   // the clearance answers for it
  EXPECT_LT(cornerRoom(true, 0.0), 0.05);   // turned off: the graph's own corner
  const double wide = cornerRoom(true, 6.0);
  const double narrow = cornerRoom(true, 3.0);
  EXPECT_GT(narrow, 0.05);
  EXPECT_GT(wide, narrow + 0.5);            // half the radius asks for less room
}

TEST(fields2cover_rp_free_space, a_connection_leaves_a_swath_along_its_own_axis) {
  // A swath ends on the border of the ground, as swath ends do: the machine is
  // still on its own line there. The graph joined that end straight to whatever
  // corner was nearest, so the first metres of the connection cut across the
  // crop's edge at an angle no machine leaving the swath can hold.
  const F2CCells crop = block(0.0, 0.0, 100.0, 100.0);
  const F2CCells ground = block(-10.0, -10.0, 110.0, 110.0).difference(crop);

  // Two swaths ending on the bottom edge, far apart, so the connection between
  // them runs the length of the headland.
  F2CSwaths sw;
  sw.emplace_back(F2CLineString({F2CPoint(20, 40), F2CPoint(20, 0)}), 2.0);
  sw.emplace_back(F2CLineString({F2CPoint(80, 0), F2CPoint(80, 40)}), 2.0);
  F2CSwathsByCells swaths;
  swaths.emplace_back(sw);

  f2c::rp::FreeSpaceRoutePlanner planner;
  planner.setClearance(6.0);
  const F2CRoute route = planner.genRoute(ground, swaths, false, 1e-4);

  // Whatever the connection does next, it has to leave the swath end going the
  // way the swath goes: the first step is along that line, not across it.
  size_t checked = 0;
  for (const F2CMultiPoint& mp : route.getConnections()) {
    if (mp.size() < 2) { continue; }
    for (const size_t i : {size_t{0}, mp.size() - 1}) {
      const F2CPoint end = mp.getGeometry(i);
      if (end.distance(F2CPoint(20, 0)) > 1e-6 &&
          end.distance(F2CPoint(80, 0)) > 1e-6) {
        continue;             // not a swath end of ours
      }
      const F2CPoint next = mp.getGeometry(i == 0 ? 1 : mp.size() - 2);
      // The swaths run along x = 20 and x = 80, so leaving along the axis means
      // the step keeps that x.
      EXPECT_NEAR(next.getX(), end.getX(), 0.5)
          << "connection leaves the swath end at " << end.getX() << ", "
          << end.getY() << " sideways, towards " << next.getX() << ", "
          << next.getY();
      ++checked;
    }
  }
  EXPECT_GT(checked, 0u);
}

TEST(fields2cover_rp_free_space, a_swath_end_buried_in_the_crop_gets_no_entry) {
  // The entry set along a swath's axis crosses the band the ground was eroded
  // by, which the swath itself covers -- half a machine width, and the machine
  // there is still on what it has just cut. A swath that stops short of the
  // border leaves its end buried in the crop instead, and then its axis finds
  // ground again only on the far side. The entry is joined to the end
  // unconditionally, so the graph gains an edge driving through the crop.
  const F2CCells crop = block(0.0, 0.0, 100.0, 100.0);
  const F2CCells ground = block(-10.0, -10.0, 110.0, 110.0).difference(crop);

  F2CSwaths sw;
  // Reaches the border, as a swath does.
  sw.emplace_back(F2CLineString({F2CPoint(50, 100), F2CPoint(50, 0)}), 2.0);
  // Stops 10 m short of it: ground lies 10 m further along the same axis.
  sw.emplace_back(F2CLineString({F2CPoint(20, 100), F2CPoint(20, 10)}), 2.0);
  F2CSwathsByCells swaths;
  swaths.emplace_back(sw);

  f2c::rp::FreeSpaceRoutePlanner planner;
  planner.setClearance(6.0);
  planner.setSampleStep(0.25);
  const F2CGraph2D g = planner.createShortestGraph(ground, swaths, 1e-4);

  size_t checked = 0;
  for (const auto& from : g.getEdges()) {
    for (const auto& e : from.second) {
      const F2CPoint a = g.indexToNode(from.first);
      const F2CPoint b = g.indexToNode(e.first);
      for (int k = 1; k < 8; ++k) {
        const double u = k / 8.0;
        const F2CPoint m(a.getX() + (b.getX() - a.getX()) * u,
                         a.getY() + (b.getY() - a.getY()) * u);
        EXPECT_FALSE(crop.isPointIn(m))
            << "edge from " << a.getX() << ", " << a.getY() << " to "
            << b.getX() << ", " << b.getY() << " drives through the crop at "
            << m.getX() << ", " << m.getY();
      }
      ++checked;
    }
  }
  EXPECT_GT(checked, 0u);
}

TEST(fields2cover_rp_free_space, a_connection_is_aligned_though_its_ends_are_off_the_ground) {
  // In the field the ground is eroded away from the crop by half a machine
  // width, so a swath ends outside it: on the crop's border, with the ground
  // starting a metre further along the swath's own axis. The corner between
  // two such swaths is the ground's own inner corner, which a turn cannot
  // round without cutting the crop, so the alignment has to move it -- and it
  // used to compute that move and then throw it away, because it sampled the
  // leg reaching the swath end and of course found it off the ground.
  const F2CCells crop = block(0.0, 0.0, 194.0, 194.0);
  const F2CCells ground =
      block(0.0, 0.0, 220.0, 220.0).difference(block(0.0, 0.0, 195.0, 195.0));

  F2CSwaths sw;
  sw.emplace_back(F2CLineString({F2CPoint(5, 100), F2CPoint(194, 100)}), 2.0);
  sw.emplace_back(F2CLineString({F2CPoint(100, 194), F2CPoint(100, 5)}), 2.0);
  F2CSwathsByCells swaths;
  swaths.emplace_back(sw);

  f2c::rp::FreeSpaceRoutePlanner planner;
  planner.setClearance(6.0);
  planner.setClearanceCost(0.0);   // the geodesic, so the corner is the ground's
  planner.setCornerTolerance(0.0);
  planner.setSampleStep(0.25);
  const F2CRoute route = planner.genRoute(ground, swaths, false, 1e-4);

  size_t checked = 0;
  for (const F2CMultiPoint& mp : route.getConnections()) {
    if (mp.size() != 5) {
      continue;                  // swath end, step in, corner, step in, end
    }
    // Both ends are swath ends, on the crop's border, and stay there.
    EXPECT_NEAR(mp.getGeometry(0).distance(crop), 0.0, 1e-6);
    EXPECT_NEAR(mp.getGeometry(4).distance(crop), 0.0, 1e-6);
    // The ground's inner corner stands 1.41 m from the crop; a quarter turn of
    // radius 6 asks for 6 (1/cos 41 - 1) = 2.0 m more along the bisector,
    // which puts the corner 3.41 m out. Left unaligned it sits at 1.41.
    const F2CPoint corner = mp.getGeometry(2);
    EXPECT_GT(corner.distance(crop), 2.5)
        << "corner left at " << corner.getX() << ", " << corner.getY()
        << ": the alignment was computed and then thrown away";
    EXPECT_LT(corner.distance(F2CPoint(195.0, 195.0)), 3.0)
        << "corner moved further than the turn asked for";
    ++checked;
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
