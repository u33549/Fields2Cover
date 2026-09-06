//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <gtest/gtest.h>
#include <cmath>
#include "fields2cover/path_planning/constrained_turn.h"
#include "fields2cover/path_planning/dubins_curves.h"
#include "fields2cover/types.h"

namespace {

F2CRobot makeRobot() {
  F2CRobot robot(3.0, 6.0);
  robot.setCruiseVel(2.0);
  robot.setMaxCurv(1.0 / 6.0);
  return robot;
}

F2CCell rectangle(double x0, double y0, double x1, double y1) {
  F2CLinearRing ring;
  ring.addPoint(F2CPoint(x0, y0));
  ring.addPoint(F2CPoint(x1, y0));
  ring.addPoint(F2CPoint(x1, y1));
  ring.addPoint(F2CPoint(x0, y1));
  ring.addPoint(F2CPoint(x0, y0));
  return F2CCell(ring);
}

// The band that wraps the outside of a right-angled crop corner: the crop
// fills the first quadrant, the band is W wide around it.
F2CCells cornerBand(double width, double arm) {
  F2CCells band;
  band.addGeometry(rectangle(-width, -width, arm, 0.0));
  band.addGeometry(rectangle(-width, -width, 0.0, arm));
  return band.unionCascaded();
}

double lengthOutside(const F2CPath& path, const F2CCells& free_space) {
  double out = 0.0;
  const auto& states = path.getStates();
  for (size_t i = 0; i + 1 < states.size(); ++i) {
    const F2CPoint a = states[i].point, b = states[i + 1].point;
    const double len = a.distance(b);
    if (len <= 0.0) { continue; }
    const int n = std::max(1, static_cast<int>(std::ceil(len / 0.05)));
    for (int j = 0; j < n; ++j) {
      const double t = (j + 0.5) / n;
      if (!free_space.isPointIn(F2CPoint(
          a.getX() + (b.getX() - a.getX()) * t,
          a.getY() + (b.getY() - a.getY()) * t))) {
        out += len / n;
      }
    }
  }
  return out;
}

double minRadius(const F2CPath& path) {
  double r = 1e18;
  const auto& states = path.getStates();
  for (size_t i = 0; i + 1 < states.size(); ++i) {
    const double ds = states[i].point.distance(states[i + 1].point);
    const double da = F2CPoint::getAngleDiffAbs(
        states[i].angle, states[i + 1].angle);
    if (da > 1e-6 && ds > 1e-9) { r = std::min(r, ds / da); }
  }
  return r;
}

}  // namespace

TEST(fields2cover_pp_constrained_turn, corners_of_a_convex_region) {
  // A rectangle turns outward at every vertex, so it has nothing to wrap.
  const F2CCells box {rectangle(-50.0, 0.0, 50.0, 18.0)};
  EXPECT_TRUE(f2c::pp::ConstrainedTurn::concaveCorners(box).empty());
}

TEST(fields2cover_pp_constrained_turn, corners_of_a_wrapped_corner) {
  // The band around a right-angled crop corner turns inward exactly once.
  const auto corners = f2c::pp::ConstrainedTurn::concaveCorners(
      cornerBand(12.0, 80.0));
  ASSERT_EQ(corners.size(), 1);
  EXPECT_NEAR(corners[0].first.getX(), 0.0, 1e-6);
  EXPECT_NEAR(corners[0].first.getY(), 0.0, 1e-6);
  // The bisector points away from the crop, into the band.
  EXPECT_NEAR(F2CPoint::getAngleDiffAbs(corners[0].second, 5.0 * M_PI_4),
      0.0, 1e-6);
}

TEST(fields2cover_pp_constrained_turn, a_turn_that_fits_is_left_alone) {
  F2CRobot robot = makeRobot();
  f2c::pp::DubinsCurves dubins;
  f2c::pp::ConstrainedTurn constrained;
  // Two swath ends 18 m apart on one face of an 18 m band: a plain u-turn.
  const F2CCells band {rectangle(-60.0, 0.0, 60.0, 18.0)};
  const F2CPoint start(0.0, 0.0), end(18.0, 0.0);

  const F2CPath plain = dubins.createTurn(
      robot, start, M_PI_2, end, -M_PI_2);
  const auto res = constrained.createTurn(
      robot, band, start, M_PI_2, end, -M_PI_2, dubins);

  EXPECT_TRUE(res.inside);
  EXPECT_FALSE(res.used_waypoint);
  EXPECT_NEAR(res.path.length(), plain.length(), 1e-9);
}

TEST(fields2cover_pp_constrained_turn, a_turn_that_cuts_a_corner_is_repaired) {
  F2CRobot robot = makeRobot();
  f2c::pp::DubinsCurves dubins;
  f2c::pp::ConstrainedTurn constrained;
  const F2CCells band = cornerBand(12.0, 80.0);
  // One pose on each face of the crop corner, each perpendicular to its face:
  // out of the crop on the bottom, back into it on the left.
  const F2CPoint start(24.0, 0.0), end(0.0, 24.0);
  const double start_angle = -M_PI_2, end_angle = 0.0;

  const F2CPath plain = dubins.createTurn(
      robot, start, start_angle, end, end_angle);
  ASSERT_GT(lengthOutside(plain, band), 1.0)
      << "the plain turn is expected to cut across the crop";

  const auto res = constrained.createTurn(
      robot, band, start, start_angle, end, end_angle, dubins);

  EXPECT_TRUE(res.inside);
  EXPECT_TRUE(res.used_waypoint);
  EXPECT_LT(lengthOutside(res.path, band), 0.05);
  // Repairing must not cost drivability: the radius still holds, and the two
  // halves meet at one pose, so the ends are the ones asked for.
  EXPECT_GT(minRadius(res.path), 0.95 * robot.getMinTurningRadius());
  EXPECT_NEAR(res.path.getStates().front().point.distance(start), 0.0, 1e-6);
  EXPECT_NEAR(res.path.getStates().back().point.distance(end), 0.0, 1e-6);
  EXPECT_GT(res.path.length(), plain.length());
}

TEST(fields2cover_pp_constrained_turn, no_waypoint_can_widen_a_narrow_band) {
  F2CRobot robot = makeRobot();
  f2c::pp::DubinsCurves dubins;
  f2c::pp::ConstrainedTurn constrained;
  // Adjacent swaths (6 m apart) need a teardrop 13.94 m deep; the band is 12.
  // A rectangle has no corner to wrap, so the answer is honest failure.
  const F2CCells band {rectangle(-60.0, 0.0, 60.0, 12.0)};
  const auto res = constrained.createTurn(
      robot, band, F2CPoint(0.0, 0.0), M_PI_2,
      F2CPoint(6.0, 0.0), -M_PI_2, dubins);

  EXPECT_FALSE(res.inside);
  EXPECT_FALSE(res.used_waypoint);
  EXPECT_GT(res.length_outside, 0.05);
  EXPECT_NEAR(res.deepest_outside, 13.937 - 12.0, 0.05);
}

TEST(fields2cover_pp_constrained_turn, settings_round_trip) {
  f2c::pp::ConstrainedTurn turn;
  turn.setWaypointOffset(-1.5);
  turn.setSampleStep(-0.1);
  turn.setTolerance(-0.2);
  EXPECT_NEAR(turn.getWaypointOffset(), 1.5, 1e-9);
  EXPECT_NEAR(turn.getSampleStep(), 0.1, 1e-9);
  EXPECT_NEAR(turn.getTolerance(), 0.2, 1e-9);
}
