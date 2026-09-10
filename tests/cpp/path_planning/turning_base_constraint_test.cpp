//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <gtest/gtest.h>
#include <cmath>
#include "fields2cover/path_planning/turning_base.h"
#include "fields2cover/path_planning/dubins_curves.h"
#include "fields2cover/path_planning/dubins_curves_cc.h"
#include "fields2cover/path_planning/reeds_shepp_curves.h"
#include "fields2cover/path_planning/reeds_shepp_curves_hc.h"
#include <memory>
#include <string>
#include <utility>
#include <vector>
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


// The rectangle with a square hole punched out of it. The hole ring is wound
// as asked for: GDAL hands a ring back the way it was given, so the two
// windings are a real pair of inputs, not the same one twice.
F2CCells holedRegion(bool hole_clockwise) {
  F2CLinearRing hole;
  const double c[4][2] {{-10, -10}, {10, -10}, {10, 10}, {-10, 10}};
  for (int i = 0; i < 5; ++i) {
    const int j = hole_clockwise ? (4 - i % 4) % 4 : i % 4;
    hole.addPoint(F2CPoint(c[j][0], c[j][1]));
  }
  F2CCell cell = rectangle(-50.0, -50.0, 50.0, 50.0);
  cell.addRing(hole);
  F2CCells region;
  region.addGeometry(cell);
  return region;
}

bool hasCorner(const std::vector<std::pair<F2CPoint, double>>& corners,
    const F2CPoint& p, double bisector) {
  for (auto&& c : corners) {
    if (c.first.distance(p) < 1e-6 &&
        F2CPoint::getAngleDiffAbs(c.second, bisector) < 1e-6) {
      return true;
    }
  }
  return false;
}

int countBackward(const F2CPath& path) {
  int n = 0;
  for (auto&& s : path) {
    if (s.dir == f2c::types::PathDirection::BACKWARD) { ++n; }
  }
  return n;
}

// One fresh instance of each planner: what the base class learned has to
// reach all four, not only the two it was written against.
std::vector<std::unique_ptr<f2c::pp::TurningBase>> allPlanners() {
  std::vector<std::unique_ptr<f2c::pp::TurningBase>> v;
  v.push_back(std::make_unique<f2c::pp::DubinsCurves>());
  v.push_back(std::make_unique<f2c::pp::DubinsCurvesCC>());
  v.push_back(std::make_unique<f2c::pp::ReedsSheppCurves>());
  v.push_back(std::make_unique<f2c::pp::ReedsSheppCurvesHC>());
  return v;
}

const char* const kPlannerNames[] {"DubinsCurves", "DubinsCurvesCC",
    "ReedsSheppCurves", "ReedsSheppCurvesHC"};

}  // namespace

TEST(fields2cover_pp_turn_constraint, corners_of_a_convex_region) {
  // A rectangle turns outward at every vertex, so it has nothing to wrap.
  const F2CCells box {rectangle(-50.0, 0.0, 50.0, 18.0)};
  EXPECT_TRUE(f2c::pp::TurningBase::concaveCorners(box).empty());
}

TEST(fields2cover_pp_turn_constraint, corners_of_a_wrapped_corner) {
  // The band around a right-angled crop corner turns inward exactly once.
  const auto corners = f2c::pp::TurningBase::concaveCorners(
      cornerBand(12.0, 80.0));
  ASSERT_EQ(corners.size(), 1);
  EXPECT_NEAR(corners[0].first.getX(), 0.0, 1e-6);
  EXPECT_NEAR(corners[0].first.getY(), 0.0, 1e-6);
  // The bisector points away from the crop, into the band.
  EXPECT_NEAR(F2CPoint::getAngleDiffAbs(corners[0].second, 5.0 * M_PI_4),
      0.0, 1e-6);
}

TEST(fields2cover_pp_turn_constraint, a_turn_that_fits_is_left_alone) {
  F2CRobot robot = makeRobot();
  f2c::pp::DubinsCurves dubins;
  // Two swath ends 18 m apart on one face of an 18 m band: a plain u-turn.
  const F2CCells band {rectangle(-60.0, 0.0, 60.0, 18.0)};
  const F2CPoint start(0.0, 0.0), end(18.0, 0.0);

  const F2CPath plain = dubins.createTurn(robot, start, M_PI_2, end, -M_PI_2);
  dubins.setFreeSpace(band);
  f2c::pp::TurnReport report;
  const F2CPath path = dubins.createTurn(
      robot, start, M_PI_2, end, -M_PI_2, &report);

  EXPECT_TRUE(report.inside);
  EXPECT_FALSE(report.used_waypoint);
  EXPECT_NEAR(path.length(), plain.length(), 1e-9);
}

TEST(fields2cover_pp_turn_constraint, a_turn_that_cuts_a_corner_is_repaired) {
  F2CRobot robot = makeRobot();
  f2c::pp::DubinsCurves dubins;
  const F2CCells band = cornerBand(12.0, 80.0);
  // One pose on each face of the crop corner, each perpendicular to its face:
  // out of the crop on the bottom, back into it on the left.
  const F2CPoint start(24.0, 0.0), end(0.0, 24.0);
  const double start_angle = -M_PI_2, end_angle = 0.0;

  const F2CPath plain = dubins.createTurn(
      robot, start, start_angle, end, end_angle);
  ASSERT_GT(lengthOutside(plain, band), 1.0)
      << "the plain turn is expected to cut across the crop";

  dubins.setFreeSpace(band);
  f2c::pp::TurnReport report;
  const F2CPath path = dubins.createTurn(
      robot, start, start_angle, end, end_angle, &report);

  EXPECT_TRUE(report.inside);
  EXPECT_TRUE(report.used_waypoint);
  EXPECT_LT(lengthOutside(path, band), 0.05);
  // Staying on the ground must not cost drivability: the radius still holds,
  // and the two halves meet at one pose, so the ends are the ones asked for.
  EXPECT_GT(minRadius(path), 0.95 * robot.getMinTurningRadius());
  EXPECT_NEAR(path.getStates().front().point.distance(start), 0.0, 1e-6);
  EXPECT_NEAR(path.getStates().back().point.distance(end), 0.0, 1e-6);
  EXPECT_GT(path.length(), plain.length());
}

TEST(fields2cover_pp_turn_constraint, no_waypoint_can_widen_a_narrow_band) {
  F2CRobot robot = makeRobot();
  f2c::pp::DubinsCurves dubins;
  // Adjacent swaths (6 m apart) need a teardrop 13.94 m deep; the band is 12.
  // A rectangle has no corner to wrap, so the answer is honest failure.
  const F2CCells band {rectangle(-60.0, 0.0, 60.0, 12.0)};
  dubins.setFreeSpace(band);
  f2c::pp::TurnReport report;
  dubins.createTurn(robot, F2CPoint(0.0, 0.0), M_PI_2,
      F2CPoint(6.0, 0.0), -M_PI_2, &report);

  EXPECT_FALSE(report.inside);
  EXPECT_FALSE(report.used_waypoint);
  EXPECT_GT(report.length_outside, 0.05);
  EXPECT_NEAR(report.deepest_outside, 13.937 - 12.0, 0.05);
}

TEST(fields2cover_pp_turn_constraint, inside_beats_shorter_across_planners) {
  F2CRobot robot = makeRobot();
  f2c::pp::ReedsSheppCurves reeds_shepp;
  // A u-turn between two ends of one face of an 18 m band. Reeds-Shepp reaches
  // it by reversing the whole way, which leaves the band; the forward turn is
  // the same length and stays in. Length cannot tell them apart, containment
  // can, and the forward turn is one Reeds-Shepp can drive itself.
  const F2CCells band {rectangle(-60.0, 0.0, 60.0, 18.0)};
  const F2CPoint start(0.0, 0.0), end(12.0, 0.0);
  const double start_angle = M_PI_2, end_angle = -M_PI_2;

  const F2CPath rs_alone = reeds_shepp.createTurn(
      robot, start, start_angle, end, end_angle);
  ASSERT_GT(lengthOutside(rs_alone, band), 1.0)
      << "Reeds-Shepp alone is expected to leave the band here";

  reeds_shepp.setFreeSpace(band);
  f2c::pp::TurnReport report;
  const F2CPath path = reeds_shepp.createTurn(
      robot, start, start_angle, end, end_angle, &report);

  EXPECT_TRUE(report.inside);
  EXPECT_FALSE(report.used_waypoint);
  EXPECT_LT(lengthOutside(path, band), 0.05);
  // Keeping the contained answer costs nothing in length here.
  EXPECT_NEAR(path.length(), rs_alone.length(), 0.5);
}

TEST(fields2cover_pp_turn_constraint, a_planner_that_already_fits_is_kept) {
  F2CRobot robot = makeRobot();
  f2c::pp::DubinsCurves dubins;
  // Given no ground to keep to, a planner answers exactly as it always did.
  const F2CCells band {rectangle(-60.0, 0.0, 60.0, 18.0)};
  const F2CPath before = dubins.createTurn(
      robot, F2CPoint(0.0, 0.0), M_PI_2, F2CPoint(24.0, 0.0), -M_PI_2);
  dubins.setFreeSpace(band);
  f2c::pp::TurnReport report;
  const F2CPath after = dubins.createTurn(robot,
      F2CPoint(0.0, 0.0), M_PI_2, F2CPoint(24.0, 0.0), -M_PI_2, &report);
  EXPECT_TRUE(report.inside);
  EXPECT_FALSE(report.used_waypoint);
  EXPECT_NEAR(after.length(), before.length(), 1e-9);
}

TEST(fields2cover_pp_turn_constraint, settings_round_trip) {
  f2c::pp::DubinsCurves turn;
  turn.setWaypointOffset(-1.5);
  turn.setSwathWidth(-6.0);
  EXPECT_NEAR(turn.getWaypointOffset(), 1.5, 1e-9);
  EXPECT_NEAR(turn.getSwathWidth(), 6.0, 1e-9);
  EXPECT_TRUE(turn.getFreeSpace().isEmpty());
  turn.setFreeSpace(F2CCells{rectangle(0.0, 0.0, 1.0, 1.0)});
  EXPECT_FALSE(turn.getFreeSpace().isEmpty());
}

TEST(fields2cover_pp_turn_constraint, every_planner_wraps_the_corner) {
  F2CRobot robot = makeRobot();
  const F2CCells band = cornerBand(12.0, 80.0);
  const F2CPoint start(24.0, 0.0), end(0.0, 24.0);
  const double start_angle = -M_PI_2, end_angle = 0.0;

  auto planners = allPlanners();
  for (size_t i = 0; i < planners.size(); ++i) {
    SCOPED_TRACE(kPlannerNames[i]);
    f2c::pp::TurningBase& planner = *planners[i];

    const F2CPath plain = planner.createTurn(
        robot, start, start_angle, end, end_angle);
    ASSERT_GT(lengthOutside(plain, band), 1.0)
        << "the plain turn is expected to cut across the crop";

    planner.setFreeSpace(band);
    f2c::pp::TurnReport report;
    const F2CPath path = planner.createTurn(
        robot, start, start_angle, end, end_angle, &report);

    EXPECT_TRUE(report.inside);
    EXPECT_TRUE(report.used_waypoint);
    EXPECT_LT(lengthOutside(path, band), 0.05);
    EXPECT_GT(minRadius(path), 0.95 * robot.getMinTurningRadius());
    EXPECT_NEAR(path.getStates().front().point.distance(start), 0.0, 1e-6);
    EXPECT_NEAR(path.getStates().back().point.distance(end), 0.0, 1e-6);
  }
}

TEST(fields2cover_pp_turn_constraint, no_free_space_leaves_planners_alone) {
  F2CRobot robot = makeRobot();
  // Short and long, aligned and skewed: whatever the turn, an empty free
  // space has to leave it untouched, point for point.
  const double cases[4][6] {
      {0.0, 0.0, M_PI_2, 12.0, 0.0, -M_PI_2},
      {0.0, 0.0, M_PI_2, 6.0, 0.0, -M_PI_2},
      {3.0, -2.0, 0.4, 27.0, 9.0, 2.7},
      {-5.0, 11.0, 3.0, 14.0, -7.0, -0.8}};

  auto planners = allPlanners();
  for (size_t i = 0; i < planners.size(); ++i) {
    f2c::pp::TurningBase& planner = *planners[i];
    for (size_t c = 0; c < 4; ++c) {
      SCOPED_TRACE(std::string(kPlannerNames[i]) + " case " +
          std::to_string(c));
      const F2CPoint start(cases[c][0], cases[c][1]);
      const F2CPoint end(cases[c][3], cases[c][4]);
      const F2CPath before = planner.createTurn(
          robot, start, cases[c][2], end, cases[c][5]);
      ASSERT_GT(before.size(), 1);

      // The other two settings are on as well: with no ground to keep to,
      // they have nothing to say either.
      planner.setFreeSpace(F2CCells());
      planner.setSwathWidth(6.0);
      planner.setWaypointOffset(2.0);
      f2c::pp::TurnReport report;
      const F2CPath after = planner.createTurn(
          robot, start, cases[c][2], end, cases[c][5], &report);

      ASSERT_EQ(after.size(), before.size());
      EXPECT_NEAR(after.length(), before.length(), 1e-12);
      for (size_t k = 0; k < before.size(); ++k) {
        ASSERT_NEAR(after[k].point.distance(before[k].point), 0.0, 1e-12);
        ASSERT_NEAR(after[k].angle, before[k].angle, 1e-12);
      }
      // A planner that was never told where the crop is does not get to
      // claim the turn ran over it.
      EXPECT_TRUE(report.inside);
      EXPECT_FALSE(report.used_waypoint);
      EXPECT_NEAR(report.length_outside, 0.0, 1e-12);
      EXPECT_NEAR(report.length_in_swath, 0.0, 1e-12);
    }
  }
}

TEST(fields2cover_pp_turn_constraint, a_report_is_optional) {
  F2CRobot robot = makeRobot();
  // Asking without a report is the default call, and it has to survive every
  // way createTurn can answer.
  auto answers_the_same = [&robot](f2c::pp::TurningBase& planner,
      const F2CPoint& start, double start_angle,
      const F2CPoint& end, double end_angle, f2c::pp::TurnReport* report) {
    const F2CPath with = planner.createTurn(
        robot, start, start_angle, end, end_angle, report);
    const F2CPath without = planner.createTurn(
        robot, start, start_angle, end, end_angle);
    EXPECT_EQ(without.size(), with.size());
    EXPECT_NEAR(without.length(), with.length(), 1e-9);
    EXPECT_GT(without.size(), 1);
  };

  const F2CCells band {rectangle(-60.0, 0.0, 60.0, 18.0)};
  f2c::pp::TurnReport report;

  // The plain turn fits.
  f2c::pp::DubinsCurves fits;
  fits.setFreeSpace(band);
  answers_the_same(fits, F2CPoint(0.0, 0.0), M_PI_2,
      F2CPoint(18.0, 0.0), -M_PI_2, &report);
  EXPECT_TRUE(report.inside);
  EXPECT_FALSE(report.used_waypoint);

  // Another of the planner's own turns fits.
  f2c::pp::ReedsSheppCurves alternative;
  alternative.setFreeSpace(band);
  answers_the_same(alternative, F2CPoint(0.0, 0.0), M_PI_2,
      F2CPoint(12.0, 0.0), -M_PI_2, &report);
  EXPECT_TRUE(report.inside);
  EXPECT_FALSE(report.used_waypoint);

  // A waypoint is needed.
  f2c::pp::DubinsCurves wrapped;
  wrapped.setFreeSpace(cornerBand(12.0, 80.0));
  answers_the_same(wrapped, F2CPoint(24.0, 0.0), -M_PI_2,
      F2CPoint(0.0, 24.0), 0.0, &report);
  EXPECT_TRUE(report.used_waypoint);

  // Nothing fits.
  f2c::pp::DubinsCurves cornered;
  cornered.setFreeSpace(F2CCells{rectangle(-60.0, 0.0, 60.0, 12.0)});
  answers_the_same(cornered, F2CPoint(0.0, 0.0), M_PI_2,
      F2CPoint(6.0, 0.0), -M_PI_2, &report);
  EXPECT_FALSE(report.inside);
}

TEST(fields2cover_pp_turn_constraint, a_planner_with_one_answer_has_no_other) {
  F2CRobot robot = makeRobot();
  f2c::pp::DubinsCurves dubins;
  // Dubins drives one shortest turn and nothing else; the base class must not
  // invent a second for it.
  EXPECT_TRUE(dubins.alternativeTurns(robot, F2CPoint(0.0, 0.0), M_PI_2,
      F2CPoint(12.0, 0.0), -M_PI_2).empty());
}

TEST(fields2cover_pp_turn_constraint, reeds_shepp_offers_its_forward_turn) {
  F2CRobot robot = makeRobot();
  const F2CPoint start(0.0, 0.0), end(12.0, 0.0);
  const double start_angle = M_PI_2, end_angle = -M_PI_2;
  f2c::pp::DubinsCurves dubins;
  const F2CPath forward = dubins.createTurn(
      robot, start, start_angle, end, end_angle);

  f2c::pp::ReedsSheppCurves reeds_shepp;
  f2c::pp::ReedsSheppCurvesHC reeds_shepp_hc;
  std::vector<f2c::pp::TurningBase*> planners {&reeds_shepp, &reeds_shepp_hc};
  for (size_t i = 0; i < planners.size(); ++i) {
    SCOPED_TRACE(i == 0 ? "ReedsSheppCurves" : "ReedsSheppCurvesHC");
    const auto alternatives = planners[i]->alternativeTurns(
        robot, start, start_angle, end, end_angle);
    ASSERT_EQ(alternatives.size(), 1);
    EXPECT_NEAR(alternatives[0].length(), forward.length(), 1e-9);
    // Between these two poses Reeds-Shepp reverses the whole way; what it
    // offers instead is the same turn driven forwards.
    EXPECT_EQ(countBackward(alternatives[0]), 0);
    EXPECT_GT(countBackward(planners[i]->createTurn(
        robot, start, start_angle, end, end_angle)), 0);
  }
}

TEST(fields2cover_pp_turn_constraint, corners_of_a_region_with_a_hole) {
  const auto corners = f2c::pp::TurningBase::concaveCorners(
      holedRegion(false));
  // The outer ring turns outward everywhere; the hole is what the ground has
  // to be driven around, and each of its vertices turns inward.
  ASSERT_EQ(corners.size(), 4);
  EXPECT_TRUE(hasCorner(corners, F2CPoint(10.0, 10.0), M_PI_4));
  EXPECT_TRUE(hasCorner(corners, F2CPoint(-10.0, 10.0), 3.0 * M_PI_4));
  EXPECT_TRUE(hasCorner(corners, F2CPoint(-10.0, -10.0), -3.0 * M_PI_4));
  EXPECT_TRUE(hasCorner(corners, F2CPoint(10.0, -10.0), -M_PI_4));
}

TEST(fields2cover_pp_turn_constraint, corners_do_not_depend_on_the_winding) {
  const auto ccw = f2c::pp::TurningBase::concaveCorners(holedRegion(false));
  const auto cw = f2c::pp::TurningBase::concaveCorners(holedRegion(true));
  // The same hole wound the other way is the same hole. Reading the winding
  // would flip every verdict; asking the region which side it is on does not.
  ASSERT_EQ(cw.size(), ccw.size());
  ASSERT_EQ(ccw.size(), 4);
  for (auto&& c : ccw) {
    EXPECT_TRUE(hasCorner(cw, c.first, c.second))
        << "reversed ring lost the corner at ("
        << c.first.getX() << ", " << c.first.getY() << ")";
  }
}

TEST(fields2cover_pp_turn_constraint, ground_a_swath_covers_is_not_crop) {
  F2CRobot robot = makeRobot();
  // The drivable band starts a metre above the swath ends, so the turn's
  // first and last metre are off it -- and lie straight along the swaths it
  // joins, which get driven either way.
  const F2CCells band {rectangle(-60.0, 1.0, 60.0, 30.0)};
  const F2CPoint start(0.0, 0.0), end(12.0, 0.0);

  f2c::pp::DubinsCurves strict;
  strict.setFreeSpace(band);
  f2c::pp::TurnReport strict_rep;
  const F2CPath strict_path = strict.createTurn(
      robot, start, M_PI_2, end, -M_PI_2, &strict_rep);
  EXPECT_FALSE(strict_rep.inside);
  EXPECT_NEAR(strict_rep.length_outside, 2.0, 0.1);
  EXPECT_NEAR(strict_rep.length_in_swath, 0.0, 1e-9);
  EXPECT_NEAR(strict_rep.deepest_outside, 1.0, 0.05);

  f2c::pp::DubinsCurves aware;
  aware.setFreeSpace(band);
  aware.setSwathWidth(6.0);
  f2c::pp::TurnReport aware_rep;
  const F2CPath aware_path = aware.createTurn(
      robot, start, M_PI_2, end, -M_PI_2, &aware_rep);
  EXPECT_TRUE(aware_rep.inside);
  EXPECT_NEAR(aware_rep.length_outside, 0.0, 1e-9);
  EXPECT_NEAR(aware_rep.length_in_swath, 2.0, 0.1);
  EXPECT_NEAR(aware_rep.deepest_outside, 0.0, 1e-9);
  // Same turn either way: what changed is the verdict on it.
  EXPECT_NEAR(aware_path.length(), strict_path.length(), 1e-9);
}

TEST(fields2cover_pp_turn_constraint, a_waypoint_leaves_no_seam) {
  F2CRobot robot = makeRobot();
  f2c::pp::DubinsCurves dubins;
  dubins.setFreeSpace(cornerBand(12.0, 80.0));
  f2c::pp::TurnReport report;
  const F2CPath path = dubins.createTurn(
      robot, F2CPoint(24.0, 0.0), -M_PI_2, F2CPoint(0.0, 24.0), 0.0, &report);
  ASSERT_TRUE(report.used_waypoint);

  double max_step = 0.0, max_swing = 0.0, to_waypoint = 1e18;
  const auto& states = path.getStates();
  for (size_t i = 0; i + 1 < states.size(); ++i) {
    max_step = std::max(max_step,
        states[i].point.distance(states[i + 1].point));
    max_swing = std::max(max_swing, F2CPoint::getAngleDiffAbs(
        states[i].angle, states[i + 1].angle));
    to_waypoint = std::min(to_waypoint,
        states[i].point.distance(report.waypoint));
  }
  // Both halves are planned to the same pose, so the join is a step like any
  // other along the turn: the vehicle neither jumps nor spins there.
  EXPECT_LT(to_waypoint, 1e-6);
  EXPECT_LT(max_step, 1.5 * dubins.getDiscretization());
  EXPECT_LT(max_swing,
      2.0 * dubins.getDiscretization() / robot.getMinTurningRadius());
}

TEST(fields2cover_pp_turn_constraint, a_long_run_down_a_swath_is_still_crop) {
  F2CRobot robot = makeRobot();
  const double R = robot.getMinTurningRadius();
  // Both poses on one swath line, thirty metres apart, facing the same way:
  // Reeds-Shepp answers by reversing the whole way down the line. Perfectly
  // aligned with the swath, and still thirty metres driven over it.
  const F2CCells elsewhere {rectangle(-60.0, 5.0, 60.0, 40.0)};
  const F2CPoint start(0.0, 0.0), end(0.0, -30.0);

  f2c::pp::ReedsSheppCurves planner;
  planner.setFreeSpace(elsewhere);
  planner.setSwathWidth(6.0);
  f2c::pp::TurnReport report;
  const F2CPath path = planner.createTurn(
      robot, start, M_PI_2, end, M_PI_2, &report);
  ASSERT_NEAR(path.length(), 30.0, 1e-6);

  // Settling onto a swath costs about a turning radius at each end; past
  // that the alignment excuses nothing. Excusing the whole line would report
  // this turn as clean.
  EXPECT_FALSE(report.inside);
  EXPECT_NEAR(report.length_in_swath, 2.0 * R, 0.1);
  EXPECT_NEAR(report.length_outside, 30.0 - 2.0 * R, 0.1);
}
