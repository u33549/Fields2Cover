#==============================================================================
#     Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
#                      Author: Gonzalo Mier
#                         BSD-3 License
#==============================================================================

import pytest
import fields2cover as f2c
import math


def makeRobot():
  robot = f2c.Robot(3.0, 6.0);
  robot.setCruiseVel(2.0);
  robot.setMaxCurv(1.0 / 6.0);
  return robot;

def rectangle(x0, y0, x1, y1):
  ring = f2c.LinearRing();
  for x, y in [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]:
    ring.addPoint(f2c.Point(x, y));
  return f2c.Cell(ring);

def cornerBand(width, arm):
  band = f2c.Cells();
  band.addGeometry(rectangle(-width, -width, arm, 0.0));
  band.addGeometry(rectangle(-width, -width, 0.0, arm));
  return band.unionCascaded();

def allPlanners():
  return [f2c.PP_DubinsCurves(), f2c.PP_DubinsCurvesCC(),
      f2c.PP_ReedsSheppCurves(), f2c.PP_ReedsSheppCurvesHC()];


def test_fields2cover_pp_constraint_no_free_space_changes_nothing():
  robot = makeRobot();
  start = f2c.Point(0.0, 0.0);
  end = f2c.Point(12.0, 0.0);
  for turn in allPlanners():
    before = turn.createTurn(robot, start, 0.5 * math.pi, end, -0.5 * math.pi);
    turn.setFreeSpace(f2c.Cells());
    report = f2c.TurnReport();
    after = turn.createTurn(
        robot, start, 0.5 * math.pi, end, -0.5 * math.pi, report);
    assert after.size() == before.size();
    assert abs(after.length() - before.length()) < 1e-12;
    # A planner never told where the crop is cannot claim to have driven on it.
    assert report.inside;
    assert not report.used_waypoint;
    assert report.length_outside == 0.0;


def test_fields2cover_pp_constraint_corner_is_wrapped():
  robot = makeRobot();
  band = cornerBand(12.0, 80.0);
  start = f2c.Point(24.0, 0.0);
  end = f2c.Point(0.0, 24.0);

  turn = f2c.PP_DubinsCurves();
  plain = turn.createTurn(robot, start, -0.5 * math.pi, end, 0.0);
  assert any([s.point for s in plain.getStates() if not band.isPointIn(s.point)]), \
      "the plain turn is expected to cut across the crop";

  turn.setFreeSpace(band);
  report = f2c.TurnReport();
  path = turn.createTurn(robot, start, -0.5 * math.pi, end, 0.0, report);
  assert report.inside;
  assert report.used_waypoint;
  assert band.isPointIn(report.waypoint);
  assert path.length() > plain.length();
  assert min([s.point.distance(report.waypoint) for s in path.getStates()]) < 1e-6;


def test_fields2cover_pp_constraint_concave_corners():
  # A rectangle turns outward at every vertex; the band around a crop corner
  # turns inward exactly once, and points away from the crop.
  assert len(f2c.PP_Turning_base_class_concaveCorners(
      f2c.Cells(rectangle(-50.0, 0.0, 50.0, 18.0)))) == 0;
  corners = f2c.PP_Turning_base_class_concaveCorners(cornerBand(12.0, 80.0));
  assert len(corners) == 1;
  assert corners[0][0].distance(f2c.Point(0.0, 0.0)) < 1e-6;
  assert f2c.Point.getAngleDiffAbs(corners[0][1], 1.25 * math.pi) < 1e-6;


def test_fields2cover_pp_constraint_alternative_turns():
  robot = makeRobot();
  start = f2c.Point(0.0, 0.0);
  end = f2c.Point(12.0, 0.0);
  # Dubins drives one shortest turn and nothing else. Reeds-Shepp reaches
  # these poses by reversing, and offers the forward turn as well.
  assert len(f2c.PP_DubinsCurves().alternativeTurns(
      robot, start, 0.5 * math.pi, end, -0.5 * math.pi)) == 0;
  forward = f2c.PP_DubinsCurves().createTurn(
      robot, start, 0.5 * math.pi, end, -0.5 * math.pi);
  for turn in [f2c.PP_ReedsSheppCurves(), f2c.PP_ReedsSheppCurvesHC()]:
    alternatives = turn.alternativeTurns(
        robot, start, 0.5 * math.pi, end, -0.5 * math.pi);
    assert len(alternatives) == 1;
    assert abs(alternatives[0].length() - forward.length()) < 1e-9;


def test_fields2cover_pp_constraint_swath_width():
  robot = makeRobot();
  # The band starts a metre above the swath ends, so the turn's first and last
  # metre are off it, straight along the swaths it joins.
  band = f2c.Cells(rectangle(-60.0, 1.0, 60.0, 30.0));
  start = f2c.Point(0.0, 0.0);
  end = f2c.Point(12.0, 0.0);

  strict = f2c.PP_DubinsCurves();
  strict.setFreeSpace(band);
  strict_report = f2c.TurnReport();
  strict.createTurn(robot, start, 0.5 * math.pi, end, -0.5 * math.pi,
      strict_report);
  assert not strict_report.inside;
  assert abs(strict_report.length_outside - 2.0) < 0.1;

  aware = f2c.PP_DubinsCurves();
  aware.setFreeSpace(band);
  aware.setSwathWidth(6.0);
  assert aware.getSwathWidth() == 6.0;
  aware_report = f2c.TurnReport();
  aware.createTurn(robot, start, 0.5 * math.pi, end, -0.5 * math.pi,
      aware_report);
  # Same turn, different verdict: that ground is covered by the swaths anyway.
  assert aware_report.inside;
  assert aware_report.length_outside == 0.0;
  assert abs(aware_report.length_in_swath - 2.0) < 0.1;
