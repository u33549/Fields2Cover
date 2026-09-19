//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#pragma once
#ifndef FIELDS2COVER_PATH_PLANNING_TURNING_BASE_H_
#define FIELDS2COVER_PATH_PLANNING_TURNING_BASE_H_

#include <memory>
#include <map>
#include <vector>
#include <limits>
#include <functional>
#include "fields2cover/types.h"
#include "fields2cover/utils/prepared_area.h"
#include "fields2cover/utils/random.h"

namespace f2c::pp {

/// @brief What a turn had to do to stay on drivable ground.
/// @details Only filled in when the planner has been given a free space; a
/// turn planner without one answers in free space and cannot know.
struct TurnReport {
  /// The turn stays on drivable ground.
  bool inside {true};
  /// The turn had to be routed through a waypoint to stay there.
  bool used_waypoint {false};
  /// Length driven over the crop [m].
  double length_outside {0.0};
  /// Length outside the drivable ground but within a swath's own strip [m].
  double length_in_swath {0.0};
  /// Deepest the turn reaches into the crop [m].
  double deepest_outside {0.0};
  /// The waypoint used, when one was.
  F2CPoint waypoint;
  /// The turn also stays on the ground the route travels through.
  bool in_preferred {true};
  /// Length driven off that ground but still clear of the crop [m].
  double length_off_preferred {0.0};
};

/// Base class for turn planners
class TurningBase {
 public:
  /// @brief Create a turn that goes from one point with a certain angle to
  /// another point.
  /// @details Start and end point are connected with a line that creates the
  /// inferior border of the turn.
  /// @param start_pos Start point
  /// @param start_angle Start angle
  /// @param end_pos End point
  /// @param end_angle End angle
  /// @return Path with the computed turn
  /// @param report Filled in with what the turn had to do, when asked for
  F2CPath createTurn(const F2CRobot& robot,
      const F2CPoint& start_pos, double start_angle,
      const F2CPoint& end_pos, double end_angle,
      TurnReport* report = nullptr);

  /// @brief Other turns this planner could drive between the same two poses.
  /// @details Asked for the shortest, a planner answers with one path. Asked
  /// to stay on drivable ground, it may need a longer one it can also drive:
  /// Reeds-Shepp reverses through the crop when the forward turn costs the
  /// same, and the forward turn is one of its own. Planners with nothing else
  /// to offer return nothing, and then the shortest is all there is.
  virtual std::vector<F2CPath> alternativeTurns(const F2CRobot& robot,
      const F2CPoint& start_pos, double start_angle,
      const F2CPoint& end_pos, double end_angle);

  /// @brief Concave (reflex) corners of a region, with the bisector pointing
  /// into it.
  /// @details A turn that leaves drivable ground leaves it at one of these:
  /// no three-part curve wraps a corner. A region without any cannot be helped
  /// by a waypoint, and the turn simply does not fit.
  static std::vector<std::pair<F2CPoint, double>> concaveCorners(
      const F2CCells& region, double eps = 1e-3);

  /// @brief Generate a turn if it has not been computed before.
  /// @param dist_start_pos Distance between start and end point
  /// @param start_angle Angle when going into the headland
  /// (0 deg is the angle of the headland)
  /// @param end_angle Angle when going out of the headland
  F2CPath createTurnIfNotCached(const F2CRobot& robot, double dist_start_pos,
      double start_angle, double end_angle);

  /// @brief Create a turn.
  /// @param dist_start_pos Distance between start and end point
  /// @param start_angle Angle when going into the headland
  /// (0 deg is the angle of the headland)
  /// @param end_angle Angle when going out of the headland
  virtual F2CPath createSimpleTurn(const F2CRobot& robot, double dist_start_pos,
      double start_angle, double end_angle) = 0;

  /// Check if the turns keep the curvature continuous.
  virtual bool hasContinuousCurvature() const;

  /// @brief Transform the turn parameters representation from two points with
  /// two angles to one distance and two angles.
  /// @param start_pos Start point
  /// @param start_angle Start angle
  /// @param end_pos End point
  /// @param end_angle End angle
  /// @return Vector with the values of the new representation
  static std::vector<double> transformToNormalTurn(const F2CPoint& start_pos,
      double start_angle, const F2CPoint& end_pos, double end_angle);

  /// Check if turn is valid
  static bool isTurnValid(const F2CPath& path, double dist_start_end,
      double end_angle, double max_dist_error = 0.05,
      double max_rot_error = 0.1);


  /// Get discretization distance from points in the turn
  double getDiscretization() const;
  /// Set discretization distance from points in the turn
  void setDiscretization(double d);

  /// Get the ground the turns may be driven on. Empty means anywhere.
  const F2CCells& getFreeSpace() const;
  /// @brief Set the ground the turns may be driven on.
  /// @details A turn planner answers in free space: it is optimal there and
  /// keeps its radius, but it does not know where the crop is, and it leaves
  /// the drivable ground exactly where that ground turns inward. Given the
  /// ground, it keeps the shortest turn that stays on it, routing through a
  /// waypoint at a concave corner when the plain turn will not. Left empty,
  /// nothing changes.
  void setFreeSpace(const F2CCells& free_space);

  /// Get the ground the turns should keep to when they can. Empty means the
  /// whole free space is equally good.
  const F2CCells& getPreferredSpace() const;
  /// @brief Set the ground the turns should keep to when they can.
  /// @details Not all drivable ground is equally welcome. A route travels
  /// down the middle of the headland, so a turn that swings out to the field
  /// edge is clear of the crop but nowhere the route planned to be. This is
  /// a preference, not a limit: among the turns that stay off the crop the
  /// one that also stays here is kept, and if none does the shortest that
  /// stays off the crop is. Making it the free space instead would rank the
  /// field edge and the crop as equally forbidden, and a turn with nowhere
  /// to go would then choose the crop. Left empty, nothing changes.
  void setPreferredSpace(const F2CCells& preferred);

  /// Get the width of the swaths the turns join [m].
  double getSwathWidth() const;
  /// @brief Set the width of the swaths the turns join [m].
  /// @details Ground within half a swath of the line the vehicle came in on,
  /// or of the one it leaves on, is ground those swaths cover anyway, so a
  /// turn that clips it early is not driving over the crop. Counting it as
  /// such rejects turns that are fine. Left at zero, every metre off the
  /// drivable ground counts.
  void setSwathWidth(double width);

  /// Get how far inside a corner a waypoint sits, in turning radii.
  double getWaypointOffset() const;
  /// Set how far inside a corner a waypoint sits, in turning radii.
  void setWaypointOffset(double offset);

  /// Get if turns are being cached or not.
  bool getUsingCache() const;
  /// Set if cache should be used when planning same turn as before.
  void setUsingCache(bool c);

  virtual ~TurningBase() = default;


 private:
  /// The turn this planner would drive with nothing in the way.
  F2CPath plainTurn(const F2CRobot& robot,
      const F2CPoint& start_pos, double start_angle,
      const F2CPoint& end_pos, double end_angle);

  static void correctPath(F2CPath& path,
      const F2CPoint& start_pos,
      const F2CPoint& end_pos,
      float max_error_dist = 0.05);

 protected:
  // To prevent memory consumption and comparative errors because of doubles
  // ints are used multiplied by 1000.
  std::map<std::vector<int>, F2CPath> path_cache_;
  double discretization {0.01};
  bool using_cache {false};
  F2CCells free_space_;
  /// The same ground, read once so the sampling below can ask it cheaply.
  f2c::PreparedArea free_area_;
  /// Where a turn may leave the ground, worked out once per free space.
  std::vector<std::pair<F2CPoint, double>> free_corners_;
  F2CCells preferred_space_;
  f2c::PreparedArea preferred_area_;
  double swath_width_ {0.0};
  double waypoint_offset_ {1.0};
};

}  // namespace f2c::pp

#endif  // FIELDS2COVER_PATH_PLANNING_TURNING_BASE_H_
