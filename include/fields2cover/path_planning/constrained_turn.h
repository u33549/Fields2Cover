//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#pragma once
#ifndef FIELDS2COVER_PATH_PLANNING_CONSTRAINED_TURN_H_
#define FIELDS2COVER_PATH_PLANNING_CONSTRAINED_TURN_H_

#include <vector>
#include "fields2cover/types.h"
#include "fields2cover/path_planning/turning_base.h"

namespace f2c::pp {

/// @brief Turn between two poses without leaving a given free space.
/// @details A turn planner answers in free space: Dubins is optimal there and
/// stays at the minimum radius, but it does not know where the crop is. Its
/// answer leaves the drivable ground exactly where that ground turns inward,
/// around a concave corner, and no three-part curve wraps a corner.
///
/// So when the plain turn leaves the free space, this retries through one
/// waypoint, placed a turning radius inside each concave corner and headed
/// along it, and keeps the shortest attempt that stays in. Both halves are
/// ordinary turns and they meet at the same pose, so the seam is continuous by
/// construction and the curvature limit is the turn planner's to keep.
class ConstrainedTurn {
 public:
  /// @brief What a constrained turn came back with.
  struct Result {
    /// Path from start to end pose. Empty only if no turn exists at all.
    F2CPath path;
    /// The path stays within the free space.
    bool inside {false};
    /// A waypoint was needed; the plain turn left the free space.
    bool used_waypoint {false};
    /// Deepest distance the path reaches outside the free space [m].
    double deepest_outside {0.0};
    /// Length of path outside the free space that the swaths do not cover [m].
    double length_outside {0.0};
    /// Length outside the free space but within a swath's own strip [m].
    double length_in_swath {0.0};
    /// Index of the turn planner whose answer was kept.
    size_t planner {0};
    /// The waypoint the turn was routed through, when one was needed.
    F2CPoint waypoint;
    /// Heading at that waypoint [rad].
    double waypoint_angle {0.0};
  };

  /// @brief Turn from one pose to another, staying inside free_space.
  /// @param robot Robot
  /// @param free_space Ground the vehicle may drive on
  /// @param start_pos Start point
  /// @param start_angle Start angle
  /// @param end_pos End point
  /// @param end_angle End angle
  /// @param turn Turn planner used for each half
  /// @return The turn, and whether it stays inside
  Result createTurn(const F2CRobot& robot, const F2CCells& free_space,
      const F2CPoint& start_pos, double start_angle,
      const F2CPoint& end_pos, double end_angle, TurningBase& turn) const;

  /// @brief Turn from one pose to another, staying inside free_space, choosing
  /// between several turn planners.
  /// @details Asked of one planner, a turn is whatever that planner thinks is
  /// shortest, and shortest is not the same as inside: Reeds-Shepp will drive a
  /// u-turn backwards through the crop when the forward turn costs the same
  /// length. Given several planners, the answer that stays inside wins, and
  /// only then the shorter one. Nothing is penalised for reversing — reversing
  /// buys room where the room is needed, and the containment test already
  /// rejects it where it is not.
  /// @param robot Robot
  /// @param free_space Ground the vehicle may drive on
  /// @param start_pos Start point
  /// @param start_angle Start angle
  /// @param end_pos End point
  /// @param end_angle End angle
  /// @param turns Turn planners to choose between, most preferred first
  /// @return The turn, whether it stays inside, and which planner gave it
  Result createTurn(const F2CRobot& robot, const F2CCells& free_space,
      const F2CPoint& start_pos, double start_angle,
      const F2CPoint& end_pos, double end_angle,
      const std::vector<TurningBase*>& turns) const;

  /// @brief Concave (reflex) corners of a region, with the bisector that
  /// points into the region.
  /// @details Public because the corners are the diagnosis: a turn that leaves
  /// the free space leaves it at one of these, and a region without any cannot
  /// be helped by a waypoint.
  /// @param free_space Region to inspect
  /// @return Corner points paired with the inward bisector angle [rad]
  static std::vector<std::pair<F2CPoint, double>> concaveCorners(
      const F2CCells& free_space, double eps = 1e-3);

  /// Get how far inside a corner the waypoint sits, in turning radii.
  double getWaypointOffset() const;
  /// Set how far inside a corner the waypoint sits, in turning radii.
  void setWaypointOffset(double offset);

  /// Get the step used to sample a path against the free space [m].
  double getSampleStep() const;
  /// Set the step used to sample a path against the free space [m].
  void setSampleStep(double step);

  /// Get the width of the swaths the turn joins [m]. Zero counts every metre
  /// outside the free space against the turn.
  double getSwathWidth() const;
  /// @brief Set the width of the swaths the turn joins [m].
  /// @details Ground within half a swath of the line the vehicle came in on,
  /// or of the one it leaves on, is ground those swaths cover anyway. A turn
  /// that clips it early is not driving over the crop, and counting it as
  /// such rejects turns that are perfectly fine — measured on the connection
  /// bench, most of what Reeds-Shepp appears to trample is this and nothing
  /// else. Left at zero, every metre outside counts.
  void setSwathWidth(double width);

  /// Get the length outside the free space still counted as inside [m].
  double getTolerance() const;
  /// Set the length outside the free space still counted as inside [m].
  void setTolerance(double tol);

 private:
  double waypoint_offset_ {1.0};
  double swath_width_ {0.0};
  double sample_step_ {0.25};
  double tolerance_ {0.05};
};

}  // namespace f2c::pp

#endif  // FIELDS2COVER_PATH_PLANNING_CONSTRAINED_TURN_H_
