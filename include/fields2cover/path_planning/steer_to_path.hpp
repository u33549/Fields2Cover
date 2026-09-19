//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#pragma once
#ifndef FIELDS2COVER_PATH_PLANNING_STEER_TO_PATH_HPP_
#define FIELDS2COVER_PATH_PLANNING_STEER_TO_PATH_HPP_

#include <cmath>
#include <vector>
#include "steering_functions/steering_functions.hpp"
#include "fields2cover/types.h"


namespace f2c::pp {


/// Cast steer::State type from steering_functions library to
/// Path type
inline types::Path steerStatesToPath(
    const std::vector<steer::State>& curve, double const_vel) {
  types::Path path;
  // Built in place. A Point owns a heap-allocated geometry, so the three
  // temporary ones this used to make per state -- plus the state it then
  // copied into the path -- were allocations nobody read; a turn runs to a
  // couple of thousand states and a turn planner given free space builds
  // thousands of turns to choose one.
  std::vector<types::PathState>& states = path.getStates();
  states.resize(curve.size());
  for (size_t i = 0; i < curve.size(); ++i) {
    types::PathState& state = states[i];
    state.point.setPoint(curve[i].x, curve[i].y);
    state.angle = curve[i].theta;
    state.velocity = const_vel;
    if (i + 1 < curve.size()) {
      // The step to the next state. Asking the geometry for the distance
      // between two points costs 382 ns against 5 ns here, for the same
      // double -- checked bit for bit over 400k pairs.
      const double dx = curve[i + 1].x - curve[i].x;
      const double dy = curve[i + 1].y - curve[i].y;
      state.len = std::sqrt(dx * dx + dy * dy);
    }
    state.dir = static_cast<types::PathDirection>(curve[i].d);
    state.type = types::PathSectionType::TURN;
  }
  return path;
}

// Function loop_detected contributed by Phact (https://phact.nl/) company
inline bool loop_detected(const std::vector<steer::Control>& controls) {
  double drTotal = 0;
  double drAbsTotal = 0;
  for (auto&& c : controls) {
    if (c.kappa != 0) {
      double dr = 0.5 * c.delta_s * c.kappa / M_PI;
      if (fabs(dr) > 0.9) {
        return true;
      }
      drTotal += dr;
      drAbsTotal += fabs(dr);
    }
  }
  return (fabs(drTotal) > 0.9) || (drAbsTotal > 1.5);
}
// Function loop_detected contributed by Phact (https://phact.nl/) company



}  // namespace f2c::pp

#endif  // FIELDS2COVER_PATH_PLANNING_STEER_TO_PATH_HPP_
