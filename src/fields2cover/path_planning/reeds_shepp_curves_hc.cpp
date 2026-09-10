//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <steering_functions/hc_cc_state_space/hc00_reeds_shepp_state_space.hpp>
#include "fields2cover/path_planning/steer_to_path.hpp"
#include "fields2cover/path_planning/reeds_shepp_curves_hc.h"

namespace f2c::pp {

F2CPath ReedsSheppCurvesHC::createSimpleTurn(const F2CRobot& robot,
    double dist_start_pos, double start_angle, double end_angle) {
  steer::State start, end;

  start.x = 0.0;
  start.y = 0.0;
  start.theta = start_angle;
  start.kappa = 0.0;
  start.d = 0;

  end.x = dist_start_pos;
  end.y = 0.0;
  end.theta = end_angle;
  end.kappa = 0.0;
  end.d = 0;

  HC00_Reeds_Shepp_State_Space ss(
      robot.getMaxCurv(),
      robot.getMaxDiffCurv(),
      discretization);

  return steerStatesToPath(ss.get_path(start, end),
      robot.getTurnVel());
}

bool ReedsSheppCurvesHC::hasContinuousCurvature() const {
  return true;
}

std::vector<F2CPath> ReedsSheppCurvesHC::alternativeTurns(const F2CRobot& robot,
    const F2CPoint& start_pos, double start_angle,
    const F2CPoint& end_pos, double end_angle) {
  forward_.setDiscretization(this->getDiscretization());
  forward_.setUsingCache(this->getUsingCache());
  const F2CPath path = forward_.createTurn(
      robot, start_pos, start_angle, end_pos, end_angle);
  if (path.size() < 2) {
    return {};
  }
  return {path};
}

}  // namespace f2c::pp

