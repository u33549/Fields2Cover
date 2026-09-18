//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#pragma once
#ifndef FIELDS2COVER_PATH_PLANNING_REEDS_SHEPP_CURVES_HC_H_
#define FIELDS2COVER_PATH_PLANNING_REEDS_SHEPP_CURVES_HC_H_

#include "fields2cover/types.h"
#include <vector>
#include "fields2cover/path_planning/dubins_curves.h"
#include "fields2cover/path_planning/turning_base.h"

namespace f2c::pp {

/// Reeds-Shepp's curves planner with continuous curves
class ReedsSheppCurvesHC : public TurningBase {
 public:
  F2CPath createSimpleTurn(const F2CRobot& robot,
      double dist_start_pos, double start_angle, double end_angle) override;
  bool hasContinuousCurvature() const override;

  /// @brief The same turn driven forwards only.
  /// @details Reeds-Shepp is Dubins plus reverse, so the forward turn is one
  /// of its own answers, just not the shortest one. It matters when the
  /// shortest reverses over the crop and the forward one does not: told where
  /// it may drive, this planner can then keep the one that stays there.
  std::vector<F2CPath> alternativeTurns(const F2CRobot& robot,
      const F2CPoint& start_pos, double start_angle,
      const F2CPoint& end_pos, double end_angle) override;

 private:
  DubinsCurves forward_;
};

}  // namespace f2c::pp

#endif  // FIELDS2COVER_PATH_PLANNING_REEDS_SHEPP_CURVES_HC_H_
