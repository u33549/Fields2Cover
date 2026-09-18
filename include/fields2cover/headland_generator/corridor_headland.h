//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#pragma once
#ifndef FIELDS2COVER_HEADLAND_GENERATOR_CORRIDOR_HEADLAND_H_
#define FIELDS2COVER_HEADLAND_GENERATOR_CORRIDOR_HEADLAND_H_

#include <vector>
#include "fields2cover/types.h"
#include "fields2cover/headland_generator/headland_generator_base.h"
#include "fields2cover/path_planning/turning_base.h"
#include "fields2cover/swath_generator/swath_generator_base.h"
#include "fields2cover/objectives/sg_obj/sg_objective.h"

namespace f2c::hg {

/// How CorridorHL splits the corridor between two cells that share a border.
enum class CorridorShareMode {
  /// The smaller cell gives the whole corridor; ties split evenly. Default.
  ASYMMETRIC,
  /// Every shared border splits evenly, regardless of the cells' size.
  SYMMETRIC,
};

/// One shared border between two cells, and how the corridor over it is split.
///
/// A pair of touching cells gives two of these, one seen from each cell. The
/// pair that shares the border evenly gives 0.5 twice; otherwise the smaller
/// cell gives the whole corridor and the larger one keeps a row with a share
/// of 0, so a report can tell which side of the border the corridor came from.
struct CorridorShare {
  /// Cell the corridor is taken out of.
  size_t cell_i {0};
  /// Neighbour on the other side of the border.
  size_t cell_k {0};
  /// Perimeter of cell \a cell_i, the size the rule compares.
  double perimeter_i {0.0};
  /// Perimeter of cell \a cell_k.
  double perimeter_k {0.0};
  /// True when both perimeters are equal up to the tolerance of the rule.
  bool same_size {false};
  /// Part of the corridor width taken out of \a cell_i: 0, 0.5 or 1.
  double share {0.0};
  /// Length of the border the two cells share.
  double shared_length {0.0};
  /// Border segments the corridor is opened along.
  F2CMultiLineString shared_border;
};

/// Class to open a corridor where cells border each other, leaving the edges
/// that face the outer boundary or a void untouched.
///
/// Meant for cells that came out of a decomposition and already have a
/// headland around the field: shrinking every border again would take a second
/// headland off ground that has one.
class CorridorHL : public HeadlandGeneratorBase {
 public:
  using HeadlandGeneratorBase::generateHeadlands;

  /// Open a corridor of the given width between cells that share a border.
  /// @param field Cells that share borders, usually from a decomposition.
  /// @param dist_headland Width of the corridor.
  /// @return Mainland area
  F2CCells generateHeadlands(
    const F2CCells& field, double dist_headland) override;

  /// Split the corridor between each pair of cells that share a border.
  ///
  /// This is the rule generateHeadlands() applies, on its own: which cell
  /// gives the corridor, how much of it, and over which part of the border.
  /// @param field Cells that share borders, usually from a decomposition.
  /// @param mode Rule to split the corridor with.
  /// @return One share per ordered pair of cells that touch along a border.
  std::vector<CorridorShare> corridorShares(
      const F2CCells& field,
      CorridorShareMode mode = CorridorShareMode::ASYMMETRIC) const;

  /// Rule generateHeadlands() uses to split the corridor. Defaults to
  /// ASYMMETRIC.
  void setShareMode(CorridorShareMode mode);

  /// Rule generateHeadlands() currently uses to split the corridor.
  CorridorShareMode getShareMode() const;

  /// Open a corridor as wide as the turn the planner actually makes.
  ///
  /// The corridor is where a turn at the end of a swath is driven, so it has
  /// to be as deep as that turn reaches. Twice the turning radius is only a
  /// bound: with the swaths far enough apart the turn reaches one radius out,
  /// and with them closer together than twice the radius it has to loop and
  /// reaches further than the bound allows for.
  ///
  /// A piece left narrower than the robot's coverage width is given to the
  /// corridor instead of returned: nothing can cover it without the implement
  /// hanging over the corridor anyway.
  /// @param field Cells that share borders, usually from a decomposition.
  /// @param robot Robot doing the coverage.
  /// @param turn Planner that will drive the turns on this field.
  /// @return Mainland area
  F2CCells generateHeadlands(
    const F2CCells& field, const F2CRobot& robot, f2c::pp::TurningBase& turn);

  /// How far a turn between two neighbouring swaths reaches past their ends.
  ///
  /// This is what a border the swaths end on has to leave room for. It is not
  /// twice the turning radius: with room to spare the turn only reaches one
  /// radius out, and with the swaths closer together than that the turn has to
  /// loop and reaches further than two. A turn that backs up instead of
  /// driving round reaches no distance at all.
  /// @param robot Robot doing the coverage.
  /// @param turn Planner that will drive the turns on this field.
  /// @return Distance the turn reaches past the end of the swaths
  double turnExtent(const F2CRobot& robot, f2c::pp::TurningBase& turn) const;

  /// How far a turn reaches past swaths that meet the border at an angle.
  ///
  /// Swaths not square to a border cross it covWidth / sin(angle) apart, so
  /// the turn between them is a different one. pi/2 is the square case, which
  /// is what turnExtent(robot, turn) asks.
  /// @param robot Robot doing the coverage.
  /// @param turn Planner that will drive the turns on this field.
  /// @param track_border_angle Angle between the swath track and the border.
  /// @return Distance the turn reaches past the end of the swaths
  double turnExtent(const F2CRobot& robot, f2c::pp::TurningBase& turn,
    double track_border_angle) const;

  /// Open a corridor whose depth follows how each cell's swaths meet the
  /// border, instead of one depth for the whole field.
  ///
  /// turnExtent() answers the worst case: swaths ending square on the
  /// border. Each of the two cells is asked instead whether its swaths end
  /// on this border at all, and if they do how far a turn there reaches; the
  /// deeper answer is the depth.
  ///
  /// Whether they end there is a count, covWidth / sin(angle) apart along the
  /// border, and a count below one is not the same as none: over 31 fields it
  /// read "none" on 207 borders swaths did end on. This is the narrowest
  /// corridor worth trying, not one certainly wide enough.
  /// @param field Cells that share borders, usually from a decomposition.
  /// @param robot Robot doing the coverage.
  /// @param turn Planner that will drive the turns on this field.
  /// @param angs Swath track angle per cell, in \a field's order. Take them
  ///        off a mainland already carved at turnExtent()'s depth, not off
  ///        the bare cells -- that is where the swaths are generated.
  /// @return Mainland area
  F2CCells generateHeadlands(
    const F2CCells& field, const F2CRobot& robot, f2c::pp::TurningBase& turn,
    const std::vector<double>& angs);

  /// Ask the angles again on the cells an uncut border leaves joined.
  ///
  /// generateHeadlands(field, robot, turn, angs) reads \a angs cell by cell,
  /// but a border it opens no corridor on is a border the two cells are one
  /// piece across, and the swaths on that piece do not have to run the way
  /// either half ran alone: two tall cells side by side make one wide cell.
  /// Every border is measured again on the joined cells, which can leave a
  /// border uncut that cell by cell looked like a corridor.
  ///
  /// The joined pieces are then asked whether the swath generator can sweep
  /// them, one swath per track line. A piece it cuts into more swaths than
  /// there are lines has a throat in it, and the borders that made the throat
  /// are opened back up to turnExtent()'s depth.
  /// @param field Cells that share borders, usually from a decomposition.
  /// @param robot Robot doing the coverage.
  /// @param turn Planner that will drive the turns on this field.
  /// @param angs Swath track angle per cell, in \a field's order, as
  ///        generateHeadlands(field, robot, turn, angs) takes them. They
  ///        decide which borders are uncut, and so what is joined.
  /// @param obj Objective the swath generators below are asked against.
  /// @param sg_angle Asked which way the swaths run on a joined cell. Only
  ///        the angle is taken, so a coarse step is enough.
  /// @param sg_check Asked for the swaths themselves, to tell a piece with a
  ///        throat from one that sweeps. May be the same as \a sg_angle.
  /// @return Mainland area
  F2CCells generateHeadlands(
    const F2CCells& field, const F2CRobot& robot, f2c::pp::TurningBase& turn,
    const std::vector<double>& angs, f2c::obj::SGObjective& obj,
    f2c::sg::SwathGeneratorBase& sg_angle,
    f2c::sg::SwathGeneratorBase& sg_check);

  /// Open a corridor wide enough for \a n_swaths passes.
  /// @param field Borders of the field and the obstacles on it.
  /// @param swath_width Width of each headland swath.
  /// @param n_swaths Number of headland swaths.
  /// @return Mainland area
  F2CCells generateHeadlandArea(
    const F2CCells& field, double swath_width, int n_swaths) override;

  /// Widen the corridor one swath at a time.
  /// @param field Borders of the field and the obstacles on it.
  /// @param swath_width Width of each headland swath.
  /// @param n_swaths Number of headland swaths.
  /// @param dir_out2in When true, the widest corridor comes first.
  /// @return Vector of size \a n_swaths, each one carved a swath wider.
  std::vector<F2CCells> generateHeadlandSwaths(
    const F2CCells& field, double swath_width, int n_swaths,
    bool dir_out2in = true) override;

 private:
  CorridorShareMode share_mode_ {CorridorShareMode::ASYMMETRIC};

  /// Floor on sin(angle) so a track along the border does not put the two
  /// swath ends turnExtent() plans between infinitely far apart.
  double min_track_sin_ {1e-2};

  /// Part of turnExtent() a cell has to be wider than for its swath angle to
  /// say anything about which of its borders turns reach.
  double thin_cell_share_ {0.25};

  /// Tolerance the neighbour is buffered by to find the shared border.
  double tol_ {1e-3};
  /// Width the zero-width spur a difference can leave behind is opened by.
  double spur_ {1e-9};
  /// Tolerance two perimeters are compared with to count as the same size.
  double same_size_tol_ {1e-9};
  /// Shortest border kept as real: buffering the neighbour by tol_ turns a
  /// shared corner into a piece a few millimetres long on each edge that
  /// reaches it, and a border that short is a corner, not a corridor.
  double min_border_ {1e-2};
};

}  // namespace f2c::hg


#endif  // FIELDS2COVER_HEADLAND_GENERATOR_CORRIDOR_HEADLAND_H_
