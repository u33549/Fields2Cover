//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>
#include "fields2cover/headland_generator/corridor_headland.h"

namespace f2c::hg {

std::vector<CorridorShare> CorridorHL::corridorShares(
    const F2CCells& field, CorridorShareMode mode) const {
  // The corridor comes out of the smaller cell, so the larger neighbour keeps
  // its shape: taking half from each side notches both, and a notched cell
  // costs the swath generator a pass. In SYMMETRIC mode every border splits
  // evenly instead, regardless of size, for comparison against that rule.
  std::vector<CorridorShare> shares;
  for (size_t i = 0; i < field.size(); ++i) {
    const F2CLinearRing ring = field.getCellBorder(i);
    const double perimeter = ring.length();
    for (size_t k = 0; k < field.size(); ++k) {
      if (k == i) {
        continue;
      }
      CorridorShare share;
      share.cell_i = i;
      share.cell_k = k;
      share.perimeter_i = perimeter;
      share.perimeter_k = field.getCellBorder(k).length();
      // Cells that are the same size share the corridor, each giving half.
      // The comparison has to be loose: perimeters that are equal in theory
      // differ in the last bits, and letting that noise pick a winner leaves
      // some borders with a full corridor and others with none.
      share.same_size = std::abs(perimeter - share.perimeter_k) <=
          same_size_tol_ * perimeter;
      share.share = mode == CorridorShareMode::SYMMETRIC ? 0.5 :
          (share.same_size ? 0.5 :
              (perimeter < share.perimeter_k ? 1.0 : 0.0));
      // Keep the part of each border edge that the neighbour actually touches.
      const F2CCells neighbour = F2CCells::buffer(field.getGeometry(k), tol_);
      for (size_t e = 0; e + 1 < ring.size(); ++e) {
        F2CMultiLineString edge;
        edge.addGeometry(
            F2CLineString({ring.getGeometry(e), ring.getGeometry(e + 1)}));
        const F2CMultiLineString shared = edge.intersection(neighbour);
        for (size_t j = 0; j < shared.size(); ++j) {
          const F2CLineString part = shared.getGeometry(j);
          if (part.size() > 1 && part.length() > min_border_) {
            share.shared_border.addGeometry(part);
            share.shared_length += part.length();
          }
        }
      }
      // Cells that only meet at a corner, or not at all, are not neighbours.
      if (share.shared_length > 0.0) {
        shares.emplace_back(share);
      }
    }
  }
  return shares;
}

void CorridorHL::setShareMode(CorridorShareMode mode) {
  share_mode_ = mode;
}

CorridorShareMode CorridorHL::getShareMode() const {
  return share_mode_;
}

namespace {
// Take each cell's part of every corridor out of it. width_of gives the
// full width of one border segment, before the cell's share of it is taken.
F2CCells carveCorridors(const F2CCells& field,
    const std::vector<CorridorShare>& shares,
    const std::function<double(const CorridorShare&, size_t)>& width_of,
    double spur) {
  F2CCells carved;
  for (size_t i = 0; i < field.size(); ++i) {
    F2CCells cell {field.getGeometry(i)};
    for (const CorridorShare& share : shares) {
      if (share.cell_i != i || share.share <= 0.0) {
        continue;  // this cell gives no corridor on that border
      }
      for (size_t j = 0; j < share.shared_border.size(); ++j) {
        cell = cell.difference(F2CCells::buffer(
            share.shared_border.getGeometry(j),
            width_of(share, j) * share.share));
      }
    }
    // The difference can leave a zero-width spur behind. Opening the result by
    // a hair drops it without moving any real edge.
    const F2CCells clean = cell.buffer(-spur).buffer(spur);
    for (size_t j = 0; j < clean.size(); ++j) {
      carved.addGeometry(clean.getGeometry(j));
    }
  }
  return carved;
}

// Nothing can cover a piece narrower than the implement without it hanging
// over the corridor, and a swath generator handed one cuts it into fragments
// the route then joins with a turn each. Leave it to the corridor instead.
F2CCells dropNarrowerThan(const F2CCells& carved, double half_width) {
  F2CCells wide;
  for (size_t i = 0; i < carved.size(); ++i) {
    const F2CCell piece = carved.getGeometry(i);
    if (F2CCells::buffer(piece, -half_width).area() > 0.0) {
      wide.addGeometry(piece);
    }
  }
  return wide;
}
}  // namespace

double CorridorHL::turnExtent(
    const F2CRobot& robot, f2c::pp::TurningBase& turn) const {
  return turnExtent(robot, turn, M_PI_2);
}

double CorridorHL::turnExtent(const F2CRobot& robot,
    f2c::pp::TurningBase& turn, double track_border_angle) const {
  // Lay the border on the x axis with two swath ends on it, as far apart as
  // swaths at this angle cross it, and turn from one into the other.
  const double sin_ang = std::max(
      std::abs(std::sin(track_border_angle)), min_track_sin_);
  const double ang = std::asin(std::min(1.0, sin_ang));
  const F2CPath path = turn.createTurn(robot,
      F2CPoint(0.0, 0.0), ang,
      F2CPoint(robot.getCovWidth() / sin_ang, 0.0), ang + M_PI);
  // A planner that cannot join the two swaths says nothing about the width,
  // so fall back on the bound taken without one. A planner that joins them
  // without reaching past their ends is a different matter: a turn that backs
  // up instead of driving round needs no room there, and saying so is the
  // point of asking.
  if (path.size() == 0) {
    return 2.0 * robot.getMinTurningRadius();
  }
  double extent = 0.0;
  for (const auto& state : path.getStates()) {
    extent = std::max(extent, state.point.getY());
  }
  return extent;
}

F2CCells CorridorHL::generateHeadlands(
    const F2CCells& field, const F2CRobot& robot,
    f2c::pp::TurningBase& turn) {
  const double width = turnExtent(robot, turn) + 0.5 * robot.getWidth();
  const F2CCells carved = carveCorridors(field,
      corridorShares(field, share_mode_),
      [width](const CorridorShare&, size_t) { return width; }, spur_);
  return dropNarrowerThan(carved, 0.5 * robot.getCovWidth());
}

F2CCells CorridorHL::generateHeadlands(
    const F2CCells& field, const F2CRobot& robot,
    f2c::pp::TurningBase& turn, const std::vector<double>& angs) {
  if (angs.size() != field.size()) {
    throw std::invalid_argument(
        "CorridorHL::generateHeadlands: angs needs one entry per cell");
  }
  const double half_w = 0.5 * robot.getWidth();
  const double square_reach = turnExtent(robot, turn);
  // turnExtent() plans a turn, and every border asks twice.
  std::vector<double> reach(91);
  for (int deg = 0; deg <= 90; ++deg) {
    reach[deg] = turnExtent(robot, turn, deg * M_PI / 180.0);
  }
  // A cell too narrow to hold a turn spills onto its own borders whatever
  // angle its swaths run at, so the count below cannot be trusted there.
  std::vector<bool> thin(field.size());
  for (size_t i = 0; i < field.size(); ++i) {
    thin[i] = F2CCells::buffer(field.getGeometry(i),
        -thin_cell_share_ * square_reach).area() <= 0.0;
  }
  const auto reach_at = [&reach](double ang) {
    const double deg = std::asin(std::min(1.0, std::abs(std::sin(ang))))
        * 180.0 / M_PI;
    const size_t lo = static_cast<size_t>(std::floor(deg));
    const size_t hi = std::min<size_t>(90, lo + 1);
    const double f = deg - lo;
    return reach[lo] * (1.0 - f) + reach[hi] * f;
  };
  // Ask each cell separately: do your swaths end here, and if so how far does
  // a turn here reach. Asking both at once lets a cell that never turns on
  // this border set its depth.
  const auto width_of = [&angs, &reach_at, &thin, half_w, &robot](
      const CorridorShare& share, size_t j) {
    const F2CLineString seg = share.shared_border.getGeometry(j);
    double border_ang = 0.0;
    if (seg.size() > 1) {
      border_ang =
          (seg.getGeometry(1) - seg.getGeometry(0)).getAngleFromPoint();
    }
    const double len = seg.length();
    const bool trust = !thin[share.cell_i] && !thin[share.cell_k];
    const auto demand = [&](size_t cell) {
      const double ang = angs[cell] - border_ang;
      // Swaths cross covWidth / sin(ang) apart: a shorter border gets none.
      if (trust && len * std::abs(std::sin(ang)) < robot.getCovWidth()) {
        return 0.0;
      }
      return reach_at(ang);
    };
    const double reach =
        std::max(demand(share.cell_i), demand(share.cell_k));
    // Neither cell turns here, so nobody drives here either: cutting the
    // implement's half width would only leave a strip nothing works.
    return reach > 0.0 ? half_w + reach : 0.0;
  };
  // Where the corridor is nothing the two pieces meet along a line, which is
  // not a polygon a caller can use. Joining them also undoes a split the
  // decomposition made where both sides wanted the same swath direction.
  const F2CCells carved = carveCorridors(
      field, corridorShares(field, share_mode_), width_of, spur_)
      .unionCascaded();
  // A piece has to hold two passes to be worth keeping: one that only holds
  // the implement has nowhere to turn between them.
  return dropNarrowerThan(carved, robot.getCovWidth());
}

F2CCells CorridorHL::generateHeadlands(
    const F2CCells& field, const F2CRobot& robot, f2c::pp::TurningBase& turn,
    const std::vector<double>& angs, f2c::obj::SGObjective& obj,
    f2c::sg::SwathGeneratorBase& sg_angle,
    f2c::sg::SwathGeneratorBase& sg_check) {
  return generateHeadlands(field, robot, turn, angs);
}

F2CCells CorridorHL::generateHeadlands(
    const F2CCells& field, double dist_headland) {
  return carveCorridors(field, corridorShares(field, share_mode_),
      [dist_headland](const CorridorShare&, size_t) { return dist_headland; },
      spur_);
}

F2CCells CorridorHL::generateHeadlandArea(
    const F2CCells& field, double swath_width, int n_swaths) {
  return generateHeadlands(field, swath_width * n_swaths);
}

std::vector<F2CCells> CorridorHL::generateHeadlandSwaths(
    const F2CCells& field, double swath_width, int n_swaths, bool dir_out2in) {
  std::vector<F2CCells> hl;
  for (int i = 0; i < n_swaths; ++i) {
    const int n = dir_out2in ? (n_swaths - i) : (i + 1);
    hl.emplace_back(generateHeadlands(field, swath_width * (n - 0.5)));
  }
  return hl;
}


}  // namespace f2c::hg
