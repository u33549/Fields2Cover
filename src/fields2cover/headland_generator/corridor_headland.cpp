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
F2CCells dropNarrowerThanImplement(const F2CCells& carved, double cov_width) {
  F2CCells wide;
  for (size_t i = 0; i < carved.size(); ++i) {
    const F2CCell piece = carved.getGeometry(i);
    if (F2CCells::buffer(piece, -0.5 * cov_width).area() > 0.0) {
      wide.addGeometry(piece);
    }
  }
  return wide;
}
}  // namespace

double CorridorHL::turnExtent(
    const F2CRobot& robot, f2c::pp::TurningBase& turn) const {
  // Lay the ends of two neighbouring swaths on the x axis and turn from one
  // into the other. How far the turn reaches is the highest point it drives
  // to, which is the room a border the swaths end on has to leave.
  const F2CPath path = turn.createTurn(robot,
      F2CPoint(0.0, 0.0), M_PI_2,
      F2CPoint(robot.getCovWidth(), 0.0), -M_PI_2);
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
  return dropNarrowerThanImplement(carved, robot.getCovWidth());
}

F2CCells CorridorHL::generateHeadlands(
    const F2CCells& field, const F2CRobot& robot,
    f2c::pp::TurningBase& turn, const std::vector<double>& angs) {
  if (angs.size() != field.size()) {
    throw std::invalid_argument(
        "CorridorHL::generateHeadlands: angs needs one entry per cell");
  }
  const double half_w = 0.5 * robot.getWidth();
  const double raw_extent = turnExtent(robot, turn);
  // Where a border meets a cell's swaths head-on, the turn there needs
  // turnExtent()'s full reach; where the swaths run along it instead,
  // nothing turns there and only the implement's width matters. Between the
  // two, scale by how far the track is from parallel to the border, and take
  // whichever of the two cells on it asks for more room.
  const auto width_of = [&angs, half_w, raw_extent](
      const CorridorShare& share, size_t j) {
    const F2CLineString seg = share.shared_border.getGeometry(j);
    double border_ang = 0.0;
    if (seg.size() > 1) {
      border_ang =
          (seg.getGeometry(1) - seg.getGeometry(0)).getAngleFromPoint();
    }
    const double v = std::max(
        std::abs(std::sin(angs[share.cell_i] - border_ang)),
        std::abs(std::sin(angs[share.cell_k] - border_ang)));
    return half_w + v * raw_extent;
  };
  const F2CCells carved = carveCorridors(
      field, corridorShares(field, share_mode_), width_of, spur_);
  return dropNarrowerThanImplement(carved, robot.getCovWidth());
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
