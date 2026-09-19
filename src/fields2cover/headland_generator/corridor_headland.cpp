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
// depth actually cut on one border segment, the cell's share included.
F2CCells carveCorridors(const F2CCells& field,
    const std::vector<CorridorShare>& shares,
    const std::function<double(size_t, size_t)>& width_of,
    double spur) {
  F2CCells carved;
  for (size_t i = 0; i < field.size(); ++i) {
    F2CCells cell {field.getGeometry(i)};
    for (size_t s = 0; s < shares.size(); ++s) {
      if (shares[s].cell_i != i || shares[s].share <= 0.0) {
        continue;  // this cell gives no corridor on that border
      }
      for (size_t j = 0; j < shares[s].shared_border.size(); ++j) {
        const double width = width_of(s, j);
        if (width <= 0.0) {
          continue;  // no corridor here, and nothing to take out for it
        }
        cell = cell.difference(F2CCells::buffer(
            shares[s].shared_border.getGeometry(j), width));
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

// Where a cell sits, cheaply: a point inside it that does not move when the
// cell is carved, so the piece it came back as can be found again.
F2CPoint bboxCentre(const F2CCell& cell) {
  const F2CLinearRing ring = cell.getExteriorRing();
  double xmin = 1e18, xmax = -1e18, ymin = 1e18, ymax = -1e18;
  for (size_t i = 0; i < ring.size(); ++i) {
    const F2CPoint p = ring.getGeometry(i);
    xmin = std::min(xmin, p.getX()); xmax = std::max(xmax, p.getX());
    ymin = std::min(ymin, p.getY()); ymax = std::max(ymax, p.getY());
  }
  return F2CPoint(0.5 * (xmin + xmax), 0.5 * (ymin + ymax));
}

// turnExtent() plans a turn, and every border asks twice, so it is asked once
// per whole degree up front and read off between.
double reachAt(const std::vector<double>& reach, double ang) {
  const double deg = std::asin(std::min(1.0, std::abs(std::sin(ang))))
      * 180.0 / M_PI;
  const size_t lo = static_cast<size_t>(std::floor(deg));
  const size_t hi = std::min<size_t>(90, lo + 1);
  const double f = deg - lo;
  return reach[lo] * (1.0 - f) + reach[hi] * f;
}

// A cell too narrow to hold a turn spills onto its own borders whatever angle
// its swaths run at, so the count in borderWidth() cannot be trusted there.
std::vector<bool> thinCells(const F2CCells& field, double erode) {
  std::vector<bool> thin(field.size());
  for (size_t i = 0; i < field.size(); ++i) {
    thin[i] = F2CCells::buffer(field.getGeometry(i), -erode).area() <= 0.0;
  }
  return thin;
}

// Ask each cell separately: do your swaths end here, and if so how far does a
// turn here reach. Asking both at once lets a cell that never turns on this
// border set its depth. The cell's share of the answer is not taken here.
double borderWidth(const CorridorShare& share, size_t j,
    const std::vector<double>& angs, const std::vector<bool>& thin,
    const std::vector<double>& reach, double half_w, double cov_width) {
  const F2CLineString seg = share.shared_border.getGeometry(j);
  double border_ang = 0.0;
  if (seg.size() > 1) {
    border_ang = (seg.getGeometry(1) - seg.getGeometry(0)).getAngleFromPoint();
  }
  const double len = seg.length();
  const bool trust = !thin[share.cell_i] && !thin[share.cell_k];
  const auto demand = [&](size_t cell) {
    const double ang = angs[cell] - border_ang;
    // Swaths cross covWidth / sin(ang) apart: a shorter border gets none.
    if (trust && len * std::abs(std::sin(ang)) < cov_width) {
      return 0.0;
    }
    return reachAt(reach, ang);
  };
  const double req = std::max(demand(share.cell_i), demand(share.cell_k));
  // Neither cell turns here, so nobody drives here either: cutting the
  // implement's half width would only leave a strip nothing works.
  return req > 0.0 ? half_w + req : 0.0;
}

// The generator is asked for the swaths on a piece and counted twice: how
// many swaths came back, and how many distinct track lines they sit on. One
// line broken into two swaths is a throat -- something the corridor left too
// narrow to drive through, so the piece cannot be swept in one pass per line.
bool sweepsOnePassPerLine(const F2CCell& piece,
    f2c::sg::SwathGeneratorBase& sg, f2c::obj::SGObjective& obj,
    double cov_width, double line_tol) {
  const F2CSwaths swaths = sg.generateBestSwaths(obj, cov_width, piece);
  if (swaths.size() == 0) {
    return false;  // nothing could be swept on it at all
  }
  const double ang = swaths.at(0).getInAngle();
  const double sin_a = std::sin(ang), cos_a = std::cos(ang);
  std::vector<double> offsets;
  for (auto&& swath : swaths) {
    const F2CPoint a = swath.startPoint(), b = swath.endPoint();
    offsets.push_back(-0.5 * (a.getX() + b.getX()) * sin_a
                     + 0.5 * (a.getY() + b.getY()) * cos_a);
  }
  std::sort(offsets.begin(), offsets.end());
  size_t lines = 1;
  for (size_t i = 1; i < offsets.size(); ++i) {
    if (offsets[i] - offsets[i - 1] > line_tol * cov_width) {
      ++lines;
    }
  }
  return swaths.size() <= lines;
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
  const std::vector<CorridorShare> shares = corridorShares(field, share_mode_);
  const F2CCells carved = carveCorridors(field, shares,
      [width, &shares](size_t s, size_t) {
        return width * shares[s].share;
      }, spur_);
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
  std::vector<double> reach(91);
  for (int deg = 0; deg <= 90; ++deg) {
    reach[deg] = turnExtent(robot, turn, deg * M_PI / 180.0);
  }
  const std::vector<bool> thin =
      thinCells(field, thin_cell_share_ * turnExtent(robot, turn));
  const std::vector<CorridorShare> shares = corridorShares(field, share_mode_);
  // Where the corridor is nothing the two pieces meet along a line, which is
  // not a polygon a caller can use. Joining them also undoes a split the
  // decomposition made where both sides wanted the same swath direction.
  const F2CCells carved = carveCorridors(field, shares,
      [&](size_t s, size_t j) {
        return borderWidth(shares[s], j, angs, thin, reach, half_w,
            robot.getCovWidth()) * shares[s].share;
      }, spur_).unionCascaded();
  // A piece has to hold two passes to be worth keeping: one that only holds
  // the implement has nowhere to turn between them.
  return dropNarrowerThan(carved, robot.getCovWidth());
}

F2CCells CorridorHL::generateHeadlands(
    const F2CCells& field, const F2CRobot& robot, f2c::pp::TurningBase& turn,
    const std::vector<double>& angs, f2c::obj::SGObjective& obj,
    f2c::sg::SwathGeneratorBase& sg_angle,
    f2c::sg::SwathGeneratorBase& sg_check) {
  if (angs.size() != field.size()) {
    throw std::invalid_argument(
        "CorridorHL::generateHeadlands: angs needs one entry per cell");
  }
  const double half_w = 0.5 * robot.getWidth();
  const double cov_width = robot.getCovWidth();
  const double square_reach = turnExtent(robot, turn);
  const double full_width = square_reach + half_w;
  std::vector<double> reach(91);
  for (int deg = 0; deg <= 90; ++deg) {
    reach[deg] = turnExtent(robot, turn, deg * M_PI / 180.0);
  }
  const double erode = thin_cell_share_ * square_reach;

  // Read off the angles as given which borders no corridor is opened on. The
  // cells on either side of one of those are one piece, so join them.
  const std::vector<CorridorShare> shares = corridorShares(field, share_mode_);
  const std::vector<bool> thin = thinCells(field, erode);
  std::vector<size_t> root(field.size());
  for (size_t i = 0; i < root.size(); ++i) {
    root[i] = i;
  }
  const auto find = [&root](size_t i) {
    while (root[i] != i) {
      i = root[i] = root[root[i]];
    }
    return i;
  };
  for (size_t s = 0; s < shares.size(); ++s) {
    if (shares[s].share <= 0.0) {
      continue;
    }
    bool uncut = true;
    for (size_t j = 0; j < shares[s].shared_border.size(); ++j) {
      if (borderWidth(shares[s], j, angs, thin, reach, half_w, cov_width)
          * shares[s].share > 0.0) {
        uncut = false;
      }
    }
    if (uncut) {
      root[find(shares[s].cell_i)] = find(shares[s].cell_k);
    }
  }
  std::vector<F2CCells> group(field.size());
  for (size_t i = 0; i < field.size(); ++i) {
    group[find(i)].addGeometry(field.getGeometry(i));
  }
  F2CCells joined;
  for (size_t g = 0; g < group.size(); ++g) {
    if (group[g].size() == 0) {
      continue;
    }
    const F2CCells one = group[g].unionCascaded();
    for (size_t i = 0; i < one.size(); ++i) {
      joined.addGeometry(one.getGeometry(i));
    }
  }

  // Ask every angle again on the joined cells, off a mainland carved at the
  // depth turnExtent() asks for -- where the swaths are generated, which is
  // the same place the caller was told to take angs off.
  const F2CCells uniform = generateHeadlands(joined, robot, turn);
  std::vector<double> joined_angs(joined.size(), 0.0);
  for (size_t i = 0; i < joined.size(); ++i) {
    const F2CCell cell = joined.getGeometry(i);
    double best_area = 0.0;
    int best = -1;
    for (size_t p = 0; p < uniform.size(); ++p) {
      const F2CCell piece = uniform.getGeometry(p);
      if (!cell.isPointIn(bboxCentre(piece))) {
        continue;
      }
      if (piece.area() > best_area) {
        best_area = piece.area();
        best = static_cast<int>(p);
      }
    }
    joined_angs[i] = sg_angle.computeBestAngle(obj, cov_width,
        best >= 0 ? uniform.getGeometry(best) : cell);
  }
  const std::vector<CorridorShare> joined_shares =
      corridorShares(joined, share_mode_);
  const std::vector<bool> joined_thin = thinCells(joined, erode);
  std::vector<std::vector<double>> widths(joined_shares.size());
  for (size_t s = 0; s < joined_shares.size(); ++s) {
    widths[s].assign(joined_shares[s].shared_border.size(), 0.0);
    if (joined_shares[s].share <= 0.0) {
      continue;
    }
    for (size_t j = 0; j < widths[s].size(); ++j) {
      widths[s][j] = borderWidth(joined_shares[s], j, joined_angs, joined_thin,
          reach, half_w, cov_width) * joined_shares[s].share;
    }
  }
  const auto carve = [&]() {
    return dropNarrowerThan(carveCorridors(joined, joined_shares,
        [&widths](size_t s, size_t j) { return widths[s][j]; }, spur_),
        cov_width);
  };
  F2CCells mainland = carve();

  // The depth above is the narrowest worth trying, not one certainly wide
  // enough. Ask the generator to sweep what came back: a piece it cuts into
  // more swaths than there are track lines has a throat, and the borders that
  // reach the throat are opened back up to the depth taken without an angle.
  for (int pass = 0; pass < cert_passes_; ++pass) {
    std::vector<F2CCell> throats;
    for (size_t p = 0; p < mainland.size(); ++p) {
      const F2CCell piece = mainland.getGeometry(p);
      if (!sweepsOnePassPerLine(piece, sg_check, obj, cov_width, line_tol_)) {
        throats.push_back(piece);
      }
    }
    if (throats.empty()) {
      break;
    }
    bool widened = false;
    for (size_t s = 0; s < joined_shares.size(); ++s) {
      if (joined_shares[s].share <= 0.0) {
        continue;
      }
      for (size_t j = 0; j < widths[s].size(); ++j) {
        if (widths[s][j] >= full_width - 1e-9) {
          continue;  // already as deep as it can be opened
        }
        const F2CLineString seg = joined_shares[s].shared_border.getGeometry(j);
        if (seg.size() < 2) {
          continue;
        }
        const F2CPoint a = seg.getGeometry(0);
        const F2CPoint b = seg.getGeometry(seg.size() - 1);
        const F2CPoint mid(0.5 * (a.getX() + b.getX()),
                           0.5 * (a.getY() + b.getY()));
        for (const F2CCell& throat : throats) {
          if (throat.isPointIn(mid) || mid.distance(throat) < 0.5 * cov_width) {
            widths[s][j] = full_width * joined_shares[s].share;
            widened = true;
            break;
          }
        }
      }
    }
    if (!widened) {
      break;
    }
    mainland = carve();
  }
  return mainland;
}

F2CCells CorridorHL::generateHeadlands(
    const F2CCells& field, double dist_headland) {
  const std::vector<CorridorShare> shares = corridorShares(field, share_mode_);
  return carveCorridors(field, shares,
      [dist_headland, &shares](size_t s, size_t) {
        return dist_headland * shares[s].share;
      }, spur_);
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
