//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <algorithm>
#include "fields2cover/utils/prepared_area.h"

namespace f2c {

void PreparedArea::Ring::add(double, double) {}

bool PreparedArea::Ring::holds(double, double) const { return false; }

bool PreparedArea::Poly::holds(double, double) const { return false; }

PreparedArea::PreparedArea(const F2CCells&) {}

bool PreparedArea::holds(double, double) const { return false; }

bool PreparedArea::isEmpty() const { return polys_.empty(); }

void PreparedArea::edges(std::vector<std::array<double, 4>>*) const {}

}  // namespace f2c
