# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `SingleCellSwathsOrderBase::genRoute`, which gives BOUSTROPHEDON, SNAKE, SPIRAL and CUSTOM a route whose connections are driven through the headland. `genSortedSwaths` is unchanged and still returns the bare order; a bare order ignores the boundary, so a snake or spiral skip cuts over covered ground on a field that is not convex.
- `SingleCellSwathsOrderBase::genSortedSwaths` overload taking `F2CSwathsByCells`, which orders each cell on its own so a pattern never runs across cells.
- `f2c::hg::CorridorHL::generateHeadlands` overload taking each cell's swath track angle alongside the turn planner, so the corridor's depth follows how a border meets the swaths instead of one depth for the whole field. Each of the two cells on a border is asked whether its swaths end there and, if they do, how far a turn there reaches; the deeper answer is the depth. A border neither cell turns on is not cut at all and the two cells come back joined, which undoes a split the decomposition made where both sides wanted the same swath direction. A cell too narrow to hold a turn spills onto its borders whatever angle its swaths run at, so the question is not asked there and the border gets the full depth. Over 31 real fields it recovers 5.1% (small robot) to 11.2% (large) more mainland than the uniform depth and never less on any field, but it also drives further into the mainland on 16 of 31 fields and leaves more mainland pieces, so it is not a drop-in replacement for the uniform depth.
- `f2c::hg::CorridorHL::turnExtent` overload taking the angle the swaths meet the border at. Swaths not square to a border cross it `covWidth / sin(angle)` apart, so the turn between them is a different one; the reach that comes back barely moves with the angle and never falls below the turning radius, unlike scaling `turnExtent()` by `sin(angle)`.

### Fixed
- `F2CSwathsByCells::flatten` gave every cell's swaths the ids they had inside that cell, so the ids repeated across cells. `F2CSwaths::sort` orders on the id alone, so ordering a flattened multi-cell set interleaved the cells instead of covering them one after another. Flattened swaths are now numbered from 0.
- `f2c::rp::CustomOrder` accepted an order with duplicate or out-of-range values and only reported it when the order was used. It is now checked when it is set, and the order must be a permutation of `0..n-1`.
- `Point::intersectionOfLines` no longer sends the result arbitrarily far away for two lines a fraction of a degree apart. The parallel check compared the determinant to exactly 0, but two real-world collinear borders -- a redundant vertex on an otherwise straight, digitized field edge -- produce a determinant that is 0 only up to rounding, dividing by which is what actually blew up.
- `ReqHL::generateHeadlands` no longer collapses a mainland to nearly nothing on a heavily digitized border. Offsetting hundreds of segments by widely different amounts crosses the resulting ring many times over; cleaning that up fell to `LinearRing::filterSelfIntersections`, which resolves one crossing at a time and drops every point between the two segments involved -- most of a real border, if that crossing is a distant one. The exterior ring's offset is now cleaned up with a GEOS buffer instead, which can also come back as more than one mainland cell where the ring legitimately splits.
- `LinearRing::getParallelLine` no longer sends a corner's offset vertex arbitrarily far from the polygon when its two sides need very different offset distances. The plain line-line miter join has no limit on how far the intersection can land as the corner sharpens; past 4 times the larger of the two offsets it now bevels (the two offset segment endpoints, joined directly) instead of following the miter out. On real fields this let `ReqHL`'s mainland reach up to 273 m outside the field it came from, occasionally taking a coverage swath with it.

## [2.1.0] - 2026-09-03

### Added
- `f2c::hg::ReqHL`, a headland generator that sizes each border on its own: a border the swaths run along is only entered, while a border they end on takes a whole turn. The difference is left to the mainland instead of being given up on every border.
- `HeadlandGeneratorBase::generateHeadlands` overloads taking a robot and the track angles, and `maxHLWidthRequired`.
- `F2CLinearRing::getParallelLine`, `bufferOutwards`, `bufferInwards`, `filterSelfIntersections`, `removePoint`, `getSegment`, `getLastSegment`, `segmentLength` and `segmentAng`, to offset each segment of a ring by its own distance.
- `Geometry::contains`, the other side of `Geometry::within`.
- `f2c::hg::CorridorHL`, a headland generator that opens a corridor where cells border each other instead of shrinking every border. Only the part of an edge a neighbour actually touches is cut, and the corridor comes out of the smaller cell so the larger neighbour keeps its shape; cells of the same size split it evenly. Edges facing the outer boundary or a void are left untouched.
- `f2c::hg::CorridorHL::corridorShares`, the rule that generator applies, on its own: for each pair of cells that share a border it reports the two perimeters, whether they count as the same size, how much of the corridor each cell gives, and the border they share.
- `f2c::hg::CorridorShareMode`, to split every corridor evenly (`SYMMETRIC`) instead of giving it all to the smaller cell (`ASYMMETRIC`, the default). `CorridorHL::corridorShares` takes it as an argument, and `CorridorHL::setShareMode`/`getShareMode` set the mode `generateHeadlands` applies.
- `f2c::hg::CorridorHL::turnExtent`, how far a turn between two neighbouring swaths reaches past their ends, and a `generateHeadlands` overload taking the turn planner that will drive the field, which opens a corridor that deep. Twice the turning radius is only a bound: with the swaths far enough apart the turn reaches one radius out, and with them closer together than twice the radius it has to loop and reaches **further** than the bound allows for, so a corridor sized on the bound is too shallow there. A turn that backs up instead of driving round reaches no distance at all. A piece left narrower than the robot's coverage width is given to the corridor instead of returned: nothing can cover it without the implement hanging over the corridor, and a swath generator handed one cuts it into fragments the route then joins with a turn each.
- Python module is built as a proper package with scikit-build-core (`pip install .`); version is taken from `CMakeLists.txt` and exposed as `fields2cover.__version__`.
- Source distribution published to PyPI (`pip install fields2cover`).

### Fixed
- `f2c::hg::CorridorHL` no longer opens a corridor where two cells only meet at a corner. The neighbour is buffered by a tolerance to find the shared border, which turned a single shared point into a few millimetres of "border" on every edge reaching it; in a field of cells meeting at one point that carved a disc out of the middle and made every slice a neighbour of every other.
- `Cells::splitByLine` no longer throws `std::invalid_argument` when a split leaves a piece touching itself at a single point. Reinflating each split piece went through `Cell::buffer`, which only accepts a single polygon back; a positive buffer on a pinched piece can separate it into two. A `MultiLineString` split also no longer silently keeps only the first line's cut: reinflating after every individual line let GEOS collapse the next cut into a no-op, so every line is now buffered and cut in one pass instead of one after another.
- `F2CCells::getCellBorder`, `getInteriorRing` and `addRing` no longer segfault on an empty polygon or an out-of-range index; they throw `std::out_of_range` like `getGeometry` does.
- `NSwathModified::computeCost` read the wrong neighbouring point for the first edge: `(i - 1) % ring.size()` wraps a `size_t` to `SIZE_MAX % n`, which is not the previous vertex. The cost of a polygon now no longer depends on which vertex its ring starts from.
- `generateBestSwaths` no longer returns an angle that covers nothing. The objectives estimate the cost from the cell border alone, so a cell with a hairline spur could score best on an angle producing no swath at all and was silently left uncovered.

### Changed
- `CorridorHL`'s tolerances (the neighbour-buffer, spur, same-size and minimum-border thresholds) are private member variables instead of constants hidden in the .cpp file, visible directly on the class in the header.
- The decomposition tutorial carves a corridor between cells instead of running the headland generator a second time, which also shrank the outer boundary.
- `cmake --install` places the python module in the interpreter's site-packages instead of calling `setup.py install`.
- Building the python module requires CMake >= 3.18 and Python >= 3.9.

## [2.0.0] - 07-02-2024

- Route planner travelling through the headlands

## [1.3.0] - 21-04-2023

- Add decomposition algorithms: trapezoidal, boustrophedon
  


## [1.2.0] - 17-10-2022
### Added
- Tests to do cover < 90% functions

### Changes
- SG use the objective function as a parameter instead of a template.
- RP do not save the swaths and modify them using the functions provided
- PP do not save the robot and use the robot params with a param on the function.

### Changes
- Objectives are split for each of the modules.
- Global objective renamed to SG objective.
- Path objective renamed to RP objective.

### Added
- PP objective
- HL objective




## [1.1.0]
### Added
- On HL module: constant headland algorithm.
- On SG module: brute force algorithm.
- On RP module: Boustrophedon, custom, snake and spiral.
- On PP module: Dubins and Reeds-Sheep with/without continuous curvature.
- Objective functions are split between global and path cost functions.














