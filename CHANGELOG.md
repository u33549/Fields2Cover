# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `SingleCellSwathsOrderBase::setStartAndEndPoint`, so BOUSTROPHEDON, SNAKE, SPIRAL and CUSTOM can start where the machine actually stands. The route begins at the cell nearest that point and returns to it, the same contract `RoutePlannerBase` already offered.
- `SingleCellSwathsOrderBase::genRoute`, which gives BOUSTROPHEDON, SNAKE, SPIRAL and CUSTOM a route whose connections are driven through the headland. `genSortedSwaths` is unchanged and still returns the bare order; a bare order ignores the boundary, so a snake or spiral skip cuts over covered ground on a field that is not convex.
- `SingleCellSwathsOrderBase::genSortedSwaths` overload taking `F2CSwathsByCells`, which orders each cell on its own so a pattern never runs across cells.
- `f2c::rp::RouteGeneratorBase`, the common interface of both route planner families. BOUSTROPHEDON, SNAKE, SPIRAL, CUSTOM and the TSP planner all answer `genRoute(cells, swaths_by_cells, d_tol)`, so a caller can hold any of them behind one type instead of branching on the mode. `RoutePlannerBase` keeps its existing `genRoute` overload, which is still the way to reach the optimizer settings.
- `f2c::pp::TurningBase::setFreeSpace`, which tells a turn planner where it may drive. A turn planner answers in free space: it is optimal there and keeps its radius, but it does not know where the crop is, so on a headland only as wide as the turn needs it leaves the drivable ground exactly where that ground turns inward. Given the ground, the shortest turn that stays on it is kept: another of the planner's own turns where it has one -- `alternativeTurns`, which Reeds-Shepp answers with its forward turn, since its shortest reverses over the crop for no gain -- and otherwise the turn split in two around a waypoint set inside a concave corner of the ground (`concaveCorners`). `createTurn` takes an optional `TurnReport` saying what the turn had to do and, when none fits, how deep the one it returns goes. `setSwathWidth` tells the planner how wide the swaths it joins are, so ground those swaths cover anyway is not counted against the turn, and `setWaypointOffset` how far inside a corner the waypoint sits. Left unset, every planner answers exactly as before.
- `f2c::rp::FreeSpaceRoutePlanner::setTurnRoom` and `getTurnRoom`, the room a connection's corners are given for the turn that follows them. The clearance is a price on the graph; this is the geometry of the turn, and the two are the same number only by habit. Zero leaves every connection where the graph put it; left unset, the clearance answers for it, so a caller that only sets a clearance sees no difference.

### Fixed
- `F2CGraph2D::getNodes` handed its nodes back in the order an `unordered_map` happened to hold them. That is not the order the node ids were handed out in, and it is not the same order on another build, yet the ids `getEdges` reports index exactly that vector -- so a caller resolving an edge through it joined the wrong pair of points, silently, and differently on another machine. The nodes now come back in id order, so `getNodes()[i]` is `indexToNode(i)`.
- `f2c::rp::FreeSpaceRoutePlanner` joined every swath end straight to whatever node it could see, so a connection left the swath sideways. A swath end stands on the border of the ground and the machine there is still on the swath's own line; a chord off that end starts across the line at an angle no machine can hold, and the turn planner then works the first metres of the connection over the crop. Where a turn room is set, each end now gets a node one step in along its own axis, joined to the end unconditionally -- the step is the machine's own line, not a shortcut the graph has to justify -- and the end leaves only that way. An end with no ground along its axis, or standing further off it than the strip its own swath covers -- buried in the crop where the ground has pinched shut, instead of resting on its border, so that the axis finds ground again only on the far side -- keeps its old edges instead, so no swath is cut off and no entry is reached by driving through the crop.
- `f2c::rp::FreeSpaceRoutePlanner` handed every connection to the path planner running through the corner of the ground it turns at. A shortest path is a geodesic and a geodesic hugs whatever it goes around, so the leg arrives at the border of the crop and leaves along it, and a turn of any radius has to cut into the crop to follow that. Each leg that runs along a border is now laid parallel to it and held `turn room (1 - cos(half the corner it turns))` clear of the border its turn cuts towards -- the offset at which an arc of that radius passes the corner instead of crossing it -- and a corner whose legs both reach a swath end, which cannot be moved, is pushed off that border itself. A leg already further out stays where it is, and a connection whose aligned legs would leave the ground is kept as it was.
- `SingleCellSwathsOrderBase::genRoute` drove the cells in the order the decomposition happened to produce them, and always entered each one from the same end. That order carries no geometric meaning, so the route crossed the field to reach a cell it could have taken on the way, and entered cells against the direction it arrived from. Cells are now ordered by what it costs to reach them -- distance measured through the headland graph, plus the turn needed to leave one cell and line up on the next -- and a cell is driven in reverse when its far end is the nearer one. A single cell has no order to choose and is untouched.
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














