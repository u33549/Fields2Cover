#==============================================================================
#     Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
#                      Author: Gonzalo Mier
#                         BSD-3 License
#==============================================================================

import pytest
import fields2cover as f2c

def near(a, b, error = 1e-7):
  assert a == pytest.approx(b, error)

def test_fields2cover_hl_corridor_gen_generateHeadlands():
  cells = f2c.Cells(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(0,0), f2c.Point(50,0), f2c.Point(50,100), f2c.Point(0,100), f2c.Point(0,0)]))));
  cells.addGeometry(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(50,0), f2c.Point(100,0), f2c.Point(100,100), f2c.Point(50,100), f2c.Point(50,0)]))));
  carved = f2c.HG_Corridor_gen().generateHeadlands(cells, 2.0);
  assert (carved.size() == 2);
  near(cells.area() - carved.area(), 200, 1e-3);

def test_fields2cover_hl_corridor_gen_corridorShares():
  cells = f2c.Cells(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(0,0), f2c.Point(20,0), f2c.Point(20,10), f2c.Point(0,10), f2c.Point(0,0)]))));
  cells.addGeometry(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(10,10), f2c.Point(20,10), f2c.Point(20,20), f2c.Point(10,20), f2c.Point(10,10)]))));
  shares = f2c.HG_Corridor_gen().corridorShares(cells);
  assert (len(shares) == 2);
  for share in shares:
    assert (not share.same_size);
    near(share.shared_length, 10, 1e-2);
  # The smaller cell gives the whole corridor, the larger one gives nothing.
  by_cell = {share.cell_i: share for share in shares};
  assert (by_cell[0].cell_k == 1 and by_cell[1].cell_k == 0);
  near(by_cell[0].share, 0.0, 1e-9);
  near(by_cell[1].share, 1.0, 1e-9);

def test_fields2cover_hl_corridor_gen_shareMode():
  corridor = f2c.HG_Corridor_gen();
  assert (corridor.getShareMode() == f2c.CorridorShareMode_ASYMMETRIC);

  cells = f2c.Cells(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(0,0), f2c.Point(20,0), f2c.Point(20,10), f2c.Point(0,10), f2c.Point(0,0)]))));
  cells.addGeometry(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(10,10), f2c.Point(20,10), f2c.Point(20,20), f2c.Point(10,20), f2c.Point(10,10)]))));

  # Sizes differ, so ASYMMETRIC would give the whole corridor to one side.
  shares = corridor.corridorShares(cells, f2c.CorridorShareMode_SYMMETRIC);
  assert (len(shares) == 2);
  for share in shares:
    near(share.share, 0.5, 1e-9);

  corridor.setShareMode(f2c.CorridorShareMode_SYMMETRIC);
  assert (corridor.getShareMode() == f2c.CorridorShareMode_SYMMETRIC);

def test_fields2cover_hl_corridor_gen_turnExtentFollowsThePlannedTurn():
  corridor = f2c.HG_Corridor_gen();
  dubins = f2c.PP_DubinsCurves();

  # Swaths far enough apart: the turn reaches one radius out, not two.
  roomy = f2c.Robot(2.0, 10.0);
  roomy.setMinTurningRadius(2.0);
  near(corridor.turnExtent(roomy, dubins), roomy.getMinTurningRadius(), 1e-2);

  # Swaths closer than twice the radius: the turn loops and reaches further
  # than the classic bound, which is the case that bound gets wrong.
  tight = f2c.Robot(2.0, 6.0);
  tight.setMinTurningRadius(5.0);
  assert (corridor.turnExtent(tight, dubins) > 2.0 * tight.getMinTurningRadius());

def test_fields2cover_hl_corridor_gen_corridorIsAsDeepAsThePlannedTurn():
  cells = f2c.Cells(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(0,0), f2c.Point(100,0), f2c.Point(100,50), f2c.Point(0,50), f2c.Point(0,0)]))));
  cells.addGeometry(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(0,50), f2c.Point(100,50), f2c.Point(100,100), f2c.Point(0,100), f2c.Point(0,50)]))));
  robot = f2c.Robot(2.0, 10.0);
  robot.setMinTurningRadius(2.0);
  corridor = f2c.HG_Corridor_gen();
  dubins = f2c.PP_DubinsCurves();

  width = corridor.turnExtent(robot, dubins) + 0.5 * robot.getWidth();
  carved = corridor.generateHeadlands(cells, robot, dubins);
  assert (carved.size() == 2);
  near(cells.area() - carved.area(), 100 * width, 1e-2);

def test_fields2cover_hl_corridor_gen_dropsAPieceNarrowerThanTheImplement():
  def band(x0, y0, y1):
    return f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
      [f2c.Point(x0,y0), f2c.Point(x0+100,y0), f2c.Point(x0+100,y1),
       f2c.Point(x0,y1), f2c.Point(x0,y0)])));
  cells = f2c.Cells(band(0, 0, 50));
  cells.addGeometry(band(0, 50, 55));     # 5 m -> 2 m left, too narrow
  cells.addGeometry(band(0, 55, 105));
  cells.addGeometry(band(200, 0, 50));
  cells.addGeometry(band(200, 50, 59));   # 9 m -> 6 m left, wide enough
  cells.addGeometry(band(200, 59, 109));
  robot = f2c.Robot(1.0, 4.0);
  robot.setMinTurningRadius(1.0);
  corridor = f2c.HG_Corridor_gen();

  carved = corridor.generateHeadlands(cells, robot, f2c.PP_DubinsCurves());
  assert (carved.size() == 5);
  near(carved.area(), 4 * 100 * 50 + 100 * 6, 1e-2);

def test_fields2cover_hl_corridor_gen_angsNarrowTheCorridorWhenSwathsRunAlongTheBorder():
  cells = f2c.Cells(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(0,0), f2c.Point(100,0), f2c.Point(100,50), f2c.Point(0,50), f2c.Point(0,0)]))));
  cells.addGeometry(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(0,50), f2c.Point(100,50), f2c.Point(100,100), f2c.Point(0,100), f2c.Point(0,50)]))));
  robot = f2c.Robot(2.0, 10.0);
  robot.setMinTurningRadius(2.0);
  corridor = f2c.HG_Corridor_gen();
  dubins = f2c.PP_DubinsCurves();

  # Both cells' swaths run along the shared (horizontal) border: no turn
  # happens there, so the corridor only has to fit the implement.
  angs = f2c.VectorDouble([0.0, 0.0]);
  carved = corridor.generateHeadlands(cells, robot, dubins, angs);
  expected_width = 0.5 * robot.getWidth();
  near(cells.area() - carved.area(), 100 * expected_width, 1e-2);
  assert (expected_width < corridor.turnExtent(robot, dubins));

def test_fields2cover_hl_corridor_gen_angsMatchTheUniformCorridorWhenSwathsMeetTheBorderHeadOn():
  cells = f2c.Cells(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(0,0), f2c.Point(100,0), f2c.Point(100,50), f2c.Point(0,50), f2c.Point(0,0)]))));
  cells.addGeometry(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(0,50), f2c.Point(100,50), f2c.Point(100,100), f2c.Point(0,100), f2c.Point(0,50)]))));
  robot = f2c.Robot(2.0, 10.0);
  robot.setMinTurningRadius(2.0);
  corridor = f2c.HG_Corridor_gen();
  dubins = f2c.PP_DubinsCurves();

  # Both cells' swaths meet the border head-on: turnExtent()'s worst case
  # already covers it, so the angle-aware corridor matches the uniform one.
  angs = f2c.VectorDouble([1.5707963267948966, 1.5707963267948966]);
  angled = corridor.generateHeadlands(cells, robot, dubins, angs);
  uniform = corridor.generateHeadlands(cells, robot, dubins);
  near(angled.area(), uniform.area(), 1e-2);

def test_fields2cover_hl_corridor_gen_angsWrongSizeThrows():
  cells = f2c.Cells(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(0,0), f2c.Point(100,0), f2c.Point(100,50), f2c.Point(0,50), f2c.Point(0,0)]))));
  cells.addGeometry(f2c.Cell(f2c.LinearRing(f2c.VectorPoint(
    [f2c.Point(0,50), f2c.Point(100,50), f2c.Point(100,100), f2c.Point(0,100), f2c.Point(0,50)]))));
  robot = f2c.Robot(2.0, 10.0);
  robot.setMinTurningRadius(2.0);
  corridor = f2c.HG_Corridor_gen();

  with pytest.raises(Exception) as e_info:
    corridor.generateHeadlands(cells, robot, f2c.PP_DubinsCurves(),
        f2c.VectorDouble([0.0]));
