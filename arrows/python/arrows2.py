import vtk
import numpy as np
from math import sqrt

def create_graphics_geometry(geom):
    """ Create geometry for display.
    """
    mapper = vtk.vtkPolyDataMapper()
    #mapper.SetInputConnection(geom)
    mapper.SetInputData(geom)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    #actor.GetProperty().SetLineWidth(3)
    #actor.GetProperty().EdgeVisibilityOn()
    return actor

def createData():

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName("bct.vtp")
    reader.Update()
    surface = reader.GetOutput()
    num_points = surface.GetPoints().GetNumberOfPoints()
    surf_points = surface.GetPoints()

    points = vtk.vtkPoints()
    points.SetNumberOfPoints(num_points)
    lines = vtk.vtkCellArray()
    lines.InsertNextCell(num_points+1)
    n = 0

    for i in range(num_points):
        pt = surf_points.GetPoint(i)
        points.SetPoint(n, pt[0], pt[1], pt[2])
        lines.InsertCellPoint(n)
        n += 1

    lines.InsertCellPoint(0)
    poly = vtk.vtkPolyData()
    poly.SetPoints(points)
    poly.SetLines(lines)

    arrow = vtk.vtkArrowSource()
    arrow.Update()

    return poly, arrow

# Create renderer and graphics window.
renderer = vtk.vtkRenderer()
renderer_win = vtk.vtkRenderWindow()
renderer_win.AddRenderer(renderer)
renderer.SetBackground(0.8, 0.8, 0.8)
renderer_win.SetSize(500, 500)

# Create input.
#    poly: vtk.vtkPolyData()
#    vel:  np.ndarray with shape (nPoints, 3)
#    mag:  np.ndarray with shape (nPoints,)
poly, arrow = createData()
num_points = poly.GetNumberOfPoints()

## Create velocity.
velocity = vtk.vtkDoubleArray()
velocity.SetName("velocity")
velocity.SetNumberOfComponents(3)
velocity.SetNumberOfTuples(num_points)
for i in range(num_points):
    vel = [0.0, 0.0, 10.0]
    velocity.SetTuple(i, vel)

## Create magnitude.
magnitude = vtk.vtkDoubleArray()
magnitude.SetNumberOfValues(num_points)
magnitude.SetName("magnitude")
for i in range(num_points):
    vel = velocity.GetTuple(i)
    mag = sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2])
    magnitude.SetValue(i, mag)

# Add to point data array.
poly.GetPointData().AddArray(velocity)
poly.GetPointData().AddArray(magnitude)
poly.GetPointData().SetActiveScalars("magnitude")
poly.GetPointData().SetActiveVectors("velocity")

# Create glyph.
glyph = vtk.vtkGlyph3D()
glyph.SetInputData(poly)
glyph.SetSourceConnection(arrow.GetOutputPort())
glyph.SetScaleFactor(0.1)
# glyph.OrientOn()
glyph.SetVectorModeToUseVector()
glyph.SetColorModeToColorByScalar()
glyph.Update()

actor = create_graphics_geometry(poly)
renderer.AddActor(actor)

actor = create_graphics_geometry(glyph.GetOutput())
renderer.AddActor(actor)

interactor = vtk.vtkRenderWindowInteractor()
interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
interactor.SetRenderWindow(renderer_win)
interactor.Start()







