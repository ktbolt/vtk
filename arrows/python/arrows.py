import vtk
import numpy as np

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
    source = vtk.vtkSphereSource()
    source.SetPhiResolution(30)
    source.SetThetaResolution(30)
    #source.SetCenter(10.0,1.0,30.0)
    source.Update()
    poly = source.GetOutput()
    arrow = vtk.vtkArrowSource()
    arrow.Update()
    # Compute some velocity field. This data you already got.
    # I'm creating here some data myself.
    nPoints = poly.GetNumberOfPoints()
    points = [poly.GetPoint(i) for i in range(nPoints)]
    points = np.asarray(points)
    x,y,z = points.T
    vel = np.c_[x*np.cos(y*np.pi), y*np.sin(2*x*np.pi), np.cos(2*z*np.pi)]
    mag = np.linalg.norm(vel, axis=1)
    return poly, vel, mag, arrow

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
poly, vel, mag, arrow = createData()
nPoints = poly.GetNumberOfPoints()

# Create data arrays.
velocity = vtk.vtkDoubleArray()
velocity.SetName("velocity")
velocity.SetNumberOfComponents(3)
velocity.SetNumberOfTuples(nPoints)
for i in range(nPoints):
    velocity.SetTuple(i, list(vel[i]))
magnitude = vtk.vtkDoubleArray()
# Similar as SetNumberOfValues(1) + SetNumberOfTuples(nPoints):
magnitude.SetNumberOfValues(nPoints)
magnitude.SetName("magnitude")
for i in range(nPoints):
    magnitude.SetValue(i, mag[i])

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







