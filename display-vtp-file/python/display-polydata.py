#!/usr/bin/env python

"""
This script is an example of how to create a cylinder aligned with an axix.
"""

import math
import vtk

# Read vtp file.
file_name = "polydata.vtp" 
reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName(file_name)
reader.Update()

polydata = reader.GetOutput()
print("Number of points: " + str(polydata.GetNumberOfPoints()))
print("Number of cells: " + str(polydata.GetNumberOfCells()))

# Visualize
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(reader.GetOutputPort())

actor = vtk.vtkActor()
actor.SetMapper(mapper)
#actor.GetProperty().SetColor(colors.GetColor3d("BurlyWood"))

renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)
renderWindowInteractor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

renderer.AddActor(actor)
#renderer.SetBackground(colors.GetColor3d('Teal'))
renderer.GetActiveCamera().Pitch(90)
renderer.GetActiveCamera().SetViewUp(0, 0, 1)
renderer.ResetCamera()

renderWindow.SetSize(600, 600)
renderWindow.Render()
renderWindow.SetWindowName('ReadPolyData')
renderWindowInteractor.Start()


