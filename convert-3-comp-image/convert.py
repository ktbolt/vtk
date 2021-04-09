# Convert a 3-component VTI to single component.
import vtk

# Read the file.
file_name = 'image.vti'
reader = vtk.vtkXMLImageDataReader()
reader.SetFileName(file_name)
reader.Update()

# Get the image data.
image = reader.GetOutput()
dims = image.GetDimensions()
extent = image.GetExtent()
origin = image.GetOrigin()
spacing = image.GetSpacing()
scalars = image.GetPointData().GetScalars()

# Convert RGB to grayscale.
luminance_filter = vtk.vtkImageLuminance()
luminance_filter.SetInputData(image)
luminance_filter.Update()

# Write converted VTI.
file_name = 'converted.vti'
writer = vtk.vtkXMLImageDataWriter()
writer.SetFileName(file_name)
writer.SetInputData(luminance_filter.GetOutput())
writer.Write()

