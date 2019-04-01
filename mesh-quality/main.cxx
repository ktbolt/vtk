
// This is an example of using the vtkMeshQuality class to check the quaility of an unstructured mesh. 
//
// The 'SetTetQualityMeasureToVolume' estimator is used to check element volumes. Elements that have a 
// volume <= 0 are counted. A graphics window is opened and the mesh is disaplayed with bad elements outlined in red.
//
// There are many estimators: SetTetQualityMeasureToEdgeRatio, SetTetQualityMeasureToAspectRatio, 
// SetTetQualityMeasureToRadiusRatio, etc. See https://vtk.org/doc/nightly/html/classvtkMeshQuality.html for
// a complete list.
//
// Usage:
//
//     mesh-quality <FileName>.vtu
//

#include <string>
#include <bits/stdc++.h> 

#include <vtkActor.h>
#include <vtkCellData.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellPicker.h>
#include <vtkCellType.h>
#include <vtkGenericCell.h>
#include <vtkGeometryFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkExtractSelection.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkMeshQuality.h>
#include <vtkNamedColors.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkShrinkFilter.h>
#include <vtkSmartPointer.h>
#include <vtkThreshold.h>
#include <vtkTriangleFilter.h>
#include "vtkUnsignedCharArray.h"
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

//---------------------- 
// MouseInteractorStyle
//---------------------- 
// This class is used to select cells.
//
class MouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
  static MouseInteractorStyle* New();

  MouseInteractorStyle()
  {
    selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    selectedActor = vtkSmartPointer<vtkActor>::New();
  }

  virtual void OnKeyPress() override
  //virtual void OnLeftButtonDown() override
  {
    vtkRenderWindowInteractor* rwi = this->Interactor;
    std::string key = rwi->GetKeySym();
    std::cout << "Pressed " << key << std::endl;
      
    if (key != "s") {
        return;
    }

    vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();

    // Get the location of the click (in window coordinates)
    int* pos = this->GetInteractor()->GetEventPosition();

    vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(0.0005);

    // Pick from this location.
    picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

    double* worldPosition = picker->GetPickPosition();
    std::cout << "Cell id is: " << picker->GetCellId() << std::endl;

    if (picker->GetCellId() != -1)
    {

      std::cout << "Pick position is: " << worldPosition[0] << " "
                << worldPosition[1] << " " << worldPosition[2] << endl;

      vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
      ids->SetNumberOfComponents(1);
      ids->InsertNextValue(picker->GetCellId());

      vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
      selectionNode->SetFieldType(vtkSelectionNode::CELL);
      selectionNode->SetContentType(vtkSelectionNode::INDICES);
      selectionNode->SetSelectionList(ids);

      vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
      selection->AddNode(selectionNode);

      vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
      extractSelection->SetInputData(0, this->Data);
      extractSelection->SetInputData(1, selection);
      extractSelection->Update();

      // In selection
      vtkSmartPointer<vtkUnstructuredGrid> selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
      selected->ShallowCopy(extractSelection->GetOutput());

      std::cout << "There are " << selected->GetNumberOfPoints() << " points in the selection." << std::endl;
      std::cout << "There are " << selected->GetNumberOfCells() << " cells in the selection." << std::endl;
      selectedMapper->SetInputData(selected);
      selectedMapper->ScalarVisibilityOff();
      selectedActor->SetMapper(selectedMapper);
      //selectedActor->GetProperty()->EdgeVisibilityOn();
      selectedActor->GetProperty()->SetRepresentationToWireframe();
      selectedActor->GetProperty()->SetColor(colors->GetColor3d("Green").GetData());
      selectedActor->GetProperty()->SetLineWidth(4);

      this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);
    }
    // Forward events
    //vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
  }

  vtkSmartPointer<vtkUnstructuredGrid> Data;
  //vtkSmartPointer<vtkPolyData> Data;
  vtkSmartPointer<vtkDataSetMapper> selectedMapper;
  vtkSmartPointer<vtkActor> selectedActor;
};

vtkStandardNewMacro(MouseInteractorStyle);

//------
// main
//------

int main(int argc, char* argv[])
{
  if(argc != 2) {
    std::cout << "Usage: " << argv[0] << " <FileName>.vtu" << std::endl;
    return EXIT_FAILURE;
  }

  // Read VTK unstructured mesh (.vtu) file.
  //
  std::string file_name = argv[1];
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(file_name.c_str());
  reader->Update();

  // Get the mesh. 
  auto mesh = reader->GetOutput();

  // Convert mesh to polydata.
  /*
  vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
  geometryFilter->SetInputData(mesh);
  geometryFilter->Update(); 
  vtkPolyData* meshPolyData = geometryFilter->GetOutput();

  vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
  triangleFilter->SetInputData(meshPolyData);
  triangleFilter->Update();
  */

  // Print mesh information.
  //
  vtkIdType numCells = mesh->GetNumberOfCells();
  vtkCellArray* cells = mesh->GetCells();
  std::cout << "Cells size " << cells->GetSize() << std::endl;
  vtkUnsignedCharArray* cellTypes = mesh->GetCellTypesArray();

  int numTri = 0;
  int numTet = 0;
  int numHex = 0;
  int numWedge = 0;
  int numQuad = 0;
  int numLine = 0;

  vtkGenericCell* cell = vtkGenericCell::New();
  // Comment out the following line to enable printing to cout.
  std::cout.setstate(std::ios_base::badbit);     

  for (vtkIdType cellId = 0; cellId < numCells; cellId++) {
    mesh->GetCell(cellId, cell);
    auto dim = cell->GetCellDimension();
    auto numPts = cell->GetNumberOfPoints();
    std::cout << "cell " << cellId << "  dim " << dim;
    std::cout << "  numPts " << numPts;
    std::cout << "  topo";
    switch (cellTypes->GetValue(cellId)) {
      case VTK_TETRA:
        std::cout << " tet ";
        numTet += 1;
      break;
      case VTK_HEXAHEDRON:
        std::cout << " hex ";
        numHex += 1;
      break;
      case VTK_WEDGE:
        std::cout << " wedge ";
        numWedge += 1;
      break;
      case VTK_TRIANGLE:
        std::cout << " tri ";
        numTri++;
      break;
      case VTK_QUAD:
        std::cout << " quad ";
        numQuad += 1;
      break;
      case VTK_VERTEX:
        std::cout << " vert ";
      break;
      case VTK_LINE:
        std::cout << " line ";
        numLine += 1;
      break;
      default:
          std::cout << " *** unknown *** '" << cellTypes->GetValue(cellId) << "'";
      break;
    }

    std::cout << "  conn: ";
    for (vtkIdType pointInd = 0; pointInd < numPts; ++pointInd) {
      auto id = cell->PointIds->GetId(pointInd);
      std::cout << id << " ";
    }
    std::cout << std::endl;
  }

  std::cout.clear();
  std::cout << std::endl;
  std::cout << "Number of cells " << numCells << std::endl;
  std::cout << "Number of hex cells " << numHex << std::endl;
  std::cout << "Number of tet cells " << numTet << std::endl;
  std::cout << "Number of wedge cells " << numWedge << std::endl;
  std::cout << "Number of tri cells " << numTri << std::endl;
  std::cout << "Number of quad cells " << numQuad << std::endl;
  std::cout << "Number of line cells " << numLine << std::endl;

  // Check element volumes.
  //
  vtkSmartPointer<vtkMeshQuality> qualityFilter = vtkSmartPointer<vtkMeshQuality>::New();
  qualityFilter->SetInputData(mesh);
  qualityFilter->SetTetQualityMeasureToVolume();
  qualityFilter->Update();

  vtkDataSet* qualityMesh = qualityFilter->GetOutput();
  vtkSmartPointer<vtkDoubleArray> qualityArray = 
      vtkDoubleArray::SafeDownCast(qualityMesh->GetCellData()->GetArray("Quality"));
  int numNegVolCells = 0;
  int numZeroVolCells = 0;
  int numZeroVolWedge = 0;
  std::vector<double> volumes;

  for(vtkIdType i = 0; i < qualityArray->GetNumberOfTuples(); i++) {
    mesh->GetCell(i, cell);
    auto dim = cell->GetCellDimension();
    if (dim != 3) { 
      continue;
    }
    double val = qualityArray->GetValue(i);
    volumes.push_back(val);
    //std::cout << "Val " << val << std::endl;

    if (val < 0.0) {
      numNegVolCells += 1;
    } else if (val == 0.0) {
      numZeroVolCells += 1;
      if (cellTypes->GetValue(i) == VTK_WEDGE) {
        numZeroVolWedge += 1;
      }
    }
  }

  std::cout << "Number of cells with negative volume " << numNegVolCells << std::endl;
  std::cout << "Number of cells with zero volume " << numZeroVolCells << std::endl;
  std::cout << "Number of wedges with zero volume " << numZeroVolWedge << std::endl;

  std::sort(volumes.begin(), volumes.end());
  std::cout << "Minimum volume: " << volumes[0] << std::endl;
  std::cout << "Maximum volume: " << volumes.back() << std::endl;

  // Select cells with volume <= 'lower'.
  //
  double lower = 1.0e-6;
  vtkSmartPointer<vtkThreshold> selectCells = vtkSmartPointer<vtkThreshold>::New();
  selectCells->ThresholdByLower(lower);
  selectCells->SetInputArrayToProcess( 0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, 
      vtkDataSetAttributes::SCALARS );
  selectCells->SetInputData(qualityMesh);
  selectCells->Update();
  auto filteredMesh = selectCells->GetOutput();
  std::cout << "Number of cells selected for volume <= " << lower << ": " <<filteredMesh->GetNumberOfCells() << std::endl;

  // Visualize bad cells.
  //
  vtkSmartPointer<vtkDataSetMapper> selectMapper = vtkSmartPointer<vtkDataSetMapper>::New();
  selectMapper->SetInputData(filteredMesh);
  selectMapper->ScalarVisibilityOff();
  vtkSmartPointer<vtkActor> selectActor = vtkSmartPointer<vtkActor>::New();
  selectActor->SetMapper(selectMapper);
  selectActor->GetProperty()->SetColor(0.0, 1.0, 0.0);
  selectActor->GetProperty()->SetEdgeColor(0.0, 1.0, 0.0); 
  selectActor->GetProperty()->EdgeVisibilityOn();
  selectActor->GetProperty()->SetRepresentationToWireframe();

  // Visualize entire mesh.
  //
  // Shrink the mesh.
  auto shrinkFactor = 0.9;
  vtkSmartPointer<vtkShrinkFilter> shrinkFilter = vtkSmartPointer<vtkShrinkFilter>::New();
  shrinkFilter->SetInputData(mesh);
  shrinkFilter->SetShrinkFactor(shrinkFactor);
  shrinkFilter->Update();

  // Create graphics geometry.
  vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
  mapper->SetInputData(shrinkFilter->GetOutput());
  //mapper->SetInputData(mesh);
  //mapper->SetInputConnection(triangleFilter->GetOutputPort());
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  mapper->ScalarVisibilityOff();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
  //actor->GetProperty()->SetRepresentationToWireframe();
  //actor->GetProperty()->EdgeVisibilityOn();
  actor->GetProperty()->BackfaceCullingOn();

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  //renderer->AddActor(selectActor);
  renderer->AddActor(actor);
  renderer->SetBackground(1.0, 1.0, 1.0); 

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  renderWindow->SetSize(500, 500);

  // Add window interactor.
  //
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  // Trackball interactor.
  /*
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle(style);
  */

  // Cell picking interactor.
  vtkSmartPointer<MouseInteractorStyle> style = vtkSmartPointer<MouseInteractorStyle>::New();
  style->SetDefaultRenderer(renderer);
  //style->Data = reader->GetOutput();
  //style->Data = triangleFilter->GetOutput();
  //style->Data = mesh;
  style->Data = shrinkFilter->GetOutput();

  renderWindowInteractor->SetInteractorStyle(style);


  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
