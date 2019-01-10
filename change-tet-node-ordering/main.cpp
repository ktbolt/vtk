
// This is an exmaple of changing the node ordering for tetrahedral elements in an unstructured mesh. 
//
//

#include <string>

#include <vtkActor.h>
#include <vtkCellData.h>
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellIterator.h"
#include "vtkCellType.h"
#include "vtkCellLinks.h"
#include "vtkGenericCell.h"
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>

#include <vtkMeshQuality.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkThreshold.h>
#include "vtkUnsignedCharArray.h"
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

int main(int argc, char* argv[])
{
  if (argc != 2) {
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

  // Print mesh information.
  //
  vtkIdType numCells = mesh->GetNumberOfCells();
  vtkCellArray* cells = mesh->GetCells();
  vtkUnsignedCharArray* cellTypes = mesh->GetCellTypesArray();

  int numTri = 0;
  int numTet = 0;
  int numHex = 0;
  int numWedge = 0;
  int numQuad = 0;
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

  // Change node ordering for tets.
  //
  int numOrderChanged = 0;
  vtkIdType ids[6];
  vtkIdType loc = 0;
  for (vtkIdType cellId = 0; cellId < numCells; cellId++) {
    mesh->GetCell(cellId, cell);
    auto dim = cell->GetCellDimension();
    auto numPts = cell->GetNumberOfPoints();
    int i = 0;
    for (vtkIdType pointInd = 0; pointInd < numPts; ++pointInd) {
      auto id = cell->PointIds->GetId(pointInd);
      ids[i++] = id;
    }
    if ((dim == 3) && (numPts == 4)) {
      auto id0 = ids[0];
      ids[0] = ids[1];
      ids[1] = id0;
      cells->ReplaceCell(loc, numPts, ids);
      numOrderChanged += 1;
    }
    loc += numPts + 1;
  }

  std::cout << std::endl;
  std::cout << "Number of cells with changed order " << numOrderChanged << std::endl;
  mesh->Modified();

  // Check element volumes.
  //
  double lower = 0.0;
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

  for(vtkIdType i = 0; i < qualityArray->GetNumberOfTuples(); i++) {
    mesh->GetCell(i, cell);
    auto dim = cell->GetCellDimension();
    if (dim != 3) { 
      continue;
    }
    double val = qualityArray->GetValue(i);
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

  // Select bad cells.
  //
  vtkSmartPointer<vtkThreshold> selectCells = vtkSmartPointer<vtkThreshold>::New();
  selectCells->ThresholdByLower(lower);
  selectCells->SetInputArrayToProcess( 0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, 
      vtkDataSetAttributes::SCALARS );
  selectCells->SetInputData(qualityMesh);
  selectCells->Update();
  auto filteredMesh = selectCells->GetOutput();
  std::cout << "Number of cells selected for volume <= 0.0 " << filteredMesh->GetNumberOfCells() << std::endl;

  // Visualize bad cells using wireframe.
  //
  vtkSmartPointer<vtkDataSetMapper> selectMapper = vtkSmartPointer<vtkDataSetMapper>::New();
  selectMapper->SetInputData(filteredMesh);
  vtkSmartPointer<vtkActor> selectActor = vtkSmartPointer<vtkActor>::New();
  selectActor->SetMapper(selectMapper);
  selectActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
  selectActor->GetProperty()->SetEdgeColor(1.0, 0.0, 0.0); 
  selectActor->GetProperty()->EdgeVisibilityOn();
  selectActor->GetProperty()->SetRepresentationToWireframe();

  vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
  mapper->SetInputData(mesh);
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(0.0, 1.0, 1.0);

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->AddActor(selectActor);
  renderer->AddActor(actor);
  renderer->SetBackground(1.0, 1.0, 1.0); 

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
