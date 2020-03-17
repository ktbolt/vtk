'''
This script reads in an unstrcutured mesh and prints the node coordinates and element connectivity.
'''
import os
import sys
import vtk

if __name__ == '__main__':

    ## Read mesh.
    #
    file_name = sys.argv[1]
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()
    mesh = reader.GetOutput()

    ## Print node coordinates.
    num_points = mesh.GetNumberOfPoints()
    points = mesh.GetPoints()
    node_ids = mesh.GetPointData().GetArray('GlobalNodeID')
    print("Number of nodes: {0:d}".format(num_points))
    print("Nodes:") 
    pt = 3*[0.0]
    min_point_id = 1e6 
    for i in range(num_points):
        points.GetPoint(i, pt)
        x = pt[0]
        y = pt[1]
        z = pt[2]
        node_id = node_ids.GetValue(i)
        print("Node {0:d}  ID: {1:d}  Coords: {2:g} {3:g} {4:g}".format(i, node_id, x, y, x))
        if node_id < min_point_id:
            min_point_id = node_id

    ## Print element connectivity.
    num_cells = mesh.GetNumberOfCells()
    element_ids = mesh.GetCellData().GetArray('GlobalElementID')
    print(" ")
    print("Number of elements: {0:d}".format(num_cells))
    cells = mesh.GetCells()
    min_cell_id = 1e6
    for i in range(num_cells):
        elem_id = element_ids.GetValue(i)
        cell = mesh.GetCell(i)
        num_cell_nodes = cell.GetNumberOfPoints()
        cell_nodes = cell.GetPoints()
        cell_node_ids = []
        for j in range(0, num_cell_nodes):
            pid = cell.GetPointId(j)
            node_id = node_ids.GetValue(pid)
            if node_id < min_cell_id:
                min_cell_id = node_id
            cell_node_ids.append(node_id)
        print("Element {0:d}  ID: {1:d}  Node IDs: {2:s}".format(i, elem_id, ' '.join(map(str,cell_node_ids))))

    print("min_point_id: " + str(min_point_id))
    print("min_cell_id: " + str(min_cell_id))


