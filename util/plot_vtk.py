import os
import vtk

from main.structs.meshes.base_mesh import BaseMesh

# Plot mesh as vtk file
def writeMesh(m: BaseMesh, path):
    base_path = '/'.join(path.split('/')[:-1])
    if not os.path.exists(base_path):
        os.makedirs(base_path, exist_ok=True)

    sgrid = vtk.vtkStructuredGrid()
    sgrid.SetDimensions([len(m._points), len(m._points[0]), 1])
    vtkpoints = vtk.vtkPoints()
    vtkpoints.Allocate(len(m._points)*len(m._points[0]))
    counter = 0
    for x in range(len(m._points)):
        for y in range(len(m._points[0])):
            vtkpoints.InsertPoint(counter, [m._points[x][y][0], m._points[x][y][1], 0])
            counter += 1
    sgrid.SetPoints(vtkpoints)
    writer = vtk.vtkStructuredGridWriter()
    writer.SetFileName(path)
    writer.SetInputData(sgrid)
    writer.Write()
    