import vtk
import os

ARC_RESOLUTION = 8

"""
#Write facet objects as vtk files
def writeFacets(facets, savename):
    try:
        os.mkdir("advection_vtk/{}".format(savename))
    except:
        pass

    g = open("advection_vtk/{}_all.visit".format(savename), "w")

    facetnum = 0
    for facet in facets:
        if facet.name == 'linear':
            line = vtk.vtkLineSource()
            line.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
            line.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
            writer.SetInputConnection(line.GetOutputPort())
            writer.Write()
            facetnum += 1
        elif facet.name == 'arc':
            arc = vtk.vtkArcSource()
            arc.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
            arc.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
            arc.SetCenter(facet.center[0], facet.center[1], 0)
            arc.SetResolution(8)
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
            writer.SetInputConnection(arc.GetOutputPort())
            writer.Write()
            facetnum += 1
        elif facet.name == 'corner':
            if facet.radiusLeft is not None:
                arc = vtk.vtkArcSource()
                arc.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
                arc.SetPoint2(facet.corner[0], facet.corner[1], 0)
                arc.SetCenter(facet.centerLeft[0], facet.centerLeft[1], 0)
                arc.SetResolution(8)
                writer = vtk.vtkXMLPolyDataWriter()
                writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                writer.SetInputConnection(arc.GetOutputPort())
            else:
                line = vtk.vtkLineSource()
                line.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
                line.SetPoint2(facet.corner[0], facet.corner[1], 0)
                writer = vtk.vtkXMLPolyDataWriter()
                writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                writer.SetInputConnection(line.GetOutputPort())
            writer.Write()
            facetnum += 1
            if facet.radiusRight is not None:
                arc = vtk.vtkArcSource()
                arc.SetPoint1(facet.corner[0], facet.corner[1], 0)
                arc.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
                arc.SetCenter(facet.centerRight[0], facet.centerRight[1], 0)
                arc.SetResolution(8)
                writer = vtk.vtkXMLPolyDataWriter()
                writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                writer.SetInputConnection(arc.GetOutputPort())
            else:
                line = vtk.vtkLineSource()
                line.SetPoint1(facet.corner[0], facet.corner[1], 0)
                line.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
                writer = vtk.vtkXMLPolyDataWriter()
                writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                writer.SetInputConnection(line.GetOutputPort())
            writer.Write()
            facetnum += 1

    g.write("!NBLOCKS {}\n".format(facetnum))
    for i in range(facetnum):
        g.write("{}/facet_{}.vtp\n".format(savename, i))
    g.close()
"""

def writeFacets(facets, path):
    base_path = '/'.join(path.split('/')[:-1])
    if not os.path.exists(base_path):
        os.makedirs(base_path, exist_ok=True)

    vtkappend = vtk.vtkAppendPolyData()

    for facet in facets:
        if facet.name == 'linear':
            line = vtk.vtkLineSource()
            line.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
            line.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
            line.Update()
            vtkappend.AddInputData(line.GetOutput())
        elif facet.name == 'arc':
            arc = vtk.vtkArcSource()
            arc.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
            arc.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
            arc.SetCenter(facet.center[0], facet.center[1], 0)
            arc.SetResolution(ARC_RESOLUTION)
            arc.Update()
            vtkappend.AddInputData(arc.GetOutput())
        elif facet.name == 'corner':
            # Left
            if facet.centerLeft is None and facet.radiusLeft is None:
                line = vtk.vtkLineSource()
                line.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
                line.SetPoint2(facet.corner[0], facet.corner[1], 0)
                line.Update()
                vtkappend.AddInputData(line.GetOutput())
            else:
                arc = vtk.vtkArcSource()
                arc.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
                arc.SetPoint2(facet.corner[0], facet.corner[1], 0)
                arc.SetCenter(facet.centerLeft[0], facet.centerLeft[1], 0)
                arc.SetResolution(ARC_RESOLUTION)
                arc.Update()
                vtkappend.AddInputData(arc.GetOutput())
            # Right
            if facet.centerRight is None and facet.radiusRight is None:
                line = vtk.vtkLineSource()
                line.SetPoint1(facet.corner[0], facet.corner[1], 0)
                line.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
                line.Update()
                vtkappend.AddInputData(line.GetOutput())
            else:
                arc = vtk.vtkArcSource()
                arc.SetPoint1(facet.corner[0], facet.corner[1], 0)
                arc.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
                arc.SetCenter(facet.centerRight[0], facet.centerRight[1], 0)
                arc.SetResolution(ARC_RESOLUTION)
                arc.Update()
                vtkappend.AddInputData(arc.GetOutput())
        else:
            print(f"Unknown facet type: {facet.name}")
    
    vtkappend.Update()
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(path)
    writer.SetInputConnection(vtkappend.GetOutputPort())
    writer.Update()
    writer.Write()