import math
from main.geoms.geoms import getDistance, getArea, mergePolys, getPolyIntersectArea, getPolyLineArea, getPolyLineIntersects, lineIntersect, getCentroid
from main.geoms.linear_facet import getLinearFacet
from main.geoms.circular_facet import getCircleIntersectArea, getCircleCircleIntersects, getArcFacet, getArcFacetNewton, getCircleLineIntersects2, getCenter
from main.geoms.corner_facet import getPolyCornerArea, getPolyCurvedCornerArea, getCurvedCornerFacet
import numpy as np
from scipy.integrate import dblquad

from main.structs.polys.base_polygon import BasePolygon
from main.structs.meshes.base_mesh import BaseMesh

#m = BaseMesh, return areas array
def initializeCircle(m: BaseMesh, center, radius):
    areas = [[0] * len(m.polys[0]) for _ in range(len(m.polys))]

    for x in range(len(areas)):
        for y in range(len(areas[0])):
            poly: BasePolygon = m.polys[x][y]
            intersectarea, _ = getCircleIntersectArea(center, radius, poly.points)

            areas[x][y] = intersectarea/poly.getMaxArea()

    return areas

def initializePoly(m: BaseMesh, poly):
    areas = [[0] * len(m.polys[0]) for _ in range(len(m.polys))]

    for x in range(len(areas)):
        for y in range(len(areas[0])):
            mesh_poly: BasePolygon = m.polys[x][y]
            polyintersects = getPolyIntersectArea(poly, mesh_poly.points)
            for polyintersect in polyintersects:
                areas[x][y] += getArea(polyintersect)
            areas[x][y] /= mesh_poly.getMaxArea()

    return areas

#theta is angle from positive x-axis to major axis
def initializeEllipse(m: BaseMesh, major_axis, minor_axis, theta, center):
    areas = [[0] * len(m.polys[0]) for _ in range(len(m.polys))]

    circle_to_ellipse = np.array([[major_axis*math.cos(theta)**2 + minor_axis*math.sin(theta)**2, (major_axis-minor_axis)*math.cos(theta)*math.sin(theta)], [(major_axis-minor_axis)*math.cos(theta)*math.sin(theta), major_axis*math.sin(theta)**2 + minor_axis*math.cos(theta)**2]])
    ellipse_to_circle = np.linalg.inv(circle_to_ellipse)

    for x in range(len(areas)):
        for y in range(len(areas[0])):
            poly: BasePolygon = m.polys[x][y]
            centered_opoly = list(map(lambda x : [x[0]-center[0], x[1]-center[1]], poly.points))
            squished_opoly = list(map(lambda x : [ellipse_to_circle[0][0]*x[0] + ellipse_to_circle[0][1]*x[1], ellipse_to_circle[1][0]*x[0] + ellipse_to_circle[1][1]*x[1]], centered_opoly))
            intersectarea, _ = getCircleIntersectArea([0, 0], 1, squished_opoly)
            areas[x][y] = intersectarea*major_axis*minor_axis
            areas[x][y] /= poly.getMaxArea()

    return areas

#Example from Zalesak: circle with rectangle removed, use with 100x100 grid
def zalesak(m: BaseMesh):
    areas = [[0] * len(m.polys[0]) for _ in range(len(m.polys))]

    for x in range(len(areas)):
        for y in range(len(areas[0])):
            poly = m.polys[x][y]

            radiussmall = 15
            center = [50.05, 75.05]
                
            area, intersect = getCircleIntersectArea(center, radiussmall, poly.points)
            areas[x][y] += area
            
            rectangle = [[47.55, 59.5], [52.55, 59.5], [52.55, 85.05], [47.55, 85.05]]
            intersects = getPolyIntersectArea(rectangle, poly.points)
            for intersect in intersects:
                areas[x][y] -= getArea(intersect)
            areas[x][y] = max(0, areas[x][y])

            if getDistance(poly.points[0], [47.55, 60.25980054160715]) < math.sqrt(2):
                areas[x][y] = getPolyCurvedCornerArea(poly.points, [35.05, 75.05], [47.55, 60.25980054160715], [47.55, 85.05], 15, None)
            elif getDistance(poly.points[0], [52.55, 60.25980054160715]) < math.sqrt(2):
                areas[x][y] = getPolyCurvedCornerArea(poly.points, [52.55, 85.05], [52.55, 60.25980054160715], [65.05, 75.05], None, 15)

            areas[x][y] /= poly.getMaxArea()

    return areas

#x+o example: use with 100x100 grid
def xpluso(m: BaseMesh):
    areas = [[0] * len(m.polys[0]) for _ in range(len(m.polys))]

    for x in range(len(areas)):
        for y in range(len(areas[0])):
            poly = m.polys[x][y]

            #material 1 cross
            xpoints = [[4, 0], [10, 6], [16, 0], [20, 4], [14, 10], [20, 16], [16, 20], [10, 14], [4, 20], [0, 16], [6, 10], [0, 4]]
            xpoints = list(map(lambda x : [x[0]+3.005, x[1]+3.025], xpoints))
            xpolyintersects = getPolyIntersectArea(xpoints, poly.points)
            for xpolyintersect in xpolyintersects:
                areas[x][y] += abs(getArea(xpolyintersect))

            #material 2 cross
            xpoints = [[6, 0], [14, 0], [14, 6], [20, 6], [20, 14], [14, 14], [14, 20], [6, 20], [6, 14], [0, 14], [0, 6], [6, 6]]
            xpoints = list(map(lambda x : [x[0]+3.005, x[1]+28.025], xpoints))
            xpolyintersects = getPolyIntersectArea(xpoints, poly.points)
            for xpolyintersect in xpolyintersects:
                areas[x][y] += abs(getArea(xpolyintersect))

            #material 1 ring
            radius = 10
            radiussmall = 6
            center = [38.005, 13.005]
            area, _ = getCircleIntersectArea(center, radius, poly.points)
            areas[x][y] += area
            area, _ = getCircleIntersectArea(center, radiussmall, poly.points)
            areas[x][y] -= area

            #material 2 ring
            radius = 10
            radiussmall = 6
            center = [38.005, 38.005]
            area, _ = getCircleIntersectArea(center, radius, poly.points)
            areas[x][y] += area
            area, _ = getCircleIntersectArea(center, radiussmall, poly.points)
            areas[x][y] -= area

            #material 2 dot
            radius = 3
            center = [38.005, 13.005]
            area, _ = getCircleIntersectArea(center, radius, poly.points)
            areas[x][y] += area

            areas[x][y] /= poly.getMaxArea()

    return areas

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np

def plotAreas(opolys, areas):
    listpatches = []
    listareas = []
    listpartials = []
    for x in range(len(opolys)):
        for y in range(len(opolys[0])):
            opoly = opolys[x][y]
            opolyarea = getArea(opoly)
            #area_fraction = areas[x][y]/opolyarea
            area_fraction = areas[x][y]/opolyarea
            patch = Polygon(np.array(opoly), True)
            listpatches.append(patch)
            listareas.append(area_fraction)
            listpartials.append(math.ceil(area_fraction) - math.floor(area_fraction))

    #Max x and y of opolys grid
    maxX = max(list(map(lambda x : x[0], opolys[len(opolys)-1][len(opolys[0])-1])))
    maxY = max(list(map(lambda x : x[1], opolys[len(opolys)-1][len(opolys[0])-1])))

    p = PatchCollection(listpatches, cmap='jet')
    patchareas = np.array(listareas)
    p.set_array(patchareas)
    fig, ax = plt.subplots()
    ax.set_xlim(0, maxX)
    ax.set_ylim(0, maxY)
    ax.add_collection(p)
    plt.savefig("advection_plt/original_areas.png", dpi=199)
    plt.clf()

    p = PatchCollection(listpatches, cmap='jet')
    patchareas = np.array(listpartials)
    p.set_array(patchareas)
    fig, ax = plt.subplots()
    ax.set_xlim(0, maxX)
    ax.set_ylim(0, maxY)
    ax.add_collection(p)
    plt.savefig("advection_plt/original_partials.png", dpi=199)
    plt.clf()