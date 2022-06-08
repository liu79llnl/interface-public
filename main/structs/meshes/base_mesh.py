import numpy as np
import math
import vtk
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as plt_polygon
from matplotlib.collections import PatchCollection

from main.structs.polys.base_polygon import BasePolygon

class BaseMesh:

    #Parameters that can be accessed publicly: self.polys
    def __init__(self, points, fractions=None):
        #Points = [[x, y]]
        self._points = points

        #Compute extremes x and y
        self._max_x = -float('inf')
        self._max_y = -float('inf')
        self._min_x = float('inf')
        self._min_y = float('inf')

        for x in range(len(self._points)):
            for y in range(len(self._points[0])):
                #Set extremes
                if points[x][y][0] > self._max_x:
                    self._max_x = points[x][y][0]
                if points[x][y][1] > self._max_y:
                    self._max_y = points[x][y][1]
                if points[x][y][0] < self._min_x:
                    self._min_x = points[x][y][0]
                if points[x][y][1] < self._min_y:
                    self._min_y = points[x][y][1]

        self.polys = [[None] * (len(points[0])-1) for _ in range(len(points)-1)]
        self._plt_patches = []
        for x in range(len(self.polys)):
            for y in range(len(self.polys[0])):
                #Make quads
                poly = BasePolygon([points[x][y], points[x+1][y], points[x+1][y+1], points[x][y+1]])
                self.polys[x][y] = poly
                patch = plt_polygon(np.array(poly), True)
                self._plt_patches.append(patch)
        
        if fractions is None:
            self._plt_patchareas = np.array([])
            self._plt_patchpartialareas = np.array([])
            self._plt_patchinitialareas = np.array([])
        else:
            self.initializeFractions(fractions)

    def initializeFractions(self, fractions):
        patchinitialareas = []

        #Flatten array
        for x in range(len(self.polys)):
            for y in range(len(self.polys[0])):
                patchinitialareas.append(fractions[x][y])

        self._plt_patchinitialareas = np.array(patchinitialareas)
        self.setFractions(fractions)

    #TODO: make sure these are actually areas and not fractions
    def setFractions(self, fractions):
        patchareas = []
        patchpartialareas = []

        #Flatten array
        for x in range(len(self.polys)):
            for y in range(len(self.polys[0])):
                self.polys.setFraction(fractions[x][y]) #TODO double-check here
                a = fractions[x][y]
                patchareas.append(a)
                patchpartialareas.append(math.ceil(a - math.floor(a)))

        self._plt_patchareas = np.array(patchareas)
        self._plt_patchpartialareas = np.array(patchpartialareas)

    #Plot mesh and areas as plt images
    def plotPolyValues(self, values, path):
        patchcollection = PatchCollection(self._plt_patches, cmap='plasma') #jet
        patchcollection.set_array(values)
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.set_xlim(self._min_x, self._max_x)
        ax.set_ylim(self._min_y, self._max_y)
        ax.add_collection(patchcollection)
        plt.savefig(path, dpi=199)
        plt.close()
        plt.clf()

    def plotAreas(self, path):
        self.plotPolyValues(self._plt_patchareas, path)

    def plotPartialAreas(self, path):
        self.plotPolyValues(self._plt_patchpartialareas, path)

    def plotInitialAreaCompare(self, path):
        self.plotPolyValues(self._plt_patchareas - self._plt_patchinitialareas, path)

    #Plot mesh as vtk file
    def plotMesh(self, path):
        sgrid = vtk.vtkStructuredGrid()
        sgrid.SetDimensions([len(self._points), len(self._points[0]), 1])
        vtkpoints = vtk.vtkPoints()
        vtkpoints.Allocate(len(self._points)*len(self._points[0]))
        counter = 0
        for x in range(len(self._points)):
            for y in range(len(self._points[0])):
                vtkpoints.InsertPoint(counter, [self._points[x][y][0], self._points[x][y][1], 0])
                counter += 1
        sgrid.SetPoints(vtkpoints)
        writer = vtk.vtkStructuredGridWriter()
        writer.SetFileName(path)
        writer.SetInputData(sgrid)
        writer.Write()