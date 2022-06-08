from main.structs.meshes.base_mesh import BaseMesh
from main.structs.polys.neighbored_polygon import NeighboredPolygon

from main.algos.merge.greedy import greedyMerge

class PolyWithNeighbors:
    def __init__(self, poly, neighbors):
        self.

class MergeMesh(BaseMesh):

    def __init__(self, points, areas=None):
        super().__init__(self, points, areas)

        #List of merged polys, index in list is used as id
        self.merged_polys = []

        #New id
        self.merged_id = 0

        #Same shape as self.polys
        self.cartesian_to_merged_id = [[None] * (len(self.polys[0])) for _ in range(len(self.polys))]

    #TODO: figure out how to load algo by name
    def runMergeAlgo()