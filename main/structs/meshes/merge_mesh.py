from concurrent.futures import process
from main.structs.meshes.base_mesh import BaseMesh
from main.structs.polys.base_polygon import BasePolygon
from main.structs.polys.neighbored_polygon import NeighboredPolygon
from main.geoms.geoms import mergePolys, getPolyIntersectArea, getArea, getPolyLineArea
from main.geoms.circular_facet import getCircleIntersectArea
from main.structs.facets.base_facet import advectPoint

from typing import Dict

class MergeMesh(BaseMesh):

    def __init__(self, points, threshold, areas=None):
        super().__init__(points, threshold, areas)

        # List of merged polys, index in list is used as id.
        # Each element is a list of (x, y)s corresponding to the Cartesian coordinates of the polys to be merged.
        self.merge_ids_to_coords = []
        # List of neighbor ids, index matches self.merge_ids_to_coords
        self.merge_id_to_neighbor_ids = []

        # Newest id that can be used
        self.next_merge_id = 0

        # Same shape as self.polys
        self.coords_to_merge_id = [[None for _ in range(len(self.polys[0]))] for _ in range(len(self.polys))]

        # Dict of NeighboredPolygon objects, index matches self.merge_ids_to_coords
        self.merged_polys = dict()

    def _get_merge_id(self, x, y):
        return self.coords_to_merge_id[x][y]

    def _get_merge_coords(self, merge_id):
        return self.merge_ids_to_coords[merge_id]

    def _get_neighbor_ids_from_merge_id(self, merge_id):
        return self.merge_id_to_neighbor_ids[merge_id]

    def _get_neighbor_ids_from_coords(self, x, y):
        return self._get_neighbor_ids_from_merge_id(self._get_merge_id(x, y))

    def _get_num_merge_ids(self):
        return self.next_merge_id

    # Merge the polys corresponding to merge_ids
    # merge_ids = list of merge_ids
    def _merge(self, merge_ids):
        # Get coords of polys to be merged
        # Get new neighbor ids
        merge_coords = []
        neighbor_ids = set()
        for merge_id in merge_ids:
            merge_coords += self._get_merge_coords(merge_id)
            some_neighbor_ids = self._get_neighbor_ids_from_merge_id(merge_id)
            # For each neighbor
            for neighbor_id in some_neighbor_ids:
                if neighbor_id not in merge_ids:
                    neighbor_ids.add(neighbor_id)
                    # Replace old merge id with new one
                    self.merge_id_to_neighbor_ids[neighbor_id] = [(self.next_merge_id if x == merge_id else x) for x in self.merge_id_to_neighbor_ids[neighbor_id]]

        self.merge_ids_to_coords.append(merge_coords)
        self.merge_id_to_neighbor_ids.append(list(neighbor_ids))

        for merge_coord in merge_coords:
            [x, y] = merge_coord
            self.coords_to_merge_id[x][y] = self.next_merge_id
            
        self.next_merge_id += 1

    # Each time new fractions are set, run merging algorithm
    def setFractions(self, fractions):
        super().setFractions(fractions)
        
        def _helper_isMixed(x, y):
            if x < 0 or y < 0 or x >= len(self.polys) or y >= len(self.polys[0]):
                return False
            return self.polys[x][y].isMixed()

        #Set global values to initials
        self.merge_ids_to_coords = []
        self.merge_id_to_neighbor_ids = []
        self.next_merge_id = 0
        self.coords_to_merge_id = [[None for _ in range(len(self.polys[0]))] for _ in range(len(self.polys))]
        self.merged_polys = dict()

        # Set each individual poly to its own merge id
        for x in range(len(self.polys)):
            for y in range(len(self.polys[0])):
                if self.polys[x][y].isMixed():
                    self.merge_ids_to_coords.append([[x, y]])
                    self.merge_id_to_neighbor_ids.append([])
                    self.coords_to_merge_id[x][y] = self.next_merge_id
                    self.next_merge_id += 1

        dirs = [[1, 0], [0, 1], [-1, 0], [0, -1]]

        # Locate mixed neighbors and set neighbor ids
        for x in range(len(self.polys)):
            for y in range(len(self.polys[0])):
                if self.polys[x][y].isMixed():
                    merge_id = self._get_merge_id(x, y)
                    for dir in dirs:
                        neighbor_coords = [x+dir[0], y+dir[1]]
                        if _helper_isMixed(neighbor_coords[0], neighbor_coords[1]):
                            neighbor_id = self._get_merge_id(neighbor_coords[0], neighbor_coords[1])
                            self.merge_id_to_neighbor_ids[merge_id].append(neighbor_id)

    # Take cell and return 3x3 array of area fractions (for Young's)
    def get3x3Stencil(self, x, y):
        def _helper_in_bounds(i, j):
            if i < 0 or j < 0 or i >= len(self.polys) or j >= len(self.polys[0]):
                return 0
            return self.polys[i][j].getFraction()

        ret = [[0,0,0],[0,0,0],[0,0,0]]
        for x_delta in range(-1, 2):
            for y_delta in range(-1, 2):
                ret[1+x_delta][1+y_delta] = _helper_in_bounds(x+x_delta, y+y_delta)

        return ret

    # Get full/empty cell coords whose fraction is closest to being mixed
    # If all are mixed, returns [None, None]
    def getMostMixedAdjacentFullCell(self, x, y):
        def _helper_diffs(i, j):
            # Ignore mixed cells
            if i < 0 or j < 0 or i >= len(self.polys) or j >= len(self.polys[0]) or self.polys[i][j].isMixed():
                return float("inf")
            return self.polys[i][j].diffFromMixed()
        
        dirs = [[1, 0], [0, 1], [-1, 0], [0, -1]]

        ret_x = None
        ret_y = None
        min_diff = float("inf")
        for dir in dirs:
            diff = _helper_diffs(x+dir[0], y+dir[1])
            if diff < min_diff:
                ret_x = x+dir[0]
                ret_y = y+dir[1]

        diffs = list(map(lambda dir : _helper_diffs(x+dir[0], y+dir[1]), dirs))
        print(diffs)
        return [ret_x, ret_y]


    #TODO add check for case of diagonal neighbors
    def merge1Neighbors(self):
        process_queue = []

        # Find all 1 neighbors
        for x in range(len(self.polys)):
            for y in range(len(self.polys[0])):
                merge_id = self._get_merge_id(x, y)
                if merge_id is not None:
                    neighbor_ids = self._get_neighbor_ids_from_merge_id(merge_id)
                    if len(neighbor_ids) == 1:
                        process_queue.append(merge_id)

        # Merge all 1 neighbors
        while process_queue:
            merge_id = process_queue.pop(0)
            [x, y] = self._get_merge_coords(merge_id)[0]
            neighbor_ids = self._get_neighbor_ids_from_merge_id(merge_id)
            # By here, this poly should only have one neighbor. If not, low resolution: case of a long chain of mixed cells?
            if len(neighbor_ids) != 1 or self._get_merge_id(x, y) != merge_id: # second case is when merge_id has already been merged somehow
                raise ValueError(f"Error in merge1Neighbors: neighbor_ids {neighbor_ids} has length != 1. Low resolution: long chain of mixed cells?")
            # Check if neighbor only has two neighbors: if so, after merging, merged poly would again have a single neighbor, which could be a long chain of mixed cells. We want to avoid this.
            neighbor_id = neighbor_ids[0]
            if len(self._get_neighbor_ids_from_merge_id(neighbor_id)) < 3:
                # Probably a long chain of mixed cells. Instead of merging with neighbor cell, find adjacent empty/full cell that's closest to being mixed and merge with that instead.
                # For now, instead of trying to turn a full cell into a mixed cell, just skip it
                # merge_coords = self._get_merge_coords(merge_id)
                # # merge_coords should be shape [[x, y]]
                # [full_x, full_y] = self.getMostMixedAdjacentFullCell(merge_coords[0][0], merge_coords[0][1])
                # self._merge([merge_id, self._get_merge_id(full_x, full_y)])
                pass
            else:
                # No issue, merge with its single neighbor
                self._merge([merge_id, neighbor_ids[0]])

    def _helper_createMergedPolys(self, merge_coords_list):
        merge_coords = merge_coords_list.copy()
        if len(merge_coords) == 1:
            merge_points = self.polys[merge_coords[0][0]][merge_coords[0][1]].points
        else:
            merge_points = []
            i = 0
            while i < len(merge_coords):
                try:
                    merge_points = mergePolys(merge_points, self.polys[merge_coords[i][0]][merge_coords[i][1]].points)
                    merge_coords.pop(i)
                    i = 0
                except:
                    i += 1

            if len(merge_coords) > 0:
                raise ValueError(f"Error in _helper_createMergedPolys: number of polys to merge: {len(merge_coords)}")

        ret_poly = NeighboredPolygon(merge_points)
        total_area = sum(list(map(lambda x : self.polys[x[0]][x[1]].getArea(), merge_coords_list)))
        ret_poly.setArea(total_area)
        return ret_poly

    # Resets self.merged_polys according to the latest values in self.merge_ids_to_coords
    # Coordinates only, no neighbors
    def createMergedPolys(self):
        #Set global variables to initials
        self.merged_polys = dict()

        # Add all merge ids in use to the queue
        for x in range(len(self.polys)):
            for y in range(len(self.polys[0])):
                merge_id = self._get_merge_id(x, y)
                if merge_id is not None and merge_id not in self.merged_polys.keys():
                    merge_coords = self._get_merge_coords(merge_id).copy()
                    self.merged_polys[merge_id] = self._helper_createMergedPolys(merge_coords)

    def advectMergedFacets(self, velocity, t, dt, checkSize=2):

        advected_areas = [[0] * len(self.polys[0]) for _ in range(len(self.polys))]
        advected_facets = []
        #List to keep track of merge ids such that advected facet has been appended
        processed_merge_ids = []

        #Advect interface and intersect with neighbors
        for advectx in range(len(self.polys)):
            for advecty in range(len(self.polys[0])):
                if self.polys[advectx][advecty].getFraction() > self.threshold:
                    #shiftpoly
                    shiftpoly = list(map(lambda x : advectPoint(x, velocity, t, dt), self.polys[advectx][advecty].points))
                    shiftbounds = [min(list(map(lambda x : x[0], shiftpoly))), min(list(map(lambda x : x[1], shiftpoly))), max(list(map(lambda x : x[0], shiftpoly))), max(list(map(lambda x : x[1], shiftpoly)))]

                    for testx in range(-checkSize, checkSize+1):
                        for testy in range(-checkSize, checkSize+1):
                            checkx = advectx-testx
                            checky = advecty-testy
                            if checkx >= 0 and checkx < len(self.polys) and checky >= 0 and checky < len(self.polys[0]):
                                testpoly = self.polys[checkx][checky].points
                                testbounds = [min(list(map(lambda x : x[0], testpoly))), min(list(map(lambda x : x[1], testpoly))), max(list(map(lambda x : x[0], testpoly))), max(list(map(lambda x : x[1], testpoly)))]
                                if not(testbounds[2] <= shiftbounds[0] or shiftbounds[2] <= testbounds[0] or testbounds[3] <= shiftbounds[1] or shiftbounds[3] <= testbounds[1]):
                                    #bounding boxes intersect, could be nonzero intersection
                                    #TODO is this part still necessary?
                                    try:
                                        polyintersections = getPolyIntersectArea(shiftpoly, testpoly)
                                    except:
                                        print("Failed polyintersect: getPolyIntersectArea({}, {})".format(shiftpoly, testpoly))
                                        testpoly = list(map(lambda x : [x[0]+1e-13, x[1]+1e-13], testpoly))
                                        polyintersections = getPolyIntersectArea(shiftpoly, testpoly)
                                    if len(polyintersections) == 0:
                                        #No intersection
                                        continue
                                    #For each overlap region
                                    for polyintersection in polyintersections:
                                        advect_merge_id = self._get_merge_id(advectx, advecty)
                                        if advect_merge_id is None:
                                            #Full cell
                                            #TODO is this necessary?
                                            try:
                                                assert self.polys[advectx][advecty].getFraction() > 1-self.threshold
                                            except:
                                                print(self.polys[advectx][advecty].points)
                                                print(self.polys[advectx][advecty].getFraction())
                                                print(self.polys[advectx][advecty].getArea())
                                                print(1/0)
                                            advected_areas[checkx][checky] += abs(getArea(polyintersection))
                                        else:
                                            #Mixed cell with facet
                                            advectedfacet = self.merged_polys[advect_merge_id].getFacet()

                                            #Linear or arc
                                            advectedfacet = advectedfacet.advected(velocity, t, dt)

                                            #Add to return
                                            if advect_merge_id not in processed_merge_ids:
                                                processed_merge_ids.append(advect_merge_id)
                                                advected_facets.append(advectedfacet)
                                            if advectedfacet.name == 'linear':
                                                polyintersectionarea = getPolyLineArea(polyintersection, advectedfacet.pLeft, advectedfacet.pRight)
                                            elif advectedfacet.name == 'arc':
                                                polyintersectionarea, _ = getCircleIntersectArea(advectedfacet.center, advectedfacet.radius, polyintersection)
                                            advected_areas[checkx][checky] += polyintersectionarea #TODO: abs here?
                                            #TODO necessary?
                                            if polyintersectionarea < 0:
                                                print("Negative polyintersectionarea")

                                            #TODO handle corner case here
                                            """
                                                elif len(predfacets[advectx][advecty]) == 2:
                                                    #Corner
                                                    advectedfacet1 = predfacets[advectx][advecty][0].advected(velocity, t, dt)
                                                    advectedfacet1r = advectedfacet1.radius if advectedfacet1.name == 'arc' else None
                                                    plot_advectedfacets.append(advectedfacet1)
                                                    advectedfacet2 = predfacets[advectx][advecty][1].advected(velocity, t, dt)
                                                    advectedfacet2r = advectedfacet2.radius if advectedfacet2.name == 'arc' else None
                                                    plot_advectedfacets.append(advectedfacet2)
                                                    nareas[checkx][checky] += getPolyCurvedCornerArea(polyintersection, advectedfacet1.pLeft, advectedfacet1.pRight, advectedfacet2.pRight, advectedfacet1r, advectedfacet2r)
                                                else:
                                                    print("More than two facets in this cell?")
                                            """
                                            
        #Update areas
        for x in range(len(self.polys)):
            for y in range(len(self.polys[0])):
                advected_areas[x][y] = min(max(advected_areas[x][y]/self.polys[x][y].getMaxArea(), 0), 1)

        self.setFractions(advected_areas)
        return advected_facets
            
    #TODO: figure out how to load algo by name
    """
    input: MergeMesh m, side effects on m

    1. Merge single-neighbor mixed cells until all (merged) polys have 2+ neighbors.
    2. For each poly with two neighbors, assume they're the true neighbors. Try to use orientation to determine left or right.
    Also add poly as a neighbor of those neighbors.
    3. For polys with 3+ neighbors, check if potential neighbors already have two neighbors set. If so, remove as potential neighbor. Readd to queue.

    Unclear if this will run into issues, especially when the interface exits the same edge twice.

    Returns merge ids #TODO poly objects instead?
    """
    def findOrientations(self):

        class MergeIdWithNeighbors:

            def __init__(self, merge_id, neighbor_ids: list):
                self.merge_id = merge_id
                self.neighbor_ids = neighbor_ids.copy()
                self.left_neighbor_id = None
                self.right_neighbor_id = None

            def has_left(self):
                return (self.left_neighbor_id is not None)

            def has_right(self):
                return (self.right_neighbor_id is not None)

            def set_left(self, neighbor_id, set_neighbor=True):
                if neighbor_id not in self.neighbor_ids:
                    raise ValueError(f"Error in set_left of findOrientations: neighbor_id {neighbor_id} is not in {self.neighbor_ids}")
                if self.has_left():
                    raise ValueError(f"Error in set_left of findOrientations: left neighbor is already {self.left_neighbor_id}")
                self.neighbor_ids.remove(neighbor_id)
                self.left_neighbor_id = neighbor_id
                if self.fully_oriented():
                    for id in self.neighbor_ids:
                        self.remove_neighbor_id(id)
                if set_neighbor:
                    neighbor_obj = merge_id_to_obj[neighbor_id]
                    neighbor_obj.set_right(self.get_merge_id(), set_neighbor=False)

            def set_right(self, neighbor_id, set_neighbor=True):
                if neighbor_id not in self.neighbor_ids:
                    raise ValueError(f"Error in set_right of findOrientations: neighbor_id {neighbor_id} is not in {self.neighbor_ids}")
                if self.has_right():
                    raise ValueError(f"Error in set_right of findOrientations: right neighbor is already {self.right_neighbor_id}")
                self.neighbor_ids.remove(neighbor_id)
                self.right_neighbor_id = neighbor_id
                if self.fully_oriented():
                    for id in self.neighbor_ids:
                        self.remove_neighbor_id(id)
                if set_neighbor:
                    neighbor_obj = merge_id_to_obj[neighbor_id]
                    neighbor_obj.set_left(self.get_merge_id(), set_neighbor=False)
            
            def get_left(self):
                return self.left_neighbor_id

            def get_right(self):
                return self.right_neighbor_id

            # Unassigned neighbor ids
            def get_neighbor_ids(self):
                return self.neighbor_ids

            def get_merge_id(self):
                return self.merge_id

            # Always removes a whole neighbor interaction: when removing neighbor id, also removes its own merge id from neighbor
            def remove_neighbor_id(self, neighbor_id):
                neighbor_obj: MergeIdWithNeighbors = merge_id_to_obj[neighbor_id]
                if neighbor_id not in self.neighbor_ids:
                    raise ValueError(f"Error in remove_neighbor_id of findOrientations: neighbor_id {neighbor_id} is not in neighbors {self.neighbor_ids}")
                elif self.merge_id not in neighbor_obj.neighbor_ids:
                    raise ValueError(f"Error in remove_neighbor_id of findOrientations: merge_id {self.merge_id} is not in neighbor's neighbors {neighbor_obj.neighbor_ids}")
                self.neighbor_ids.remove(neighbor_id)
                neighbor_obj.neighbor_ids.remove(self.merge_id)

            def fully_oriented(self):
                return (self.has_left() and self.has_right())

            #TODO
            def __str__(self):
                return f"Merge id: {self.merge_id}\nLeft neighbor id: {self.left_neighbor_id}\nRight neighbor id: {self.right_neighbor_id}\nCurrent neighbor candidates: {self.neighbor_ids}\n"

        merge_id_to_obj: Dict[int, MergeIdWithNeighbors] = dict()
        process_queue = []
        processed_merge_ids = []

        # Add all merge ids in use to the queue
        for x in range(len(self.polys)):
            for y in range(len(self.polys[0])):
                merge_id = self._get_merge_id(x, y)
                if merge_id is not None and merge_id not in processed_merge_ids:
                    processed_merge_ids.append(merge_id)
                    neighbor_ids = self._get_neighbor_ids_from_merge_id(merge_id)
                    merge_id_with_neighbors = MergeIdWithNeighbors(merge_id, neighbor_ids)
                    merge_id_to_obj[merge_id] = merge_id_with_neighbors
                    process_queue.append(merge_id)

        # Add new MergeIdWithNeighbors object to merge_id_to_obj which corresponds to merging merge_id with neighbor_id
        def mergeObjs(merge_id, neighbor_id, use_neighbor_id_orientation=False):
            # If either id is not in merge_id_to_obj, throw error
            if merge_id not in merge_id_to_obj.keys():
                raise ValueError(f"{merge_id} merge id not in {merge_id_to_obj.keys()}")
            elif neighbor_id not in merge_id_to_obj.keys():
                raise ValueError(f"{neighbor_id} neighbor id not in {merge_id_to_obj.keys()}")

            merge_id_with_neighbors: MergeIdWithNeighbors = merge_id_to_obj[merge_id]
            neighbor_id_with_neighbors: MergeIdWithNeighbors = merge_id_to_obj[neighbor_id]
            # For each of merge_id's (original) neighbors, replace merge_id with self.next_merge_id
            for merge_neighbor_id in self._get_neighbor_ids_from_merge_id(merge_id):
                if merge_neighbor_id != neighbor_id:
                    neighbor_with_neighbors: MergeIdWithNeighbors = merge_id_to_obj[merge_neighbor_id]
                    if neighbor_with_neighbors.has_left() and neighbor_with_neighbors.get_left() == merge_id:
                        neighbor_with_neighbors.left_neighbor_id = self.next_merge_id
                    if neighbor_with_neighbors.has_right() and neighbor_with_neighbors.get_right() == merge_id:
                        neighbor_with_neighbors.right_neighbor_id = self.next_merge_id
                    for i in range(len(neighbor_with_neighbors.get_neighbor_ids())):
                        if neighbor_with_neighbors.get_neighbor_ids()[i] == merge_id:
                            neighbor_with_neighbors.get_neighbor_ids()[i] = self.next_merge_id
            # For each of neighbor_id's (original) neighbors, replace neighbor_id with self.next_merge_id
            for neighbor_neighbor_id in self._get_neighbor_ids_from_merge_id(neighbor_id):
                if neighbor_neighbor_id != merge_id:
                    neighbor_with_neighbors: MergeIdWithNeighbors = merge_id_to_obj[neighbor_neighbor_id]
                    if neighbor_with_neighbors.has_left() and neighbor_with_neighbors.get_left() == neighbor_id:
                        neighbor_with_neighbors.left_neighbor_id = self.next_merge_id
                    if neighbor_with_neighbors.has_right() and neighbor_with_neighbors.get_right() == neighbor_id:
                        neighbor_with_neighbors.right_neighbor_id = self.next_merge_id
                    for i in range(len(neighbor_with_neighbors.get_neighbor_ids())):
                        if neighbor_with_neighbors.get_neighbor_ids()[i] == neighbor_id:
                            neighbor_with_neighbors.get_neighbor_ids()[i] = self.next_merge_id

            # New list of neighbors is merge_id's neighbors + neighbor_id's neighbors
            # Reset left/right neighbors
            new_neighbor_ids = list(filter(lambda x : x != merge_id and x != neighbor_id, merge_id_with_neighbors.get_neighbor_ids() + neighbor_id_with_neighbors.get_neighbor_ids()))
            new_merge_id_with_neighbors = MergeIdWithNeighbors(self.next_merge_id, new_neighbor_ids)
            if use_neighbor_id_orientation:
                # Danger: potential to ruin neighbor operation symmetry
                if neighbor_id_with_neighbors.has_left() and neighbor_id_with_neighbors.get_left() != merge_id:
                    new_merge_id_with_neighbors.left_neighbor_id = neighbor_id_with_neighbors.get_left()
                if neighbor_id_with_neighbors.has_right() and neighbor_id_with_neighbors.get_right() != merge_id:
                    new_merge_id_with_neighbors.right_neighbor_id = neighbor_id_with_neighbors.get_right()
            merge_id_to_obj[self.next_merge_id] = new_merge_id_with_neighbors

            # Remove merge_id and neighbor_id from merge_id_to_obj
            merge_id_to_obj.pop(merge_id)
            merge_id_to_obj.pop(neighbor_id)
            # Merge polygons and original neighbors
            self._merge([merge_id, neighbor_id])

            return new_merge_id_with_neighbors

        def doGreedyOrientations():
            iters_without_progress = 0
            while process_queue:
                merge_id = process_queue.pop(0)
                merge_id_with_neighbors: MergeIdWithNeighbors = merge_id_to_obj[merge_id]
                # 1 neighbor
                if len(merge_id_with_neighbors.get_neighbor_ids()) == 1:
                    neighbor_id = merge_id_with_neighbors.get_neighbor_ids()[0]
                    # If merge poly only has one unassigned neighbor and only needs to assign one more neighbor, do it
                    if merge_id_with_neighbors.has_left() and not(merge_id_with_neighbors.has_right()) and not(merge_id_to_obj[neighbor_id].has_left()):
                        merge_id_with_neighbors.set_right(neighbor_id)
                        iters_without_progress = 0
                    elif merge_id_with_neighbors.has_right() and not(merge_id_with_neighbors.has_left()) and not(merge_id_to_obj[neighbor_id].has_right()):
                        merge_id_with_neighbors.set_left(neighbor_id)
                        iters_without_progress = 0
                    # There should never be a case when neither orientation is set but only one neighbor
                    elif not(merge_id_with_neighbors.has_left()) and not(merge_id_with_neighbors.has_right()):
                        print("Neither orientation is set but one neighbor: why did this happen?")
                        #TODO not sure how to handle this, merge into its neighbor for now
                        process_queue.append(merge_id)
                        iters_without_progress += 1
                        #raise ValueError("Neither orientation is set but one neighbor: why did this happen?")
                    # This case should be dealt with via implementation of remove_neighbor_id
                    elif merge_id_with_neighbors.fully_oriented():
                        print("Fully oriented but one neighbor: why did this happen?")
                        #TODO is this an actual issue? For now, assume this cell is ok and just remove from queue
                        # process_queue.append(merge_id_with_neighbors)
                        # iters_without_progress += 1
                        #raise ValueError("Fully oriented but one neighbor: why did this happen?")
                    else: # single neighbor but it's already got the neighbor we would want to assign merge_id to
                        print("Single neighbor but orientations are inconsistent, skip for now")
                        process_queue.append(merge_id)
                        iters_without_progress += 1
                # 2 neighbors
                elif len(merge_id_with_neighbors.get_neighbor_ids()) == 2:
                    if not(merge_id_with_neighbors.has_left()) and not(merge_id_with_neighbors.has_right()):
                        # If merge poly consists of only one poly and has exactly two neighbors
                        if len(self._get_merge_coords(merge_id)) == 1:
                            # Check if orientation is easy to figure out
                            [x, y] = self._get_merge_coords(merge_id)[0]
                            dirs = [[1, 0], [0, 1], [-1, 0], [0, -1]]
                            neighbor_dirs = []
                            nonneighbor_dirs = []
                            def _helper_in_bounds(x, y):
                                return not(x < 0 or y < 0 or x >= len(self.polys) or y >= len(self.polys[0]))
                            for i,dir in enumerate(dirs):
                                neighbor_coords = [x+dir[0], y+dir[1]]
                                if _helper_in_bounds(neighbor_coords[0], neighbor_coords[1]) and self.polys[neighbor_coords[0]][neighbor_coords[1]].isMixed() and self._get_merge_id(neighbor_coords[0], neighbor_coords[1]) in merge_id_with_neighbors.get_neighbor_ids():
                                    neighbor_dirs.append(i)
                                else:
                                    nonneighbor_dirs.append(i)
                            # TODO go by self._get_neighbor_ids_from_merge_ids instead?
                            # if more than two neighbors, then this cell used to have more than two mixed neighbors but some were eliminated, pass
                            if len(neighbor_dirs) > 2:
                                print("Passing on case where cell used to have more than two mixed neighbors but only two are left, still not easily orientable")
                                process_queue.append(merge_id)
                                iters_without_progress += 1
                                continue
                            elif len(neighbor_dirs) < 2:
                                # TODO this case might never happen
                                print("Passing on case where cell has fewer than two cells, not easily orientable")
                                process_queue.append(merge_id)
                                iters_without_progress += 1
                                continue
                            
                            # Verifies neighbors = coordinate neighbors
                            # print(merge_id_with_neighbors.get_neighbor_ids())
                            # print(list(map(lambda d: self._get_merge_id(x+dirs[d][0], y+dirs[d][1]), neighbor_dirs)))

                            # otherwise, two mixed neighbors are either across from each other or adjacent
                            neighbor_modes = abs(neighbor_dirs[0]-neighbor_dirs[1]) # = 1 or 2 or 3 (3 is same case as 1)
                            # 1,3 = adjacent; 2 = across
                            nonneighbor_statuses = []
                            for nonneighbor_dir in nonneighbor_dirs:
                                nonneighbor_coords = [x+dirs[nonneighbor_dir][0], y+dirs[nonneighbor_dir][1]]
                                if _helper_in_bounds(nonneighbor_coords[0], nonneighbor_coords[1]) and self.polys[nonneighbor_coords[0]][nonneighbor_coords[1]].isFull():
                                    nonneighbor_statuses.append(1)
                                else: # must be empty
                                    nonneighbor_statuses.append(0)
                            if neighbor_modes == 1 or neighbor_modes == 3:
                                clockwisemost = min(neighbor_dirs) if neighbor_modes == 1 else 3
                                [c_x, c_y] = [x+dirs[clockwisemost][0], y+dirs[clockwisemost][1]]
                                c_merge_id = self._get_merge_id(c_x, c_y)
                                [cc_x, cc_y] = [x+dirs[(clockwisemost+1)%4][0], y+dirs[(clockwisemost+1)%4][1]]
                                cc_merge_id = self._get_merge_id(cc_x, cc_y)
                                if nonneighbor_statuses == [0, 0] and not(merge_id_to_obj[c_merge_id].has_left()) and not(merge_id_to_obj[cc_merge_id].has_right()): # both empty, clockwise
                                    merge_id_with_neighbors.set_left(self._get_merge_id(cc_x, cc_y))
                                    merge_id_with_neighbors.set_right(self._get_merge_id(c_x, c_y))
                                    iters_without_progress = 0
                                elif nonneighbor_statuses == [1, 1] and not(merge_id_to_obj[c_merge_id].has_right()) and not(merge_id_to_obj[cc_merge_id].has_left()): # both full, counterclockwise
                                    merge_id_with_neighbors.set_left(self._get_merge_id(c_x, c_y))
                                    merge_id_with_neighbors.set_right(self._get_merge_id(cc_x, cc_y))
                                    iters_without_progress = 0
                                else:
                                    print(f"Error in easy orientation with two adjacent neighbors: {nonneighbor_statuses}")
                                    process_queue.append(merge_id)
                                    iters_without_progress += 1
                            else: # neighbor_modes == 2
                                if nonneighbor_statuses == [0, 1]:
                                    full_index = 1
                                elif nonneighbor_statuses == [1, 0]:
                                    full_index = 0
                                else:
                                    print(f"Error in easy orientation with two opposite neighbors: {nonneighbor_statuses}")
                                    process_queue.append(merge_id)
                                    iters_without_progress += 1
                                    break
                                [l_x, l_y] = [x+dirs[(nonneighbor_dirs[full_index]+1)%4][0], y+dirs[(nonneighbor_dirs[full_index]+1)%4][1]]
                                l_merge_id = self._get_merge_id(l_x, l_y)
                                [r_x, r_y] = [x+dirs[(nonneighbor_dirs[full_index]-1)%4][0], y+dirs[(nonneighbor_dirs[full_index]-1)%4][1]]
                                r_merge_id = self._get_merge_id(r_x, r_y)
                                if not(merge_id_to_obj[l_merge_id].has_right()) and not(merge_id_to_obj[r_merge_id].has_left()):
                                    merge_id_with_neighbors.set_left(l_merge_id)
                                    merge_id_with_neighbors.set_right(r_merge_id)
                                    iters_without_progress = 0
                                else:
                                    print("Error in easy orientation with two opposite neighbors but inconsistent orientations")
                                    process_queue.append(merge_id)
                                    iters_without_progress += 1
                        # Not easily orientable, pass
                        else:
                            print("Passing on case with two neighbors and two missing orientations")
                            process_queue.append(merge_id)
                            iters_without_progress += 1
                    elif merge_id_with_neighbors.fully_oriented():
                        # TODO this doesn't seem to happen
                        print("Fully oriented but two neighbors: why did this happen?")
                        process_queue.append(merge_id)
                        iters_without_progress += 1
                    # One more orientation to be set and two neighbors, pass
                    else:
                        print("Passing on case with two neighbors and one missing orientation")
                        process_queue.append(merge_id)
                        iters_without_progress += 1
                # 3+ neighbors
                elif len(merge_id_with_neighbors.get_neighbor_ids()) >= 3:
                    if merge_id_with_neighbors.fully_oriented():
                        raise ValueError(f"Fully oriented but {len(merge_id_with_neighbors.get_neighbor_ids())} neighbors: why did this happen?")
                    else:
                        print("Passing on case with 3+ neighbors")
                        process_queue.append(merge_id)
                        iters_without_progress += 1
                # 0 neighbors
                else:
                    # if still not fully oriented, these are problematic and we want to do something
                    if not(merge_id_with_neighbors.fully_oriented()):
                        print("Zero neighbors and not full oriented, passing")
                        process_queue.append(merge_id)
                        iters_without_progress += 1

                # iters_without_progress = length + 1 (full cycle of queue without progress)
                if iters_without_progress >= len(process_queue)+1:
                    print(f"Rest of queue cannot be resolved: length {len(process_queue)}")
                    break

        doGreedyOrientations()

        # Cases at this point:
        # 1 neighbor candidate:
            # Two unfilled neighbors
            # One unfilled neighbor but orientations are inconsistent
        # 2 neighbor candidates:
            # Not an "easy orientation" case
        # 3+ neighbor candidates:
            # All
        # 0 neighbor candidates:
            # Not fully oriented

        iters_without_progress = 0
        while process_queue:
            merge_id = process_queue.pop(0)
            merge_id_with_neighbors: MergeIdWithNeighbors = merge_id_to_obj[merge_id]
            print(merge_id)
            if len(merge_id_with_neighbors.get_neighbor_ids()) == 0:
                # No orientations set
                if not(merge_id_with_neighbors.has_left()) and not(merge_id_with_neighbors.has_right()):
                    print("Weird case 1")
                    # Choose a random mixed neighbor and merge with it
                    neighbor_ids = self._get_neighbor_ids_from_merge_id(merge_id)
                    neighbor_id = neighbor_ids[0] # TODO make random, or choose the one in the direction most normal to the interface?
                    new_merge_id_with_neighbors = mergeObjs(merge_id, neighbor_id, use_neighbor_id_orientation=True)
                    process_queue.append(new_merge_id_with_neighbors.get_merge_id())
                    doGreedyOrientations()
                    iters_without_progress = 0
                # Has left neighbor
                elif merge_id_with_neighbors.has_left():
                    # TODO am I ok to set its right neighbor to itself
                    merge_id_with_neighbors.right_neighbor_id = merge_id
                # Has right neighbor
                elif merge_id_with_neighbors.has_right():
                    # TODO am I ok to set its left neighbor to itself
                    merge_id_with_neighbors.left_neighbor_id = merge_id
            elif merge_id_with_neighbors.has_left() and not(merge_id_with_neighbors.has_right()):
                print("Weird case 2")
                left_neighbor_id = merge_id_with_neighbors.get_left()
                neighbor_ids = self._get_neighbor_ids_from_merge_id(merge_id)
                did_update = False
                for i, neighbor_id in enumerate(neighbor_ids):
                    if neighbor_id == left_neighbor_id and not(merge_id_to_obj[neighbor_ids[(i+1) % len(neighbor_ids)]].has_left()):
                        # This uses a hack: since we processed neighbors in counterclockwise order, this should give the mixed neighbor to the right of the left neighbor
                        # TODO is this neighbor guaranteed to not already have a left neighbor?
                        merge_id_with_neighbors.set_right(neighbor_ids[(i+1) % len(neighbor_ids)], set_neighbor=True)
                        did_update = True
                if did_update:
                    doGreedyOrientations()
                    iters_without_progress = 0
                else:
                    process_queue.append(merge_id_with_neighbors.get_merge_id())
                    iters_without_progress += 1
            elif not(merge_id_with_neighbors.has_left()) and merge_id_with_neighbors.has_right():
                print("Weird case 3")
                right_neighbor_id = merge_id_with_neighbors.get_right()
                neighbor_ids = self._get_neighbor_ids_from_merge_id(merge_id)
                did_update = False
                for i, neighbor_id in enumerate(neighbor_ids):
                    if neighbor_id == right_neighbor_id and not(merge_id_to_obj[neighbor_ids[(i-1) % len(neighbor_ids)]].has_right()):
                        # This uses a hack: since we processed neighbors in counterclockwise order, this should give the mixed neighbor to the left of the right neighbor
                        # TODO is this neighbor guaranteed to not already have a right neighbor?
                        merge_id_with_neighbors.set_left(neighbor_ids[(i-1) % len(neighbor_ids)], set_neighbor=True)
                        did_update = True
                if did_update:
                    doGreedyOrientations()
                    iters_without_progress = 0
                else:
                    process_queue.append(merge_id_with_neighbors.get_merge_id())
                    iters_without_progress += 1
            else:
                print("Weird case 4")
                process_queue.append(merge_id_with_neighbors.get_merge_id())
                iters_without_progress += 1
            if iters_without_progress >= len(process_queue)+1:
                print(f"Rest of meta queue cannot be resolved: length {len(process_queue)}")
                break

        # very_ambiguous_ids = process_queue.copy()

        # # Attempts to merge cells that are still unresolved
        # while process_queue:
        #     merge_id_with_neighbors: MergeIdWithNeighbors = process_queue.pop(0)
        #     merge_id = merge_id_with_neighbors.get_merge_id()
        #     # TODO handle other cases (neither left nor right neighbors are set?)
        #     print(f"Unresolved cells: num neighbors = {len(merge_id_with_neighbors.get_neighbor_ids())}")
        #     # Choose a random mixed neighbor and merge with it
        #     neighbor_ids = list(filter(lambda x : x not in list(map(lambda y : y.get_merge_id(), very_ambiguous_ids)), merge_id_with_neighbors.get_neighbor_ids()))
        #     if len(neighbor_ids) == 0:
        #         neighbor_ids = list(filter(lambda x : x not in list(map(lambda y : y.get_merge_id(), very_ambiguous_ids)), self._get_neighbor_ids_from_merge_id(merge_id)))
        #         if len(neighbor_ids) == 0:
        #             #TODO no clue if this is possible
        #             raise ValueError("Is it possible that an entirely unoriented cell also has no neighbors not in process queue here?")
        #     neighbor_id = neighbor_ids[0]
        #     # TODO make random, or choose the one in the direction most normal to the interface?
        #     mergeObjs(merge_id, neighbor_id, use_neighbor_id_orientation=True)

        # Create polygon objects

        # Find the vertices of the merged polygons
        self.createMergedPolys()

        # Merge ids that failed to be oriented
        failed_merge_ids = []
        added_merge_ids = dict()
        for merge_id in merge_id_to_obj.keys():
            merge_id_with_neighbors: MergeIdWithNeighbors = merge_id_to_obj[merge_id]
            if merge_id_with_neighbors.fully_oriented():
                merged_poly: NeighboredPolygon = self.merged_polys[merge_id]
                merged_poly.setNeighbor(self.merged_polys[merge_id_with_neighbors.get_left()], "left")
                merged_poly.setNeighbor(self.merged_polys[merge_id_with_neighbors.get_right()], "right")
                # Check if we used hack and had set left/right neighbor to itself to signify dead-end cell with single mixed neighbor
                if merge_id_with_neighbors.get_left() == merge_id or merge_id_with_neighbors.get_right() == merge_id:
                    merged_poly.setFacetType("linear")
            else:
                print("Final failed orientations:")
                print(merge_id_with_neighbors)
                # For all polys in failed merge id, add a lone polygon with a 3x3 stencil and run Young's
                # May still need this poly because it could be a neighbor of something else so can't directly remove it
                for merge_coords in self._get_merge_coords(merge_id_with_neighbors.get_merge_id()):
                    lone_base_poly = self.polys[merge_coords[0]][merge_coords[1]]
                    merged_poly = NeighboredPolygon(lone_base_poly.points)
                    merged_poly.setFraction(lone_base_poly.getFraction())
                    merged_poly.set3x3Stencil(self.get3x3Stencil(merge_coords[0], merge_coords[1]))

                    # TODO does this throw off the algorithm somehow because we're appending to merge_id_to_obj only? (Some invariant where merge_id_to_obj has to have same length as something else?)
                    self.merged_polys[self.next_merge_id] = merged_poly
                    self.coords_to_merge_id[merge_coords[0]][merge_coords[1]] = self.next_merge_id
                    self.merge_ids_to_coords.append(merge_coords)
                    self.merge_id_to_neighbor_ids.append([])
                    added_merge_ids[self.next_merge_id] = merged_poly
                    self.next_merge_id += 1
                failed_merge_ids.append(merge_id)

        # Remove failed merge ids from list of polygons
        for failed_merge_id in failed_merge_ids:
            merge_id_to_obj.pop(failed_merge_id)
        # Add Young's polygons
        for added_merge_id in added_merge_ids.keys():
            merge_id_to_obj[added_merge_id] = added_merge_ids[added_merge_id]
        return list(merge_id_to_obj.keys())

#TODO return polys instead and write different functions for only linear vs. only circular vs. corners?

    def fitFacets(self, merge_ids, setting="circular"):
        if setting == "linear":
            i = 0
            while i < len(merge_ids):
                merge_id = merge_ids[i]
                merged_poly: NeighboredPolygon = self.merged_polys[merge_id]
                if merged_poly.fullyOriented():
                    merged_poly.fitLinearFacet()
                elif merged_poly.has3x3Stencil():
                    merged_poly.runYoungs()
                else:
                    print("Skip")
                    print(self.merged_polys[merge_id])
                    merge_ids.pop(i)
                    i -= 1
                i += 1

        elif setting == "circular":
            i = 0
            while i < len(merge_ids):
                merge_id = merge_ids[i]
                merged_poly: NeighboredPolygon = self.merged_polys[merge_id]
                if merged_poly.fullyOriented():
                    if merged_poly.facet_type == "linear":
                        merged_poly.fitLinearFacet()
                    else:
                        merged_poly.fitCircularFacet()
                elif merged_poly.has3x3Stencil():
                    merged_poly.runYoungs()
                else:
                    merge_ids.pop(i)
                    i -= 1
                i += 1

        elif setting == "linear+corner":
            pass

        elif setting == "circular+corner":
            pass

        elif setting == "extra_corners":
            pass

        # Return merge_polys
        return list(map(lambda x : self.merged_polys[x], merge_ids))