from main.structs.meshes.base_mesh import BaseMesh
from main.structs.polys.base_polygon import BasePolygon
from main.structs.polys.neighbored_polygon import NeighboredPolygon
from main.geoms.geoms import mergePolys, getPolyIntersectArea, getArea, getPolyLineArea
from main.geoms.circular_facet import getCircleIntersectArea
from main.structs.facets.base_facet import advectPoint

from functools import reduce

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
            for neighbor_id in some_neighbor_ids:
                neighbor_ids.add(neighbor_id)
                # Replace old merge id with new one
                self.merge_id_to_neighbor_ids[neighbor_id] = [self.next_merge_id if x == merge_id else x for x in self.merge_id_to_neighbor_ids[neighbor_id]]

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
            # No issue, merge with its single neighbor
            self._merge([merge_id, neighbor_ids[0]])
            # If most recently created merge poly only has one neighbor, add to queue
            if len(self._get_neighbor_ids_from_merge_id(self.next_merge_id-1)) == 1:
                process_queue.append(self.next_merge_id-1)

    # Resets self.merged_polys according to the latest values in self.merge_ids_to_coords
    def createMergedPolys(self):
        #Set global variables to initials
        self.merged_polys = dict()

        # Add all merge ids in use to the queue
        for x in range(len(self.polys)):
            for y in range(len(self.polys[0])):
                merge_id = self._get_merge_id(x, y)
                if merge_id is not None and merge_id not in self.merged_polys.keys():
                    merge_coords = self._get_merge_coords(merge_id)
                    merge_points = reduce(mergePolys, list(map(lambda x : self.polys[x[0]][x[1]].points, merge_coords)), [])
                    self.merged_polys[merge_id] = NeighboredPolygon(merge_points)
                    total_area = sum(list(map(lambda x : self.polys[x[0]][x[1]].getArea(), merge_coords)))
                    self.merged_polys[merge_id].setArea(total_area)

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
                                            assert self.polys[advectx][advecty].getArea() > 1-self.threshold
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

    Returns poly objects
    """
    def findOrientations(self):
        
        merge_id_to_obj = dict()

        class MergeIdWithNeighbors:

            def __init__(self, merge_id, neighbor_ids: list):
                self.merge_id = merge_id
                self.neighbor_ids = neighbor_ids
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
                    for neighbor_id in self.neighbor_ids:
                        self.remove_neighbor_id(neighbor_id)
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
                    for neighbor_id in self.neighbor_ids:
                        self.remove_neighbor_id(neighbor_id)
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
                return f"{self.merge_id}\n{self.neighbor_ids}\n"

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
                    process_queue.append(merge_id_with_neighbors)

        iters_without_progress = 0
        while process_queue:
            merge_id_with_neighbors: MergeIdWithNeighbors = process_queue.pop(0)
            merge_id = merge_id_with_neighbors.get_merge_id()
            if len(merge_id_with_neighbors.get_neighbor_ids()) == 1:
                neighbor_id = merge_id_with_neighbors.get_neighbor_ids()[0]
                # If merge poly only has one unassigned neighbor and only needs to assign one more neighbor, do it
                if merge_id_with_neighbors.has_left() and not(merge_id_with_neighbors.has_right()):
                    merge_id_with_neighbors.set_right(neighbor_id)
                    iters_without_progress = 0
                elif merge_id_with_neighbors.has_right() and not(merge_id_with_neighbors.has_left()):
                    merge_id_with_neighbors.set_left(neighbor_id)
                    iters_without_progress = 0
                # There should never be a case when neither orientation is set but only one neighbor
                elif not(merge_id_with_neighbors.has_left()) and not(merge_id_with_neighbors.has_right()):
                    raise ValueError("Neither orientation is set but one neighbor: why did this happen?")
                # This case should be dealt with via implementation of remove_neighbor_id
                elif merge_id_with_neighbors.fully_oriented():
                    raise ValueError("Fully oriented but one neighbor: why did this happen?")
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
                            if _helper_in_bounds(neighbor_coords[0], neighbor_coords[1]) and self.polys[neighbor_coords[0]][neighbor_coords[1]].isMixed():
                                neighbor_dirs.append(i)
                            else:
                                nonneighbor_dirs.append(i)
                        # if more than two neighbors, then this cell used to have more than two mixed neighbors but some were eliminated, pass
                        if len(neighbor_dirs) > 2:
                            print("Passing on case where cell used to have more than two mixed neighbors but only two are left, still not easily orientable")
                            process_queue.append(merge_id_with_neighbors)
                            iters_without_progress += 1
                            continue
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
                            [cc_x, cc_y] = [x+dirs[(clockwisemost+1)%4][0], y+dirs[(clockwisemost+1)%4][1]]
                            if nonneighbor_statuses == [0, 0]: # both empty, clockwise
                                merge_id_with_neighbors.set_left(self._get_merge_id(cc_x, cc_y))
                                merge_id_with_neighbors.set_right(self._get_merge_id(c_x, c_y))
                                iters_without_progress = 0
                            elif nonneighbor_statuses == [1, 1]: # both full, counterclockwise
                                merge_id_with_neighbors.set_left(self._get_merge_id(c_x, c_y))
                                merge_id_with_neighbors.set_right(self._get_merge_id(cc_x, cc_y))
                                iters_without_progress = 0
                            else:
                                print(neighbor_dirs)
                                for neighbor_dir in neighbor_dirs:
                                    dir = dirs[neighbor_dir]
                                    neighbor_coords = [x+dir[0], y+dir[1]]
                                    print(self.polys[neighbor_coords[0]][neighbor_coords[1]])
                                print(nonneighbor_dirs)
                                print(merge_id_with_neighbors.get_left())
                                print(merge_id_with_neighbors.get_right())
                                raise ValueError(f"Error in easy orientation with two adjacent neighbors: {nonneighbor_statuses}")
                        else: # neighbor_modes == 2
                            if nonneighbor_statuses == [0, 1]:
                                full_index = 1
                            elif nonneighbor_statuses == [1, 0]:
                                full_index = 0
                            else:
                                raise ValueError(f"Error in easy orientation with two opposite neighbors: {nonneighbor_statuses}")
                            [l_x, l_y] = [x+dirs[(nonneighbor_dirs[full_index]+1)%4][0], y+dirs[(nonneighbor_dirs[full_index]+1)%4][1]]
                            [r_x, r_y] = [x+dirs[(nonneighbor_dirs[full_index]-1)%4][0], y+dirs[(nonneighbor_dirs[full_index]-1)%4][1]]
                            merge_id_with_neighbors.set_left(self._get_merge_id(l_x, l_y))
                            merge_id_with_neighbors.set_right(self._get_merge_id(r_x, r_y))
                            iters_without_progress = 0
                    # Not easily orientable, pass
                    else:
                        print("Passing on case with two neighbors and two missing orientations")
                        process_queue.append(merge_id_with_neighbors)
                        iters_without_progress += 1
                elif merge_id_with_neighbors.fully_oriented():
                    raise ValueError("Fully oriented but two neighbors: why did this happen?")
                # One more orientation to be set and two neighbors, pass
                else:
                    print("Passing on case with two neighbors and one missing orientation")
                    process_queue.append(merge_id_with_neighbors)
                    iters_without_progress += 1
            # 3+ neighbors
            elif len(merge_id_with_neighbors.get_neighbor_ids()) >= 3:
                if merge_id_with_neighbors.fully_oriented():
                    raise ValueError(f"Fully oriented but {len(merge_id_with_neighbors.get_neighbor_ids())} neighbors: why did this happen?")
                else:
                    print("Passing on case with 3+ neighbors")
                    process_queue.append(merge_id_with_neighbors)
                    iters_without_progress += 1

            # iters_without_progress = length + 1 (full cycle of queue without progress)
            if iters_without_progress >= len(process_queue)+1:
                print(len(process_queue))
                raise ValueError(f"Failed to merge") #TODO better error message here

        # Create polygon objects

        # Find the vertices of the merged polygons
        self.createMergedPolys()

        for merge_id in merge_id_to_obj.keys():
            merged_poly: NeighboredPolygon = self.merged_polys[merge_id]
            merged_poly.setNeighbor(self.merged_polys[merge_id_to_obj[merge_id].get_left()], "left")
            merged_poly.setNeighbor(self.merged_polys[merge_id_to_obj[merge_id].get_right()], "right")
            merged_poly.setCircularFacet()

        return list(map(lambda x : self.merged_polys[x], merge_id_to_obj.keys()))

#TODO return polys instead and write different functions for only linear vs. only circular vs. corners?
