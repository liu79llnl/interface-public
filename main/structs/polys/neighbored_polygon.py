from main.structs.polys.base_polygon import BasePolygon

class NeighboredPolygon(BasePolygon):

    def __init__(self, points):
        super().__init__(points)

        self.left_neighbor = None
        self.right_neighbor = None

        #TODO List of unoriented neighbors, maybe useful later?
        self.unoriented_neighbors = []

    #poly = NeighboredPolygon
    #orientation = "left", "right", other
    def setNeighbor(self, poly, orientation):
        if orientation == "left":
            self.left_neighbor = poly
        elif orientation == "right":
            self.right_neighbor = poly
        else:
            self.unoriented_neighbors.append(poly)

    #orientation = "left", "right", other
    def clearNeighbor(self, orientation):
        if orientation == "left":
            self.left_neighbor = None
        elif orientation == "right":
            self.right_neighbor = None
        else:
            self.unoriented_neighbors = []

    def hasLeftNeighbor(self):
        return (self.left_neighbor is not None)

    def hasRightNeighbor(self):
        return (self.right_neighbor is not None)

    def __str__(self):
        if self.hasFacet:
            return f"\nPoints: {self.points}\nFraction: {self.fraction}\nLeft neighbor: {self.left_neighbor.points if self.left_neighbor is not None else None}\nRight neighbor: {self.right_neighbor.points if self.right_neighbor is not None else None}\nFacet: {self.facet}\n"
        else:
            return f"\nPoints: {self.points}\nFraction: {self.fraction}\nLeft neighbor: {self.left_neighbor.points if self.left_neighbor is not None else None}\nRight neighbor: {self.right_neighbor.points if self.right_neighbor is not None else None}\n"