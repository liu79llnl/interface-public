from main.structs.polys.base_polygon import BasePolygon

from main.structs.facets.base_facet import LinearFacet, ArcFacet
from main.geoms.circular_facet import getArcFacet
from main.geoms.linear_facet import getLinearFacet, getPolyLineArea, getPolyLineIntersects

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

    def setCircularFacet(self):
        # If both neighbors, try linear and circular TODO
        if self.hasLeftNeighbor() and self.hasRightNeighbor():
            linearerrorthreshold = 1e-6 #if area fraction error in linear facet < linearerrorthreshold, use linear facet at this cell
            threshold = 1e-10 # in optimizations
            facetline1, facetline2 = getLinearFacet(self.left_neighbor.points, self.right_neighbor.points, self.left_neighbor.getFraction(), self.right_neighbor.getFraction(), threshold)
            if abs(self.getFraction() - getPolyLineArea(self.points, facetline1, facetline2)/self.getArea()) < linearerrorthreshold and (getPolyLineArea(self.points, facetline1, facetline2)/self.getArea() > threshold) and (getPolyLineArea(self.points, facetline1, facetline2)/self.getArea() < 1-threshold):
                intersects = getPolyLineIntersects(self.points, facetline1, facetline2)
                self.setFacet(LinearFacet(intersects[0], intersects[-1]))
            else:
                arccenter, arcradius, arcintersects = getArcFacet(self.left_neighbor.points, self.points, self.right_neighbor.points, self.left_neighbor.getFraction(), self.getFraction(), self.right_neighbor.getFraction(), threshold)
                #TODO If failed: default to linear
                if arccenter is None or arcradius is None or arcintersects is None:
                    intersects = getPolyLineIntersects(self.points, facetline1, facetline2)
                    self.setFacet(LinearFacet(intersects[0], intersects[-1]))
                else:
                    self.setFacet(ArcFacet(arccenter, arcradius, arcintersects[0], arcintersects[-1]))
        else:
            print("Failed to make facet")

    #TODO add settings for only circles / only lines / corners

    def __str__(self):
        if self.hasFacet:
            return f"\nPoints: {self.points}\nFraction: {self.fraction}\nLeft neighbor: {self.left_neighbor.points if self.left_neighbor is not None else None}\nRight neighbor: {self.right_neighbor.points if self.right_neighbor is not None else None}\nFacet: {self.facet}\n"
        else:
            return f"\nPoints: {self.points}\nFraction: {self.fraction}\nLeft neighbor: {self.left_neighbor.points if self.left_neighbor is not None else None}\nRight neighbor: {self.right_neighbor.points if self.right_neighbor is not None else None}\n"