from pydoc import doc
from main.geoms.corner_facet import getPolyCornerArea
from main.structs.polys.base_polygon import BasePolygon
from main.geoms.geoms import getArea, getDistance, lineIntersect, pointInPoly
from main.structs.facets.base_facet import LinearFacet, ArcFacet, CornerFacet
from main.geoms.circular_facet import getArcFacet
from main.geoms.linear_facet import getLinearFacet, getPolyLineArea, getPolyLineIntersects, getLinearFacetFromNormal

class NeighboredPolygon(BasePolygon):

    linearity_threshold = 1e-6 #if area fraction error in linear facet < this value, use linear facet at this cell
    optimization_threshold = 1e-10 # in optimizations
    linear_corner_match_threshold = 10 # if area fraction error in corner facet < this value, use (straight edged) corner facet at this cell

    def __init__(self, points):
        super().__init__(points)

        self.facet_type = None
        self.left_neighbor = None
        self.right_neighbor = None

        #TODO List of unoriented neighbors, maybe useful later?
        self.unoriented_neighbors = []

    def setFacetType(self, facet_type):
        self.facet_type = facet_type

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

    def getLeftNeighbor(self):
        return self.left_neighbor

    def getRightNeighbor(self):
        return self.right_neighbor

    def fullyOriented(self):
        return self.hasLeftNeighbor() and self.hasRightNeighbor()

    def fitCircularFacet(self):
        # If both neighbors, try linear and circular TODO
        if self.hasLeftNeighbor() and self.hasRightNeighbor():
            facetline1, facetline2 = getLinearFacet(self.left_neighbor.points, self.right_neighbor.points, self.left_neighbor.getFraction(), self.right_neighbor.getFraction(), NeighboredPolygon.optimization_threshold)
            if abs(self.getFraction() - getPolyLineArea(self.points, facetline1, facetline2)/self.getArea()) < NeighboredPolygon.linearity_threshold and (getPolyLineArea(self.points, facetline1, facetline2)/self.getArea() > NeighboredPolygon.optimization_threshold) and (getPolyLineArea(self.points, facetline1, facetline2)/self.getArea() < 1-NeighboredPolygon.optimization_threshold):
                intersects = getPolyLineIntersects(self.points, facetline1, facetline2)
                self.setFacet(LinearFacet(intersects[0], intersects[-1]))
            else:
                try:
                    arccenter, arcradius, arcintersects = getArcFacet(self.left_neighbor.points, self.points, self.right_neighbor.points, self.left_neighbor.getFraction(), self.getFraction(), self.right_neighbor.getFraction(), NeighboredPolygon.optimization_threshold)
                    #TODO If failed: default to linear
                    if arccenter is None or arcradius is None or arcintersects is None:
                        l1 = facetline1
                        l2 = facetline2
                        normal = [(-l2[1]+l1[1])/getDistance(l1, l2), (l2[0]-l1[0])/getDistance(l1, l2)]
                        facetline1, facetline2 = getLinearFacetFromNormal(self.points, self.getFraction(), normal, NeighboredPolygon.optimization_threshold)
                        self.setFacet(LinearFacet(facetline1, facetline2))
                    else:
                        self.setFacet(ArcFacet(arccenter, arcradius, arcintersects[0], arcintersects[-1]))
                except:
                    l1 = facetline1
                    l2 = facetline2
                    normal = [(-l2[1]+l1[1])/getDistance(l1, l2), (l2[0]-l1[0])/getDistance(l1, l2)]
                    facetline1, facetline2 = getLinearFacetFromNormal(self.points, self.getFraction(), normal, NeighboredPolygon.optimization_threshold)
                    self.setFacet(LinearFacet(facetline1, facetline2))
        else:
            print("Not enough neighbors: failed to make circular facet")

    # if doCollinearityCheck, then only set linear facet if middle area fraction matches within threshold
    def fitLinearFacet(self, doCollinearityCheck=False):
        if self.hasLeftNeighbor() and self.hasRightNeighbor():
            l1, l2 = getLinearFacet(self.left_neighbor.points, self.right_neighbor.points, self.left_neighbor.getFraction(), self.right_neighbor.getFraction(), NeighboredPolygon.optimization_threshold)
            if doCollinearityCheck:
                # Check if this facet matches middle area
                if abs(self.getFraction() - (getPolyLineArea(self.points, l1, l2)/self.getMaxArea())) < NeighboredPolygon.linearity_threshold:
                    intersects = getPolyLineIntersects(self.points, l1, l2)
                    self.setFacet(LinearFacet(intersects[0], intersects[-1]))
            else:
                # Use this as normal and fit a linear facet
                normal = [(-l2[1]+l1[1])/getDistance(l1, l2), (l2[0]-l1[0])/getDistance(l1, l2)]
                #TODO getLinearFacetFromNormal failed on areafraction=1e-8, threshold=1e-12
                try:
                    facetline1, facetline2 = getLinearFacetFromNormal(self.points, self.getFraction(), normal, NeighboredPolygon.optimization_threshold)
                except:
                    print(f"getLinearFacetFromNormal({self.points}, {self.getFraction()}, {normal}, {NeighboredPolygon.optimization_threshold})")
                    facetline1, facetline2 = getLinearFacetFromNormal(self.points, self.getFraction()*1e3, normal, NeighboredPolygon.optimization_threshold)
                
                self.setFacet(LinearFacet(facetline1, facetline2))

        else:
            print("Not enough neighbors: failed to make linear facet")

    # Extends lines l1 to l2, p1 to p2, calculates intersection if it exists, checks whether this corner matches area fraction, updates self.facet if it does
    # (l1, l2, p1, p2) = (l.l, l.r, r.r, r.l)
    def checkCornerFacet(self, l1, l2, p1, p2):
        corner, _, _ = lineIntersect(l1, l2, p1, p2)
        if corner is not None:
            # Feasible intersection, continue #TODO: do a test to see if corner point lies within the two segments
            corner = [l1, corner, p1]
            cornerareafraction = getPolyCornerArea(self.points, corner[0], corner[1], corner[2])/self.getMaxArea()
            if abs(cornerareafraction - self.getFraction()) < NeighboredPolygon.linear_corner_match_threshold:
                self.setFacet(CornerFacet(None, None, None, None, corner[0], corner[1], corner[2]))
            else:
                print("Unmatching corner area:")
                print(cornerareafraction)
                print(self.getFraction())

    #TODO add settings for only circles / only lines / corners

    def __str__(self):
        if self.hasFacet:
            return f"\nPoints: {self.points}\nFraction: {self.fraction}\nLeft neighbor: {self.left_neighbor.points if self.left_neighbor is not None else None}\nRight neighbor: {self.right_neighbor.points if self.right_neighbor is not None else None}\nFacet: {self.facet}\n"
        else:
            return f"\nPoints: {self.points}\nFraction: {self.fraction}\nLeft neighbor: {self.left_neighbor.points if self.left_neighbor is not None else None}\nRight neighbor: {self.right_neighbor.points if self.right_neighbor is not None else None}\n"
            