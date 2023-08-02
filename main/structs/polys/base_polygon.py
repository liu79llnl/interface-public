from main.geoms.geoms import getArea, getDistance, getPolyLineArea, lineIntersect, pointInPoly
from main.structs.facets.base_facet import LinearFacet, ArcFacet
from main.algos.plic_normals import getYoungsNormal, getLVIRANormal
from main.geoms.linear_facet import getPolyLineIntersects, getLinearFacetFromNormal
from main.geoms.circular_facet import matchArcArea, getCenter

class BasePolygon:

    fraction_tolerance = 1e-10
    C0_linear_tolerance = 1e-5 # should be a couple orders of magnitude higher than linearity_threshold in NeighboredPolygon

    #Invariant: self.fraction should always be between 0 and 1
    def __init__(self, points):
        self.points = points
        self.max_area = abs(getArea(self.points))
        self.fraction = None
        self.facet = None

        #TODO read from config
        #Error less than this is acceptable
        self.fraction_tolerance = BasePolygon.fraction_tolerance

        # 3x3 stencil of fractions for Young's
        self.stencil = None

    #TODO set fraction to 0 or 1 if within threshold?
    def setArea(self, area):
        self.fraction = area/self.max_area
        if self.fraction < 0 or self.fraction > 1:
            raise ValueError(f"Fraction {self.fraction} is invalid for BasePolygon after setArea")

    def setFraction(self, fraction):
        self.fraction = fraction
        if self.fraction < 0 or self.fraction > 1:
            raise ValueError(f"Fraction {self.fraction} is invalid for BasePolygon after setFraction")

    def clearFraction(self):
        self.fraction = None

    def setFractionTolerance(self, fraction_tolerance):
        self.fraction_tolerance = fraction_tolerance

    def getArea(self):
        return self.fraction*self.max_area

    def getFraction(self):
        return self.fraction

    def getFractionTolerance(self):
        return self.fraction_tolerance

    def getMaxArea(self):
        return self.max_area

    def isFull(self):
        return self.fraction > 1-self.fraction_tolerance

    def isEmpty(self):
        return self.fraction < self.fraction_tolerance

    def isMixed(self):
        return not(self.isFull() or self.isEmpty())

    def diffFromMixed(self):
        if self.isMixed():
            return 0
        else:
            return max(self.fraction_tolerance-self.fraction, self.fraction-(1-self.fraction_tolerance))

    def setFacet(self, facet):
        self.facet = facet

    def clearFacet(self):
        self.facet = None

    def getFacet(self):
        return self.facet

    def hasFacet(self):
        return (self.facet is not None)

    def set3x3Stencil(self, stencil):
        self.stencil = stencil

    def has3x3Stencil(self):
        return (self.stencil is not None)

    def runYoungs(self, ret=False):
        assert self.has3x3Stencil()
        threshold = 1e-10 # linear facet optimization
        normal = getYoungsNormal(self.stencil)
        facetline1, facetline2 = getLinearFacetFromNormal(self.points, self.getFraction(), normal, threshold)
        intersects = getPolyLineIntersects(self.points, facetline1, facetline2)
        youngsFacet = LinearFacet(intersects[0], intersects[-1], name="Youngs")
        if ret:
            return youngsFacet
        else:
            self.setFacet(youngsFacet)

    def runLVIRA(self, ret=False):
        assert self.has3x3Stencil()
        threshold = 1e-10 # linear facet optimization
        normal = getLVIRANormal(self.stencil)
        facetline1, facetline2 = getLinearFacetFromNormal(self.points, self.getFraction(), normal, threshold)
        intersects = getPolyLineIntersects(self.points, facetline1, facetline2)
        youngsFacet = LinearFacet(intersects[0], intersects[-1], name="LVIRA")
        if ret:
            return youngsFacet
        else:
            self.setFacet(youngsFacet)

    # Given two endpoints of facet and area fraction, find unique curvature satisfying those constraints
    def fitCurvature(self, pLeft, pRight, fraction_tolerance=fraction_tolerance, ret=False):
        d = getDistance(pLeft, pRight)
        lineArea = getPolyLineArea(self.points, pLeft, pRight)
        # If line area is close to target area, return linear facet
        if abs(self.getArea()-lineArea)/self.getMaxArea() < BasePolygon.C0_linear_tolerance:
            facet = LinearFacet(pLeft, pRight)
        # Otherwise, need to fit curvature
        else:
            radius = matchArcArea(d, self.getArea()-lineArea, fraction_tolerance*self.getMaxArea())
            if radius == float("inf") or radius == -float("inf"):
                facet = LinearFacet(pLeft, pRight)
            else:
                center = getCenter(pLeft, pRight, radius)
                facet = ArcFacet(center, radius, pLeft, pRight)
        if ret:
            return facet
        else:
            self.setFacet(facet)

    # Given two endpoints of facet and slopes (both pointing from endpoint toward corner), form unique linear corner and #TODO
    def checkCorner(self, pLeft, pRight, slopeLeft, slopeRight):
        p2Left = [pLeft[0]+slopeLeft[0], pLeft[1]+slopeLeft[1]]
        p2Right = [pRight[0]+slopeRight[0], pRight[1]+slopeRight[1]]
        corner, _, _ = lineIntersect(pLeft, p2Left, pRight, p2Right)
        

    def __str__(self):
        if self.hasFacet:
            return f"\nPoints: {self.points}\nFraction: {self.fraction}\nFacet: {self.facet}\n"
        else:
            return f"\nPoints: {self.points}\nFraction: {self.fraction}\n"
