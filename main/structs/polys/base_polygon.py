from main.geoms.geoms import getArea
from main.structs.facets.base_facet import LinearFacet
from main.algos.plic_normals import getYoungsNormal
from main.geoms.linear_facet import getPolyLineIntersects, getLinearFacetFromNormal

class BasePolygon:

    #Invariant: self.fraction should always be between 0 and 1
    def __init__(self, points):
        self.points = points
        self.max_area = getArea(self.points)
        self.fraction = None
        self.facet = None

        #TODO read from config
        #Error less than this is acceptable
        self.fraction_tolerance = 1e-10

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

    def runYoungs(self):
        assert self.stencil is not None
        threshold = 1e-10 # linear facet optimization
        normal = getYoungsNormal(self.stencil)
        facetline1, facetline2 = getLinearFacetFromNormal(self.points, self.getFraction(), normal, threshold)
        intersects = getPolyLineIntersects(self.points, facetline1, facetline2)
        self.setFacet(LinearFacet(intersects[0], intersects[-1]))

    def __str__(self):
        if self.hasFacet:
            return f"\nPoints: {self.points}\nFraction: {self.fraction}\nFacet: {self.facet}\n"
        else:
            return f"\nPoints: {self.points}\nFraction: {self.fraction}\n"

            