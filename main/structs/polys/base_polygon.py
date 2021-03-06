from main.geoms.geoms import getArea

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

    def setFacet(self, facet):
        self.facet = facet

    def clearFacet(self):
        self.facet = None

    def getFacet(self):
        return self.facet

    def hasFacet(self):
        return (self.facet is not None)

    def __str__(self):
        if self.hasFacet:
            return f"\nPoints: {self.points}\nFraction: {self.fraction}\nFacet: {self.facet}\n"
        else:
            return f"\nPoints: {self.points}\nFraction: {self.fraction}\n"

            