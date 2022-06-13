from main.structs.meshes.merge_mesh import MergeMesh
from util.init_points import makeFineCartesianGrid
from util.initialize_areas import initializeCircle

gridSize = 100
resolution = 1
opoints = makeFineCartesianGrid(gridSize, resolution)
m = MergeMesh(opoints)

threshold = 1e-8
fractions = initializeCircle(m, [50.01, 75.01], 15, threshold)
m.initializeFractions(fractions)
m.plotPartialAreas('plots/plt/1.png')
m.merge1Neighbors()
merged_polys = m.findOrientations()
for merged_poly in merged_polys:
    print(merged_poly)