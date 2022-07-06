import math

from main.structs.meshes.merge_mesh import MergeMesh
from util.init_points import makeFineCartesianGrid
from util.initialize_areas import initializeCircle, initializePoly, zalesak, xpluso
from util.write_facets import writeFacets

import pickle

#Settings

#Mesh settings
grid_size = 100
resolution = 32/100

#Test settings
test_type = "vortex"

#Area and facet settings
threshold = 1e-10

#Advection settings
dt = 0.01
total_t = 2

#-------------------------------------------------
#Initialize mesh
opoints = makeFineCartesianGrid(grid_size, resolution)
m = MergeMesh(opoints, threshold)
m.writeMesh('plots/vtk/mesh.vtk')

#Initialize test
if test_type == 'vortex':
    fractions = initializeCircle(m, [50.01, 75.01], 15)
    velocity = lambda t, p : [-200*math.cos(math.pi*t/total_t)*(math.sin(math.pi*p[0]/100))**2 * math.sin(math.pi*p[1]/100) * math.cos(math.pi*p[1]/100),
                              200*math.cos(math.pi*t/total_t)*math.sin(math.pi*p[0]/100)*math.cos(math.pi*p[0]/100) * (math.sin(math.pi*p[1]/100))**2]
if test_type == 'deformation':
    fractions = initializeCircle(m, [50.01, 75.01], 15)
    velocity = lambda t, p: [100*math.cos(math.pi*t/total_t)*math.sin(4*math.pi*(p[0]+50)/100)*math.sin(4*math.pi*(p[1]+50)/100),
                             100*math.cos(math.pi*t/total_t)*math.cos(4*math.pi*(p[0]+50)/100)*math.cos(4*math.pi*(p[1]+50)/100)]
elif test_type == 'zalesak':
    fractions = zalesak(m)
    velocity = lambda t, p: [(2*math.pi)/total_t*(50-p[1]), (2*math.pi)/total_t*(p[0]-50)]
elif test_type == 'xpluso':
    fractions = xpluso(m)
    velocity = lambda t, p: [6, 6]
elif test_type == 'triangle':
    fractions = initializePoly(m, [[1.5, 1.5], [9.5, 11.5], [5.5, 15.2]])
    velocity = lambda t, p: [6, 6]
elif test_type == 'rectangle':
    fractions = initializePoly(m, [[8.2, 8.2], [12.2, 5.2], [15.2, 9.2], [11.2, 12.2]])
    velocity = lambda t, p: [6, 6]

m.initializeFractions(fractions)

#Plot areas
m.plotAreas('plots/plt/areas/initial.png')
m.plotPartialAreas('plots/plt/partial_areas/initial.png')

#Merge and run interface reconstruction
m.merge1Neighbors()
merged_polys = m.findOrientations()
reconstructed_facets = [p.getFacet() for p in merged_polys]
writeFacets(reconstructed_facets, 'plots/vtk/reconstructed/initial.vtp')

#Advection
t = 0
num_iters = int(total_t/dt)+5
# pickle_f = open('plots/temp.pickle', 'rb')
# m = pickle.load(pickle_f)
for iter in range(num_iters):
    # if iter < 141:
    #     continue
    print("t = {}".format(t))

    #Advect facets and compute new areas
    advected_facets = m.advectMergedFacets(velocity, t, dt, checkSize=2)
    writeFacets(advected_facets, 'plots/vtk/advected/{}.vtp'.format(iter))

    #Plot areas
    m.plotAreas('plots/plt/areas/{}.png'.format(iter))
    m.plotPartialAreas('plots/plt/partial_areas/{}.png'.format(iter))

    #Merge and run interface reconstruction
    m.merge1Neighbors()
    merged_polys = m.findOrientations()
    reconstructed_facets = [p.getFacet() for p in merged_polys]
    writeFacets(reconstructed_facets, 'plots/vtk/reconstructed/{}.vtp'.format(iter))
    
    m.writeToPickle('plots/temp.pickle')

    t += dt
