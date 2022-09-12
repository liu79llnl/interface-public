import argparse
import math

from main.structs.meshes.merge_mesh import MergeMesh
from util.config import read_yaml
from util.init_points import makeFineCartesianGrid
from util.initialize_areas import initializeCircle, initializePoly, zalesak, xpluso
from util.write_facets import writeFacets

import pickle


def main(config_setting):

    config = read_yaml(f"config/{config_setting}.yaml")

    #Test settings
    from_temp = config["TEST"]["FROM_TEMP"]
    save_name = config["TEST"]["SAVE_NAME"]
    test_setting = config["TEST"]["SETTING"]

    #Mesh settings
    grid_size = config["MESH"]["GRID_SIZE"]
    resolution = config["MESH"]["RESOLUTION"]

    #Area and facet settings
    facet_algo = config["GEOMS"]["FACET_ALGO"]
    threshold = config["GEOMS"]["THRESHOLD"]

    #Advection settings
    dt = config["ADVECTION"]["DT"]
    t_total = config["ADVECTION"]["T_TOTAL"]

    #-----

    #Initialize mesh
    if not(from_temp):
        print("Generating mesh...")
        opoints = makeFineCartesianGrid(grid_size, resolution)
        m = MergeMesh(opoints, threshold)
        m.writeMesh(f"plots/{save_name}/vtk/mesh.vtk")

        #Initialize fractions
        if test_setting == "vortex":
            fractions = initializeCircle(m, [50.01, 75.01], 15)
        if test_setting == "deformation":
            fractions = initializeCircle(m, [50.01, 75.01], 15)
        elif test_setting == "zalesak":
            fractions = zalesak(m)
        elif test_setting == "xpluso":
            fractions = xpluso(m)
        elif test_setting == "triangle":
            fractions = initializePoly(m, [[1.5, 1.5], [9.5, 11.5], [5.5, 15.2]])
        elif test_setting == "rectangle":
            fractions = initializePoly(m, [[8.2, 8.2], [12.2, 5.2], [15.2, 9.2], [11.2, 12.2]])
            
        m.initializeFractions(fractions)

        #Plot areas
        m.plotAreas(f"plots/{save_name}/plt/areas/initial.png")
        m.plotPartialAreas(f"plots/{save_name}/plt/partial_areas/initial.png")

        #Merge and run interface reconstruction
        print("Initial interface reconstruction")
        m.merge1Neighbors()
        merge_ids = m.findOrientations()
        merged_polys = m.fitFacets(merge_ids, setting=facet_algo)
        reconstructed_facets = [p.getFacet() for p in merged_polys]
        writeFacets(reconstructed_facets, f"plots/{save_name}/vtk/reconstructed/initial.vtp")
    else:
        print("Loading mesh...")
        temp_f = open(f"plots/{save_name}/temp.pickle", "rb")
        [m, t_init] = pickle.load(temp_f)
        print(f"Skipping to timestep {t_init}...")

    #Initialize velocities
    if test_setting == "vortex":
        velocity = lambda t, p : [-200*math.cos(math.pi*t/t_total)*(math.sin(math.pi*p[0]/100))**2 * math.sin(math.pi*p[1]/100) * math.cos(math.pi*p[1]/100),
                                200*math.cos(math.pi*t/t_total)*math.sin(math.pi*p[0]/100)*math.cos(math.pi*p[0]/100) * (math.sin(math.pi*p[1]/100))**2]
    if test_setting == "deformation":
        velocity = lambda t, p: [100*math.cos(math.pi*t/t_total)*math.sin(4*math.pi*(p[0]+50)/100)*math.sin(4*math.pi*(p[1]+50)/100),
                                100*math.cos(math.pi*t/t_total)*math.cos(4*math.pi*(p[0]+50)/100)*math.cos(4*math.pi*(p[1]+50)/100)]
    elif test_setting == "zalesak":
        velocity = lambda t, p: [(2*math.pi)/t_total*(50-p[1]), (2*math.pi)/t_total*(p[0]-50)]
    elif test_setting == "xpluso":
        velocity = lambda t, p: [6, 6]
    elif test_setting == "triangle":
        velocity = lambda t, p: [6, 6]
    elif test_setting == "rectangle":
        velocity = lambda t, p: [6, 6]

    #Advection
    print("Advection test")
    t = 0
    num_iters = int(t_total/dt)+5
    for iter in range(num_iters):

        # If loading from temp, skip
        if from_temp and iter < t_init:
            t += dt
            continue

        print("t = {}".format(t))

        #Advect facets and compute new areas
        advected_facets = m.advectMergedFacets(velocity, t, dt, checkSize=2)
        writeFacets(advected_facets, f"plots/{save_name}/vtk/advected/{iter}.vtp")

        #Plot areas
        m.plotAreas(f"plots/{save_name}/plt/areas/{iter}.png")
        m.plotPartialAreas(f"plots/{save_name}/plt/partial_areas/{iter}.png")

        #Merge and run interface reconstruction
        m.merge1Neighbors()
        merge_ids = m.findOrientations()
        merged_polys = m.fitFacets(merge_ids, setting=facet_algo)
        reconstructed_facets = [p.getFacet() for p in merged_polys]
        writeFacets(reconstructed_facets, f"plots/{save_name}/vtk/reconstructed/{iter}.vtp")
        
        if from_temp:
            print(1/0)
        m.writeToPickle(f"plots/{save_name}/temp.pickle", iter)

        t += dt

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Advection tests.")
    parser.add_argument("--config", type=str, help="config setting")
    args = parser.parse_args()

    main(config_setting=args.config)