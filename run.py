import argparse
import pickle

from main.structs.meshes.merge_mesh import MergeMesh
from util.config import read_yaml
from util.initialize_points import makeFineCartesianGrid
from util.initialize_areas import initializeAreas
from util.initialize_velocity import initializeVelocity
from util.write_facets import writeFacets
from util.plot_plt import plotAreas, plotPartialAreas
from util.plot_vtk import writeMesh
from util.util import writeToPickle

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
        writeMesh(m, f"plots/{save_name}/vtk/mesh.vtk")

        #Initialize fractions
        fractions = initializeAreas(m, test_setting=test_setting)
        m.initializeFractions(fractions)

        #Plot areas
        plotAreas(m, f"plots/{save_name}/plt/areas/initial.png")
        plotPartialAreas(m, f"plots/{save_name}/plt/partial_areas/initial.png")

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

    # Initialize velocities
    velocity = initializeVelocity(m, t_total, test_setting=test_setting)

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
        plotAreas(m, f"plots/{save_name}/plt/areas/{iter}.png")
        plotPartialAreas(m, f"plots/{save_name}/plt/partial_areas/{iter}.png")

        #Merge and run interface reconstruction
        m.merge1Neighbors()
        merge_ids = m.findOrientations()
        m.updatePlts()
        merged_polys = m.fitFacets(merge_ids, setting=facet_algo)
        reconstructed_facets = [p.getFacet() for p in merged_polys]
        writeFacets(reconstructed_facets, f"plots/{save_name}/vtk/reconstructed/{iter}.vtp")
        
        if from_temp:
            print(1/0)
        writeToPickle(m, f"plots/{save_name}/temp.pickle", iter)

        t += dt
        if iter == 10:
            print(1/0)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Advection tests.")
    parser.add_argument("--config", type=str, help="config setting")
    args = parser.parse_args()

    main(config_setting=args.config)
