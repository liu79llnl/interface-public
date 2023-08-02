import argparse
import pickle
import numpy as np

from main.structs.meshes.merge_mesh import MergeMesh
from util.config import read_yaml
from util.initialize_points import makeFineCartesianGrid
from util.initialize_areas import initializeAreas
from util.initialize_velocity import initializeVelocity
# from util.write_facets import writeFacets
from util.plot_plt import plotAreas, plotPartialAreas, plotInitialAreaCompare
from util.plot_vtk import writeMesh, writePartialCells, writeFacets
from util.util import writeToPickle
from util.metrics import trueFinalAreas, L2ErrorFractions

def main(config_setting):

    config = read_yaml(f"config/{config_setting}.yaml")

    # Test settings
    from_temp = config["TEST"]["FROM_TEMP"]
    save_name = config["TEST"]["SAVE_NAME"]
    test_setting = config["TEST"]["SETTING"]

    # Mesh settings
    grid_size = config["MESH"]["GRID_SIZE"]
    resolution = config["MESH"]["RESOLUTION"]

    # Area and facet settings
    facet_algo = config["GEOMS"]["FACET_ALGO"]
    threshold = config["GEOMS"]["THRESHOLD"]
    do_c0 = config["GEOMS"]["DO_C0"]

    # Advection settings
    do_advect = config["ADVECTION"]["DO_ADVECT"]
    if do_advect:
        dt = config["ADVECTION"]["DT"]
        t_total = config["ADVECTION"]["T_TOTAL"]

    #-----

    # Initialize mesh
    if not(from_temp):
        print("Generating mesh...")
        opoints = makeFineCartesianGrid(grid_size, resolution)
        m = MergeMesh(opoints, threshold)
        writeMesh(m, f"plots/{save_name}/vtk/mesh.vtk")

        # Initialize fractions
        fractions = initializeAreas(m, test_setting=test_setting)
        m.initializeFractions(fractions)

        if do_advect:
            # Calculate final areas
            true_final_areas = trueFinalAreas(m, test_setting=test_setting, t=t_total)
    
        # Plot areas
        plotAreas(m, f"plots/{save_name}/plt/areas/initial.png")
        plotPartialAreas(m, f"plots/{save_name}/plt/partial_areas/initial.png")

        # Merge and run interface reconstruction
        print("Initial interface reconstruction")

        if facet_algo in ["Youngs", "LVIRA"]:
            m.createMergedPolys()
            plotAreas(m, f"plots/{save_name}/plt/areas/0.png")
            plotPartialAreas(m, f"plots/{save_name}/plt/partial_areas/0.png")
            writePartialCells(m, f"plots/{save_name}/vtk/reconstructed/mixed_cells/0.vtp")
            if facet_algo == "Youngs":
                m.runYoungs()
            elif facet_algo == "LVIRA":
                m.runLVIRA()
            reconstructed_facets = [p.getFacet() for p in m.merged_polys.values()]
        else:
            m.merge1Neighbors()
            merge_ids = m.findOrientations()
            # Update plot variables with merged cells
            m.updatePlots()
            plotAreas(m, f"plots/{save_name}/plt/areas/0.png")
            plotPartialAreas(m, f"plots/{save_name}/plt/partial_areas/0.png")
            writePartialCells(m, f"plots/{save_name}/vtk/reconstructed/mixed_cells/0.vtp")
            merged_polys = m.fitFacets(merge_ids, setting=facet_algo)
            reconstructed_facets = [p.getFacet() for p in merged_polys]

            # Run C0
            if do_c0:
                merged_polys = m.makeC0(merged_polys)
                C0_facets = [p.getFacet() for p in merged_polys]
                writeFacets(C0_facets, f"plots/{save_name}/vtk/reconstructed/C0_facets/0.vtp")
        
        writeFacets(reconstructed_facets, f"plots/{save_name}/vtk/reconstructed/facets/0.vtp")

    else:
        print("Loading mesh...")
        temp_f = open(f"plots/{save_name}/temp.pickle", "rb")
        [m, t_init] = pickle.load(temp_f)
        print(f"Skipping to timestep {t_init}...")

    if do_advect:

        # Initialize velocities
        velocity = initializeVelocity(m, t_total, test_setting=test_setting)

        # Advection
        print("Advection test")
        t = 0
        num_iters = int(t_total/dt)
        for iter in range(1, num_iters+1):

            # If loading from temp, skip
            if from_temp and iter < t_init:
                t += dt
                continue

            print("t = {}".format(t))

            # Advect facets and compute new areas
            advected_facets = m.advectMergedFacets(velocity, t, dt, checkSize=2)
            writeFacets(advected_facets, f"plots/{save_name}/vtk/advected/facets/{iter}.vtp")

            # Plot areas
            plotAreas(m, f"plots/{save_name}/plt/areas/{iter}.png")
            plotPartialAreas(m, f"plots/{save_name}/plt/partial_areas/{iter}.png")

            # Merge and run interface reconstruction
            if facet_algo in ["Youngs", "LVIRA"]:
                m.createMergedPolys()
                writePartialCells(m, f"plots/{save_name}/vtk/reconstructed/mixed_cells/{iter}.vtp")
                if facet_algo == "Youngs":
                    m.runYoungs()
                elif facet_algo == "LVIRA":
                    m.runLVIRA()
                reconstructed_facets = [p.getFacet() for p in m.merged_polys.values()]
            else:
                m.merge1Neighbors()
                merge_ids = m.findOrientations()
                # Update plot variables with merged cells
                m.updatePlots()
                writePartialCells(m, f"plots/{save_name}/vtk/reconstructed/mixed_cells/{iter}.vtp")
                merged_polys = m.fitFacets(merge_ids, setting=facet_algo)
                reconstructed_facets = [p.getFacet() for p in merged_polys]

                # Run C0
                if do_c0:
                    merged_polys = m.makeC0(merged_polys)
                    C0_facets = [p.getFacet() for p in merged_polys]
                    writeFacets(C0_facets, f"plots/{save_name}/vtk/reconstructed/C0_facets/{iter}.vtp")
                
            writeFacets(reconstructed_facets, f"plots/{save_name}/vtk/reconstructed/facets/{iter}.vtp")
            
            if from_temp:
                print(1/0)
                
            writeToPickle(m, f"plots/{save_name}/temp.pickle", iter)

            t += dt
            # if iter == 11:
            #     print(1/0)

        # Compute final metrics
        # Volume errors
        volume_l2_error, l2_mixed_count = L2ErrorFractions(m.getFractions(), true_final_areas)
        with open(f"plots/{save_name}/volume_l2_error.txt", "w") as f:
            f.write(f"{volume_l2_error}\n{l2_mixed_count}\n")
        plotInitialAreaCompare(m, f"plots/{save_name}/plt/{iter}_initial_compare.png")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Advection tests.")
    parser.add_argument("--config", type=str, help="config setting")
    args = parser.parse_args()

    main(config_setting=args.config)
