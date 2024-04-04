#!/bin/bash

# Settings
# settings=("vortex" "x+o" "zalesak")
settings=("vortex" "x+o" "zalesak")
# Resolutions
# resolutions=(50 100)
resolutions=(150)
# resolutions=(32 64 128)
# Algos
# algos=("youngs" "lvira" "linear" "circular" "lcorner" "ccorner")
# algos=("youngs" "lvira" "linear" "circular")
algos=("lcorner" "ccorner")

for setting in "${settings[@]}"; do
    for resolution in "${resolutions[@]}"; do
        for algo in "${algos[@]}"; do
            command="python3 run.py --config advection/${setting}/${resolution}/${setting}_${resolution}_${algo}"
            echo "Running: $command"
            eval "$command"
        
        done
    done
done
