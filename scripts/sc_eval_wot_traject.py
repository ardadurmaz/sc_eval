import sys
import wot
import numpy as np

# Read Input
TMAP_PATH = sys.argv[1]
LOCAL_TIME = sys.argv[2]
tmap_model = wot.tmap.TransportMapModel.from_directory(TMAP_PATH)
cell_sets = wot.io.read_sets('/home/arda/sc_eval/data/CellSet.gmt', as_dict=True)
pops = tmap_model.population_from_cell_sets(cell_sets, at_time=LOCAL_TIME)

trajectory_ds = tmap_model.trajectories(pops)
np.savetxt(X=trajectory_ds.X, fname='/home/arda/sc_eval/data/TrajectoryWOT.txt', delimiter='\t')
