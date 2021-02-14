import sys
import os
import numpy as np
import leidenalg
import igraph as ig
from sklearn.neighbors import NearestNeighbors

os.chdir(sys.argv[1])
dim_res = np.loadtxt('ReducedDimension.txt')
nbrs = NearestNeighbors(n_neighbors=15, algorithm='ball_tree').fit(dim_res)
local_g = ig.Graph.Adjacency(nbrs.kneighbors_graph(dim_res).toarray().tolist())
part = leidenalg.find_partition(local_g, leidenalg.ModularityVertexPartition)
np.savetxt(fname='LeidenClusters.txt', X=np.asarray(part.membership))
