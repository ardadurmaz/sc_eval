import numpy as np
import ot
from ot.gromov import gromov_wasserstein

C1 = np.loadtxt('/home/arda/sc_eval/data/DistMat_A.tsv', delimiter='\t', dtype=np.float64)
C2 = np.loadtxt('/home/arda/sc_eval/data/DistMat_B.tsv', delimiter='\t', dtype=np.float64)

p = ot.unif(C1.shape[0])
q = ot.unif(C2.shape[0])

Gg, log = gromov_wasserstein(C1, C2, p, q, loss_fun='square_loss', verbose=False, log=True)
print(log['gw_dist'])
