import sys
import numpy as np
import scanpy as sc
from scanpy.tools._utils import get_init_pos_from_paga as get_paga

wd = sys.argv[1]
adata = sc.read_text(filename="{}/NormExpr.txt".format(wd))
sc.pp.neighbors(adata, use_rep='X', n_neighbors=10, n_pcs=50)
sc.tl.leiden(adata)
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata)
sc.tl.umap(adata, init_pos=get_paga(adata), n_components=2)
np.savetxt(X=adata.obsm['X_umap'], fname='{}/UMAP_Paga.txt'.format(wd), delimiter='\t')
