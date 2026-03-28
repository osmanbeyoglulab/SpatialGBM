import anndata as ann
import json, os
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from glob import glob
from matplotlib.image import imread
from pathlib import Path


def read_visium_sge(sample_id='Paraent_Visium_Human_Glioblastoma', min_cells=5, min_counts=500):
    fpath = f'data/{sample_id}'
    if os.path.exists(f'{fpath}/filtered_feature_bc_matrix.h5'):
        adata = sq.read.visium(f'data/{sample_id}')
    else:
        adata = sc.datasets.visium_sge(sample_id)
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    adata.layers['raw'] = adata.X.copy()
    return adata


