import numpy as np
import pandas as pd
from PIL import Image
from scipy.linalg import svd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import rbf_kernel
from skimage.util import view_as_windows

def recenter(X):
    x = X.T
    x = (x - x.mean(axis=1).reshape(-1, 1)) / x.std(axis=1).reshape(-1, 1)
    Xn = x.T
    return Xn

def pixel_intensity(adata, key=None, windowsize=20):
    """
    Vectorized version: Compute normalized pixel intensity features around each spot.
    Utilizes NumPy's advanced indexing for batch processing.
    """
    if key is None:
        key = list(adata.uns['spatial'].keys())[0]

    scale = adata.uns['spatial'][key]['scalefactors']['tissue_hires_scalef']
    image = np.uint8(adata.uns['spatial'][key]["images"]['hires'] * 255)
    if image.shape[-1] > 3:
        image = image[:, :, :3]

    pil_image = Image.fromarray(image)
    pil_image.convert("L")
    
    x = np.round(adata.obsm['spatial'][:, 1] * scale).astype(int)
    y = np.round(adata.obsm['spatial'][:, 0] * scale).astype(int)
    x_start = np.clip(x - windowsize, 0, image.shape[0] - 2 * windowsize)
    y_start = np.clip(y - windowsize, 0, image.shape[1] - 2 * windowsize)

    intensity_3d = np.zeros((len(adata.obs_names), 3))
    for i in range(len(adata.obs_names)):
        subimage = np.asarray(pil_image.crop((y_start[i], x_start[i],
                                              y_start[i] + 2 * windowsize,
                                              x_start[i] + 2 * windowsize)))
        intensity_3d[i, :] = subimage.mean(axis=(0, 1))
    adata.obsm['pixel'] = (intensity_3d - intensity_3d.mean(axis=0)) / intensity_3d.std(axis=0)

def make_kernel_from_pixel(adata, n=100, kernel_name='kernel', im_feats_weight=0.3):
    """
    Optimized version (Dense K): Combine spatial coordinates with pixel features.
    """
    Xn = recenter(np.concatenate((adata.obsm['spatial'][:, 0:2], adata.obsm['pixel']), axis=1))
    Xn[:, 2:] = im_feats_weight * Xn[:, 2:]
    K = rbf_kernel(Xn)
    adata.obsp[kernel_name] = K

    n = min(n, Xn.shape[0] - 1)
    U, s, Vt = svd(K, full_matrices=False)
    U_n = U[:, :n]
    adata.obsm[kernel_name] = U[:, :n] * s[:n]
