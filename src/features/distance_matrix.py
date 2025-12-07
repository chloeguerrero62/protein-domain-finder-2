import numpy as np
from scipy.spatial.distance import cdist

def compute_distance_matrix(coords):
    if not isinstance(coords, np.ndarray):
        coords = np.array(coords)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"Expected shape (n, 3), got {coords.shape}")
    return cdist(coords, coords, metric='euclidean')
