import numpy as np
from scipy.spatial.distance import pdist, cdist


def native_masked_hamming_distance_pdist(X, pseudo=0):
    isnan_X = np.isnan(X)
    displaced_cityblock = pdist(np.nan_to_num(X, nan=0.5), metric="cityblock")
    mismatched_nan_count = pdist(isnan_X, metric="cityblock")
    cityblock = displaced_cityblock - (mismatched_nan_count / 2)
    anynan_count = (
        np.nan_to_num(pdist(isnan_X, metric="kulczynski1"), nan=0)
        * mismatched_nan_count
        + mismatched_nan_count
    )
    unmasked_positions = X.shape[1] - anynan_count
    masked_hamming_distance = (cityblock + pseudo) / (unmasked_positions + pseudo)
    return masked_hamming_distance


def native_masked_hamming_distance_cdist(X, Y, pseudo=0):
    assert X.shape[1] == Y.shape[1]
    isnan_X, isnan_Y = np.isnan(X), np.isnan(Y)
    displaced_cityblock = cdist(
        np.nan_to_num(X, nan=0.5), np.nan_to_num(Y, nan=0.5), metric="cityblock"
    )
    mismatched_nan_count = cdist(isnan_X, isnan_Y, metric="cityblock")
    cityblock = displaced_cityblock - (mismatched_nan_count / 2)
    anynan_count = (
        np.nan_to_num(cdist(isnan_X, isnan_Y, metric="kulczynski1"), nan=0)
        * mismatched_nan_count
        + mismatched_nan_count
    )
    unmasked_positions = X.shape[1] - anynan_count
    masked_hamming_distance = (cityblock + pseudo) / (unmasked_positions + pseudo)
    return masked_hamming_distance
