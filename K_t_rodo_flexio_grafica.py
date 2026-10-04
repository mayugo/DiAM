import numpy as np
from scipy.interpolate import RegularGridInterpolator

def K_t_rodo_flexio_grafica(D_d, r_d):
    r_d = np.maximum(np.minimum(r_d, 0.3), 0.025)
    D_d = np.maximum(np.minimum(D_d, 3), 1.02)

    # r_d data
    r_d_vec = np.array([0.025, 0.050, 0.075, 0.100, 0.125, 0.150,
                        0.175, 0.200, 0.225, 0.250, 0.275, 0.300])

    # D_d data flexió
    D_d_vec_f = np.array([1.02, 1.05, 1.1, 1.5, 3])

    # K_t values flexió (files: D/d ; columnes: r/d)
    data_f = np.array([
        [1.913, 1.667, 1.537, 1.450, 1.388, 1.341, 1.304, 1.268, 1.239, 1.228, 1.209, 1.198],
        [2.140, 1.804, 1.645, 1.544, 1.471, 1.410, 1.359, 1.326, 1.297, 1.275, 1.257, 1.238],
        [2.256, 1.884, 1.699, 1.594, 1.515, 1.464, 1.413, 1.366, 1.366, 1.366, 1.366, 1.366],
        [2.591, 2.075, 1.826, 1.677, 1.572, 1.507, 1.460, 1.416, 1.400, 1.400, 1.400, 1.400],
        [2.852, 2.317, 2.002, 1.803, 1.677, 1.587, 1.510, 1.452, 1.402, 1.362, 1.325, 1.307]])

    # interp2 de MATLAB (interpolació bilineal)
    interp = RegularGridInterpolator((D_d_vec_f, r_d_vec), data_f, method="linear")
    D_d, r_d = np.broadcast_arrays(np.asarray(D_d, dtype=float),
                                   np.asarray(r_d, dtype=float))
    K_t = interp(np.stack([D_d, r_d], axis=-1))
    return K_t.item() if K_t.size == 1 else K_t