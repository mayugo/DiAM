import numpy as np
from scipy.interpolate import RegularGridInterpolator


def K_ts_rodo_torsor_grafica(D_d, r_d):
    r_d = np.maximum(np.minimum(r_d, 0.3), 0.025)
    D_d = np.maximum(np.minimum(D_d, 2), 1.09)

    # r_d data
    r_d_vec = np.array([0.025, 0.050, 0.075, 0.100, 0.125, 0.150,
                        0.175, 0.200, 0.225, 0.250, 0.275, 0.300])

    # D_d data torsió
    D_d_vec_t = np.array([1.09, 1.20, 1.33, 2])

    # K_ts values torsió (files: D/d ; columnes: r/d)
    data_t = np.array([
        [1.489, 1.290, 1.209, 1.165, 1.140, 1.125, 1.114, 1.106, 1.100, 1.100, 1.100, 1.100],
        [1.915, 1.592, 1.434, 1.345, 1.283, 1.239, 1.195, 1.162, 1.140, 1.125, 1.118, 1.105],
        [2.033, 1.695, 1.518, 1.426, 1.353, 1.305, 1.254, 1.224, 1.195, 1.180, 1.173, 1.170],
        [2.154, 1.768, 1.573, 1.467, 1.397, 1.342, 1.301, 1.261, 1.228, 1.202, 1.191, 1.180]])

    # interp2 de MATLAB (interpolació bilineal)
    interp = RegularGridInterpolator((D_d_vec_t, r_d_vec), data_t, method="linear")
    D_d, r_d = np.broadcast_arrays(np.asarray(D_d, dtype=float),
                                   np.asarray(r_d, dtype=float))
    K_ts = interp(np.stack([D_d, r_d], axis=-1))
    return K_ts.item() if K_ts.size == 1 else K_ts