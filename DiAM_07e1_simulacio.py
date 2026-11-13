import numpy as np
from scipy.integrate import solve_ivp


def simulacio_transitori(m, k, zeta, y0, omega, t_end=None, n_per=60):
    """Simula el moviment vertical x(t) d'una massa m (kg) sobre una molla
    de rigidesa total k (N/m) i esmorteiment relatiu zeta, quan la base
    es mou com y = y0*sin(omega*t) (y0 en m, omega en rad/s).
    Condicions inicials de repos (x = 0, dx/dt = 0).

    t_end : temps de simulacio en s (per defecte, 6 constants de temps)
    n_per : punts de sortida per periode d'excitacio

    Retorna:
      t, x, y : temps, desplacament del motlle i de la base
      F       : forca dinamica sobre cada suport (N), F = (k*(x-y) + c*(dx-dy))/2
      X_ss    : amplitud de x a regim permanent (ultims 5 periodes)
      F_ss    : amplitud de F a regim permanent (ultims 5 periodes)"""
    wn = np.sqrt(k/m)
    c = 2*zeta*np.sqrt(k*m)
    if t_end is None:
        t_end = 6/(zeta*wn)
    T = 2*np.pi/omega
    t = np.arange(0, t_end, T/n_per)

    def f(t, z):
        x, v = z
        y = y0*np.sin(omega*t)
        yd = y0*omega*np.cos(omega*t)
        return [v, (-c*(v - yd) - k*(x - y))/m]

    sol = solve_ivp(f, [0, t[-1]], [0.0, 0.0], t_eval=t, rtol=1e-8, atol=1e-11)
    x = sol.y[0]
    v = sol.y[1]
    y = y0*np.sin(omega*t)
    yd = y0*omega*np.cos(omega*t)
    F = (k*(x - y) + c*(v - yd))/2
    ss = t >= t[-1] - 5*T
    X_ss = np.max(np.abs(x[ss]))
    F_ss = np.max(np.abs(F[ss]))
    return t, x, y, F, X_ss, F_ss
