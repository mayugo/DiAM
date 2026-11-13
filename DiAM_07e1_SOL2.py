import numpy as np
from scipy.optimize import brentq

# ---------------------------------------------------------------
# ESPECIFICACIONS (mm, N, MPa, kg, rad/s)
# ---------------------------------------------------------------
E = 210e3                       # modul d'elasticitat [MPa]
Sy = 820                        # limit elastic [MPa]
Sut = 950                       # resistencia a traccio [MPa]
zeta = 0.02                     # esmorteiment relatiu
X_lim = 5e-3                    # amplitud desitjada de x [m]
g = 9.81
SF_min = 1.8                    # factor de seguretat minim (estatic i fatiga)
tol_k = 0.05                    # tolerancia relativa de la rigidesa
L_min, L_max, dL = 100, 150, 5  # longitud de flexio [mm]

# configuracions: omega [rad/s], y0 [m], m [kg]
configs = [(60, 1.5e-3, 1.6), (120, 0.9e-3, 1.2),
           (45, 1.4e-3, 1.8), (30, 1.1e-3, 2.0)]

# cataleg de seccions (gruix t, amplada b) [mm]
cataleg = np.array([
    (0.25, 15), (0.25, 20), (0.25, 25), (0.3, 20), (0.3, 30), (0.3, 40), (0.3, 50), (0.4, 20),
    (0.5, 6), (0.5, 8), (0.5, 10), (0.5, 12), (0.5, 15), (0.5, 20), (0.5, 25), (0.5, 28),
    (0.5, 30), (0.5, 40), (0.5, 50), (0.5, 60), (0.5, 70), (0.5, 80), (0.5, 100), (0.6, 40),
    (0.6, 50), (0.8, 10), (0.8, 20), (0.8, 25), (0.8, 30), (0.8, 40), (0.8, 50), (1, 10),
    (1, 12), (1, 15), (1, 20), (1, 25), (1, 30), (1, 35), (1, 40), (1, 50),
    (1, 60), (1, 70), (1, 80), (1, 100), (1.2, 20), (1.2, 30), (1.5, 10), (1.5, 12),
    (1.5, 15), (1.5, 20), (1.5, 25), (1.5, 30), (1.5, 35), (1.5, 40), (1.5, 45), (1.5, 50),
    (1.5, 60), (1.5, 70), (1.5, 80), (1.5, 100), (2, 10), (2, 12), (2, 15), (2, 20),
    (2, 25), (2, 30), (2, 35), (2, 40), (2, 45), (2, 50), (2, 60), (2, 70),
    (2, 80), (2, 100), (2.5, 20), (2.5, 25), (2.5, 30), (2.5, 50), (3, 15), (3, 20),
    (3, 30), (3, 40), (3, 50), (3, 60), (3, 70), (3, 80), (3, 100),
])


# ---------------------------------------------------------------
# PART 1: rigidesa necessaria i forces (resum del SOL1)
# ---------------------------------------------------------------
def transmitancia(Omega, zeta):
    num = np.sqrt(1 + (2*zeta*Omega)**2)
    den = np.sqrt((1 - Omega**2)**2 + (2*zeta*Omega)**2)
    return num/den


def part1(omega, y0, m):
    """Retorna k_req [N/mm], F_m i F_a [N] per suport."""
    Omega = brentq(lambda O: transmitancia(O, zeta) - X_lim/y0, 1.0, 100.0)
    k_req = m*(omega/Omega)**2/1e3
    return k_req, m*g/2, m*omega**2*X_lim/2


# ---------------------------------------------------------------
# PART 2: formulacio del disseny (vectorial: b, t i L poden ser arrays)
# ---------------------------------------------------------------
def rigidesa(b, t, L):
    """Rigidesa total dels dos suports [N/mm] (biga en voladiu)."""
    I = b*t**3/12
    return 2*3*E*I/L**3


def tensio(F, b, t, L):
    """Tensio maxima de flexio a l'encastament [MPa]."""
    return 6*F*L/(b*t**2)


def SF_estatic(F_m, F_a, b, t, L):
    return Sy/tensio(F_m + F_a, b, t, L)


def limit_fatiga(b, t):
    """Limit d'enduriment corregit Se [MPa]."""
    ka = 1.58*Sut**(-0.085)                     # superficie rectificada
    de = 0.808*np.sqrt(b*t)                     # diametre equivalent
    kb = np.where(de < 2.79, 1.0,
                  np.where(de <= 51, (de/7.62)**(-0.107), 1.51*de**(-0.157)))
    kc = 1.0                                    # flexio
    return 0.5*Sut*ka*kb*kc


def SF_fatiga(F_m, F_a, b, t, L):
    """Goodman modificada: sigma_m > 0 (el pes dona F_m > 0)."""
    sig_a = tensio(F_a, b, t, L)
    sig_m = tensio(F_m, b, t, L)
    Se = limit_fatiga(b, t)
    return 1/(sig_a/Se + sig_m/Sut)


def dissenys_factibles(k_req, F_m, F_a):
    """Prova totes les seccions del cataleg i les L cada dL mm.
    Retorna els dissenys factibles ordenats de mes a menys SF."""
    L_llista = np.arange(L_min, L_max + dL, dL)
    t = cataleg[:, 0][:, None]
    b = cataleg[:, 1][:, None]
    L = L_llista[None, :]
    k = rigidesa(b, t, L)
    err = (k - k_req)/k_req
    SFe = SF_estatic(F_m, F_a, b, t, L)
    SFf = SF_fatiga(F_m, F_a, b, t, L)
    SF = np.minimum(SFe, SFf)
    ok = (np.abs(err) <= tol_k) & (SF >= SF_min)
    i, j = np.nonzero(ok)
    dis = np.column_stack([t[i, 0], b[i, 0], L_llista[j], k[i, j],
                           100*err[i, j], SFe[i, j], SFf[i, j], SF[i, j]])
    return dis[np.argsort(-dis[:, 7])]


# ---------------------------------------------------------------
# CALCULS
# ---------------------------------------------------------------
for cfg, (omega, y0, m) in enumerate(configs, start=1):
    k_req, F_m, F_a = part1(omega, y0, m)
    dis = dissenys_factibles(k_req, F_m, F_a)
    print(f'\nConfiguracio {cfg}: k_req = {k_req:.2f} N/mm, '
          f'F_m = {F_m:.2f} N, F_a = {F_a:.2f} N -> {len(dis)} dissenys')
    print('   t x b [mm]   L [mm]  k [N/mm]  err k [%]  SF_est  SF_fat  SF')
    for d in dis:
        print(f'   {d[0]:g} x {d[1]:.0f}'.ljust(14) +
              f'{d[2]:7.0f}  {d[3]:8.2f}  {d[4]:9.2f}  {d[5]:6.2f}  {d[6]:6.2f}  {d[7]:5.2f}')
