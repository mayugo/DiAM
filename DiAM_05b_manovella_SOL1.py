"""
Analisi estatica de la manovella (apartat 1)
Factor de seguretat estatic (von Mises) a la seccio 1 (punts a i b) i a la
seccio 2 (punts a i b), per a una distancia d donada de la forca.
J.A. Mayugo (2026)
"""
import numpy as np


# ----------------------------------------------------------------------------
# Funcions
# ----------------------------------------------------------------------------
def section_properties(b, h, t, eta_a, eta_b, alpha=1.1):
    """Propietats de la seccio en U.
    Inputs : amplada b, cantell h, gruix t [mm] i coeficients de torsio
             (alpha, eta_a, eta_b)
    Outputs: inercia a flexio Ix [mm^4], constant de torsio It [mm^4] i
             centre de tallant c [mm]
    """
    Ix = b*h**3/12 - (b - t)*(h - 2*t)**3/12
    It = alpha*(2*eta_a*b*t**3 + eta_b*(h - 2*t)*t**3)
    c = (b - t/2)**2*(h - t)**2*t/4/Ix
    return Ix, It, c


def loads(F, d, s, c, t):
    """Sollicitacions de la seccio.
    Inputs : forca F [N], distancia d de la forca al perfil [mm],
             cota de la seccio s [mm], centre de tallant c i gruix t [mm]
    Outputs: tallant V [N], moment flector M [N*mm] i moment torsor T [N*mm]
    """
    V = F
    M = F*s
    T = F*(d + t/2 - c)
    return V, M, T


def nominal_stresses(V, M, T, b, h, t, Ix, It, point):
    """Tensions nominals en el punt 'a' (fibra extrema) o 'b' (eix neutre).
    Outputs: tensio normal per flexio sig i tallant total tau [MPa]
    """
    tau_T = T*t/It                        # tallant per torsio (igual a i b)
    if point == 'a':
        sig = M*(h/2)/Ix                  # flexio
        Me = (b - t)*t*(h - t)/2          # moment estatic de l'area
    else:
        sig = 0.0                         # a l'eix neutre no hi ha flexio
        Me = b*t*(h - t)/2 + (h/2 - t)*t*(h/2 - t)/2
    tau_V = V*Me/t/Ix                     # tallant per forca tallant
    return sig, tau_V + tau_T


def von_mises(sig, tau):
    """Tensio equivalent de von Mises per a un estat pla (sig, tau)."""
    return np.sqrt(sig**2 + 3*tau**2)


def static_safety_factor(sig, tau, Sy):
    """Factor de seguretat estatic (maxima energia de distorsio)."""
    return Sy/von_mises(sig, tau)


# ----------------------------------------------------------------------------
# Dades del problema
# ----------------------------------------------------------------------------
F = 150.0                 # N, forca aplicada
d = 20.0                  # mm, distancia entre la forca i el perfil en U
Sy = 70.0                 # MPa, limit elastic de l'alumini

t = 4.0                   # mm, gruix (seccions 1 i 2)
sections = {              # s: cota de la seccio; b: amplada; eta: coef. torsio
    'Seccio 1': {'s': 180.0, 'b': 25.0, 'eta': (0.300, 0.306)},
    'Seccio 2': {'s': 15.0,  'b': 20.0, 'eta': (0.291, 0.306)},
}

# ----------------------------------------------------------------------------
# Calcul i resultats
# ----------------------------------------------------------------------------
for name, sec in sections.items():
    h = 30 + (40 - 30)/200*sec['s']                    # mm, cantell
    Ix, It, c = section_properties(sec['b'], h, t, *sec['eta'])
    for point in ('a', 'b'):
        V, M, T = loads(F, d, sec['s'], c, t)
        sig, tau = nominal_stresses(V, M, T, sec['b'], h, t, Ix, It, point)
        SFe = static_safety_factor(sig, tau, Sy)
        print(f"{name}, punt {point}: "
              f"sigma = {sig:5.1f} MPa, tau = {tau:5.1f} MPa, "
              f"sigma_eq = {von_mises(sig, tau):5.1f} MPa, SFe = {SFe:4.2f}")
