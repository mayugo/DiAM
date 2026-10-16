# -*- coding: utf-8 -*-
"""
DiAM - Exercici 04: Semieix tubular a flexió amb forat transversal i canvi
de secció.

Càlcul analític de:
  1) les sol·licitacions a les seccions del forat i del canvi de diàmetre;
  2) els factors de concentració de tensions K_t, les tensions nominals i
     màximes, la secció crítica i el coeficient de seguretat respecte a S_ut.

Unitats: N, mm, MPa (N/mm^2), N·mm.
"""

import numpy as np


# =============================================================================
# DADES
# =============================================================================
F = 20_000.0        # [N]   càrrega transversal a l'extrem lliure
L = 210.0           # [mm]  distància de la cara de la placa a la força

D = 50.0            # [mm]  diàmetre exterior del tram gran (on hi ha el forat)
d = 40.0            # [mm]  diàmetre exterior del tram petit
d_i = 20.0          # [mm]  diàmetre del forat axial (interior)
a = 10.0            # [mm]  diàmetre del forat transversal
r = 3.0             # [mm]  radi d'acord del canvi de diàmetre

x_forat = 30.0      # [mm]  posició del forat (des de la cara de la placa)
x_canvi = 60.0      # [mm]  posició del canvi de diàmetre

S_ut = 1800.0       # [MPa] resistència a tracció del 100Cr6 templat

Kt_canvi_grafica = 1.78   # K_t llegit a la gràfica A-15-9 de Shigley (per comparar)


# =============================================================================
# K_t d'un tub amb forat transversal a flexió (taula 3.1)
# =============================================================================
def K_t_forat_flexio(a_D, d_D):
    """Retorna (A, K_t) per a la flexió d'una barra rodona o un tub amb un
    forat transversal, interpolant linealment a la taula 3.1.

    a_D : relació a/D (diàmetre del forat / diàmetre exterior)
    d_D : relació d/D (diàmetre interior / diàmetre exterior); 0 si és massís

    Primer s'interpola cada columna (d/D fix) en funció de a/D i després
    s'interpola el resultat en funció de d/D.
    """
    # --- Dades de la taula 3.1 ---
    t_a_D = np.array([0.050, 0.075, 0.10, 0.125, 0.15, 0.175,
                      0.20, 0.225, 0.25, 0.275, 0.30])
    t_d_D = np.array([0.0, 0.6, 0.9])

    # Files: a/D; columnes: d/D = 0, 0.6, 0.9
    t_A = np.array([
        [0.88, 0.91, 0.92],
        [0.86, 0.88, 0.89],
        [0.83, 0.85, 0.86],
        [0.80, 0.82, 0.82],
        [0.76, 0.79, 0.79],
        [0.72, 0.75, 0.76],
        [0.68, 0.72, 0.73],
        [0.65, 0.68, 0.69],
        [0.61, 0.64, 0.67],
        [0.58, 0.61, 0.66],
        [0.54, 0.58, 0.64],
    ])
    t_Kt = np.array([
        [2.42, 2.55, 2.63],
        [2.35, 2.43, 2.55],
        [2.27, 2.36, 2.49],
        [2.20, 2.32, 2.41],
        [2.15, 2.29, 2.39],
        [2.10, 2.26, 2.38],
        [2.07, 2.23, 2.39],
        [2.04, 2.21, 2.40],
        [2.00, 2.18, 2.42],
        [1.97, 2.16, 2.48],
        [1.94, 2.14, 2.52],
    ])

    # --- Comprovació del rang ---
    if not (t_a_D[0] <= a_D <= t_a_D[-1]):
        raise ValueError(f"a/D = {a_D:.3f} fora del rang de la taula 3.1")
    if not (t_d_D[0] <= d_D <= t_d_D[-1]):
        raise ValueError(f"d/D = {d_D:.3f} fora del rang de la taula 3.1")

    # --- Interpolació lineal doble ---
    def interpola(taula):
        valors_columna = [np.interp(a_D, t_a_D, taula[:, j])
                          for j in range(len(t_d_D))]
        return np.interp(d_D, t_d_D, valors_columna)

    A = interpola(t_A)
    K_t = interpola(t_Kt)
    return A, K_t


# =============================================================================
# K_t d'un rodó amb canvi de secció a flexió (expressió aproximada)
# =============================================================================
def K_t_rodo_flexio(D_d, r_d, i_flag=1):
    # Calcula K_t d'un rodó sotmès a flexió amb canvi de secció D -> d.
    # D = diàmetre gran; d = diàmetre petit; r = radi de l'entalla
    # D_d : relació D/d; r_d : relació r/d
    # i_flag : =1 -> K_t per tensió normal,  !=1 -> K_t per von Mises
    if i_flag == 1:
        K_t = (0.632 + 0.377 * D_d**(-4.4) + r_d**(-0.5)
               * ((-0.14 - 0.363 * D_d**2 + 0.503 * D_d**4)
                  / (1 - 2.39 * D_d**2 + 3.368 * D_d**4))**0.5)
    else:
        K_t = (0.622 + 0.38 * D_d**(-4.3) + r_d**(-0.5)
               * ((-0.322 - 0.277 * D_d**2 + 0.599 * D_d**4)
                  / (1 - 2.55 * D_d**2 + 5.27 * D_d**4))**0.5)
    return K_t


# =============================================================================
# PROPIETATS DE LA SECCIÓ
# =============================================================================
def modul_resistent_tub(D_ext, D_int):
    """Mòdul resistent a flexió d'un tub: Z = pi (D^4 - d^4) / (32 D)."""
    return np.pi * (D_ext**4 - D_int**4) / (32.0 * D_ext)


def area_tub(D_ext, D_int):
    return np.pi / 4.0 * (D_ext**2 - D_int**2)


# =============================================================================
# 1) SOL·LICITACIONS
# =============================================================================
def moment_flector(x):
    """Moment flector d'un voladís amb càrrega F a x = L."""
    return F * (L - x)


V = F
M_forat = moment_flector(x_forat)
M_canvi = moment_flector(x_canvi)
M_encast = moment_flector(0.0)

# Tensió tangencial màxima (fibra neutra) a la secció petita, ordre de magnitud
tau_max = 2.0 * V / area_tub(d, d_i)


# =============================================================================
# 2) FORAT TRANSVERSAL
# =============================================================================
a_D = a / D
dD_forat = d_i / D
A_forat, Kt_forat = K_t_forat_flexio(a_D, dD_forat)

Z_net = A_forat * modul_resistent_tub(D, d_i)   # Z_net = pi A (D^4-d^4)/(32 D)
sigma0_forat = M_forat / Z_net
sigmamax_forat = Kt_forat * sigma0_forat
SF_forat = S_ut / sigmamax_forat


# =============================================================================
# 2) CANVI DE DIÀMETRE
# =============================================================================
D_d = D / d
r_d = r / d
Kt_canvi = K_t_rodo_flexio(D_d, r_d)

Z_canvi = modul_resistent_tub(d, d_i)            # secció de diàmetre petit
sigma0_canvi = M_canvi / Z_canvi
sigmamax_canvi = Kt_canvi * sigma0_canvi
SF_canvi = S_ut / sigmamax_canvi

# Mateix càlcul amb el K_t llegit a la gràfica
sigmamax_canvi_graf = Kt_canvi_grafica * sigma0_canvi
SF_canvi_graf = S_ut / sigmamax_canvi_graf


# =============================================================================
# RESULTATS
# =============================================================================
if __name__ == "__main__":
    print("=" * 66)
    print("1) SOL·LICITACIONS")
    print("=" * 66)
    print(f"Esforç tallant V                    = {V/1000:8.1f} kN")
    print(f"Moment a l'encastament  (x =   0 mm) = {M_encast/1000:8.1f} N·m")
    print(f"Moment al forat         (x = {x_forat:3.0f} mm) = {M_forat/1000:8.1f} N·m")
    print(f"Moment al canvi secció  (x = {x_canvi:3.0f} mm) = {M_canvi/1000:8.1f} N·m")
    print(f"Tensió tangencial màxima (tallant)  = {tau_max:8.1f} MPa")

    print()
    print("=" * 66)
    print("2) CONCENTRACIÓ DE TENSIONS")
    print("=" * 66)
    print(f"Forat:  a/D = {a_D:.3f},  d/D = {dD_forat:.3f}")
    print(f"        A   = {A_forat:.3f},  K_t = {Kt_forat:.3f}")
    print(f"        Z_net = {Z_net:.1f} mm^3")
    print(f"Canvi:  D/d = {D_d:.3f},  r/d = {r_d:.3f}")
    print(f"        K_t = {Kt_canvi:.3f}")
    print(f"        Z   = {Z_canvi:.1f} mm^3")

    print()
    print(f"{'':16s}{'M [N·m]':>10s}{'K_t':>8s}{'sigma_0':>10s}"
          f"{'sigma_max':>12s}{'SF':>8s}")
    print(f"{'':16s}{'':>10s}{'':>8s}{'[MPa]':>10s}{'[MPa]':>12s}{'':>8s}")
    print("-" * 64)
    for nom, M, Kt, s0, smax, SF in [
            ("Forat", M_forat, Kt_forat, sigma0_forat, sigmamax_forat, SF_forat),
            ("Canvi secció", M_canvi, Kt_canvi, sigma0_canvi, sigmamax_canvi, SF_canvi)]:
        print(f"{nom:16s}{M/1000:10.1f}{Kt:8.3f}{s0:10.1f}{smax:12.1f}{SF:8.2f}")

    print()
    if sigmamax_forat >= sigmamax_canvi:
        critica, smax_c, smax_a = "el forat", sigmamax_forat, sigmamax_canvi
    else:
        critica, smax_c, smax_a = "el canvi de secció", sigmamax_canvi, sigmamax_forat
    dif = 100.0 * (smax_c - smax_a) / smax_a
    print(f"Secció crítica: {critica} "
          f"(sigma_max {dif:.1f} % superior a l'altra discontinuïtat)")
    print(f"Coeficient de seguretat mínim: SF = {S_ut / smax_c:.2f}")

    print()
    print("-" * 66)
    print("Comparació amb el K_t llegit a la gràfica del canvi de secció")
    print("-" * 66)
    print(f"K_t funció  = {Kt_canvi:.3f}  ->  sigma_max = {sigmamax_canvi:6.1f} MPa,"
          f"  SF = {SF_canvi:.2f}")
    print(f"K_t gràfica = {Kt_canvi_grafica:.3f}  ->  sigma_max = {sigmamax_canvi_graf:6.1f} MPa,"
          f"  SF = {SF_canvi_graf:.2f}")
    print(f"Diferència en K_t: {100*(Kt_canvi/Kt_canvi_grafica-1):.1f} %")
    critica_graf = "el forat" if sigmamax_forat >= sigmamax_canvi_graf else "el canvi de secció"
    print(f"Amb el K_t de la gràfica, la secció crítica seria {critica_graf}.")