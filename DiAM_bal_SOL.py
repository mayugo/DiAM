# DiAM_bal_SOL.py
# Traducci\u00f3 del script MATLAB a Python
# DiAM, 2024 -> Python port

import numpy as np
import matplotlib.pyplot as plt


def f_coef_torsio_rec(h_b):
    """Retorna els coeficients eta1, eta2, eta3 segons la relació h/b.
    Es limita h/b a 100 com al script original."""
    h_b_ = np.array([1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 100])
    eta_1 = np.array([1, 0.858, 0.796, 0.775, 0.753, 0.743, 0.743, 0.743, 0.743, 0.743, 0.743, 0.743, 0.743])
    eta_2 = np.array([0.208, 0.238, 0.256, 0.269, 0.278, 0.290, 0.298, 0.303, 0.307, 0.310, 0.312, 0.314, 0.333])
    eta_3 = np.array([0.141, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.299, 0.303, 0.307, 0.310, 0.313, 0.333])

    hb = min(h_b, 100)
    eta1 = np.interp(hb, h_b_, eta_1)
    eta2 = np.interp(hb, h_b_, eta_2)
    eta3 = np.interp(hb, h_b_, eta_3)

    return eta1, eta2, eta3


def f_invariants2D(sigX, sigY, tauXY):
    """Calcula els principals i altres invariants per un estat pla de tensions.
    Retorna: sig_A (max), sig_B (min), tau_max, sig_VM (Von Mises).
    """
    S = np.array([[sigX, tauXY], [tauXY, sigY]])
    eigs = np.linalg.eigvals(S)
    sigs = np.sort(eigs)[::-1]  # descending
    sig_A = np.real(sigs[0])
    sig_B = np.real(sigs[1])
    tau_max = 0.5 * (sig_A - sig_B)
    sig_VM = np.sqrt(sig_A**2 + sig_B**2 - sig_A*sig_B)
    return sig_A, sig_B, tau_max, sig_VM


def f_cercle_Mohr2D(sigma_A, sigma_B, desc_punt):
    """Dibuixa el cercle de Mohr per un estat pla donat sigma_A i sigma_B.
    Les tensions s'espera que estiguin en N/mm^2 (MPa).
    """
    theta = np.linspace(0, 2*np.pi, 400)
    center = 0.5 * (sigma_A + sigma_B)
    radius = 0.5 * (sigma_A - sigma_B)

    fig, ax = plt.subplots(figsize=(6,6))
    ax.plot(radius * np.cos(theta) + center,
            radius * np.sin(theta), linewidth=1.6)
    # eixos
    ax.axvline(0, color='k', linewidth=1.2)
    ax.axhline(0, color='k', linewidth=1.2)
    # punts
    ax.plot([sigma_A, sigma_B], [0,0], 'o', markersize=8)
    # etiquetes
    y_lim = ax.get_ylim()
    gap = (y_lim[1] - y_lim[0]) / 24
    ax.text(sigma_A, gap, r'$\sigma_A$', fontsize=14)
    ax.text(sigma_B, gap, r'$\sigma_B$', fontsize=14)
    ax.set_xlabel(r'$\sigma$ [MPa]')
    ax.set_ylabel(r'$-\tau$ [MPa]')
    ax.set_title(f'Cercle de Mohr del punt {desc_punt}')
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True)
    plt.tight_layout()
    plt.show()

# %% Dades del problema (units: mm i N)

l1 = 105.0
l2 = 65.0
lA = 80.0
lB = 45.0

h = 22.0
b = 10.0

F1 = 1000.0

# Sollicitacions a la seccio A-A
N = 0.0
V_y = F1
V_z = 0.0
M_x = lA * F1        # N*mm
M_y = 0.0
T = lB * F1          # N*mm (torsor)

# Propietats de la seccio A-A (rectangular, quantitats en mm)
I = (1.0/12.0) * b * h**3  # mm^4

eta1, eta2, eta3 = f_coef_torsio_rec(h / b)

# %% Estudi punts a, b, c

sig_A = np.zeros(3)
sig_B = np.zeros(3)
tau_max = np.zeros(3)
sig_VM = np.zeros(3)

# Punt a (superior)
sig_M = M_x * (h/2.0) / I   # N/mm^2 = MPa
tau_T = 0.0
sig_A[0], sig_B[0], tau_max[0], sig_VM[0] = f_invariants2D(sig_M, 0.0, tau_T)

# Punt b (inferior?): moment invers i torsor
sig_M = -M_x * (h/2.0) / I
tau_T = eta1/eta2 * T / (b**2 * h)  # torsional shear estimate (N/mm^2)
sig_A[1], sig_B[1], tau_max[1], sig_VM[1] = f_invariants2D(sig_M, 0.0, tau_T)

# Punt c (centre pàgina) - maxima de torsor
sig_M = 0.0
tau_T = 1.0/eta2 * T / (b**2 * h)
sig_A[2], sig_B[2], tau_max[2], sig_VM[2] = f_invariants2D(sig_M, 0.0, tau_T)

punts = ['a', 'b', 'c']
for i in range(3):
    print(f"Punt {punts[i]}: sigma_A = {sig_A[i]:.3f} MPa, sigma_B = {sig_B[i]:.3f} MPa, tau_max = {tau_max[i]:.3f} MPa, sigma_VM = {sig_VM[i]:.3f} MPa")

# Dibuixar cercles de Mohr
for i in range(3):
    f_cercle_Mohr2D(sig_A[i], sig_B[i], punts[i])