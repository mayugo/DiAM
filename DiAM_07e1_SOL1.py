import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from DiAM_07e1_simulacio import simulacio_transitori

# ---------------------------------------------------------------
# DADES (SI: m, kg, N, rad/s)
# ---------------------------------------------------------------
w_tab = [60, 120, 45, 30]                 # rad/s
y0_tab = [1.5, 0.9, 1.4, 1.1]             # mm
m_tab = [1.6, 1.2, 1.8, 2.0]              # kg
X_lim = 5e-3                              # amplitud desitjada de x
zeta = 0.02
g = 9.81


# ---------------------------------------------------------------
# FORMULACIO
# ---------------------------------------------------------------
def transmitancia(Omega, zeta):
    """Transmitancia de desplacament X/y0 amb esmorteiment."""
    num = np.sqrt(1 + (2*zeta*Omega)**2)
    den = np.sqrt((1 - Omega**2)**2 + (2*zeta*Omega)**2)
    return num/den


def rigidesa_necessaria(omega, m, zeta, X, y0):
    """Rigidesa total k tal que la transmitancia sigui X/y0 amb Omega>1."""
    Omega = brentq(lambda Omega: transmitancia(Omega, zeta) - X/y0, 1.0, 100.0)
    wn = omega/Omega
    return m*wn**2, Omega


def forces_suport(omega, m, zeta, k, X, y0):
    """Forces sobre cada suport: mitjana (pes) i alterna (dinamica)."""
    F_mitja = m*g/2
    F_alt = m*omega**2*X/2          # forca d'inercia del motlle / 2
    # comprovacio: forca de la molla i l'esmorteidor, a partir de x-y
    Omega = omega/np.sqrt(k/m)
    Z = y0*Omega**2/np.sqrt((1 - Omega**2)**2 + (2*zeta*Omega)**2)
    c = 2*zeta*np.sqrt(k*m)
    F_alt_comp = np.sqrt(k**2 + (omega*c)**2)*Z/2
    return F_mitja, F_alt, F_alt_comp


# ---------------------------------------------------------------
# CALCULS I SIMULACIO PER A LES 4 CONFIGURACIONS
# ---------------------------------------------------------------
for cfg in range(1, 5):
    omega = w_tab[cfg-1]
    y0 = y0_tab[cfg-1]*1e-3
    m = m_tab[cfg-1]

    k, Omega = rigidesa_necessaria(omega, m, zeta, X_lim, y0)
    wn = np.sqrt(k/m)
    Tr = transmitancia(Omega, zeta)
    F_mitja, F_alt, F_alt_comp = forces_suport(omega, m, zeta, k, X_lim, y0)

    t, x, y, F, X_ss, F_ss = simulacio_transitori(m, k, zeta, y0, omega)

    print(f'\nConfiguracio {cfg}: w = {omega} rad/s, y0 = {y0*1e3:.2f} mm, m = {m} kg')
    print(f'  Omega = {Omega:.3f}   wn = {wn:.1f} rad/s   T = {Tr:.2f}')
    print(f'  k total = {k/1e3:.2f} N/mm   k per suport = {k/2e3:.2f} N/mm')
    print(f'  F_m = {F_mitja:.2f} N   F_a = {F_alt:.2f} N (comprovacio {F_alt_comp:.2f} N)')
    print(f'  F_max = {F_mitja + F_alt:.2f} N   F_min = {F_mitja - F_alt:.2f} N')
    print(f'  Simulacio: X_ss = {X_ss*1e3:.2f} mm, max|x| transitori = {np.max(np.abs(x))*1e3:.2f} mm')
    print(f'             F_ss = {F_ss:.2f} N, max|F| transitori = {np.max(np.abs(F)):.2f} N, '
          f'F total max/min = {F_mitja + np.max(F):.2f} / {F_mitja + np.min(F):.2f} N')

    # figura nomes per a la configuracio 1
    if cfg == 1:
        plt.rcParams.update({'font.size': 12})
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
        ax1.plot(t, x*1e3, label='x (motlle)')
        ax1.plot(t, y*1e3, label='y (base)', lw=0.8)
        ax1.axhline(X_lim*1e3, color='k', ls='--', lw=0.8)
        ax1.axhline(-X_lim*1e3, color='k', ls='--', lw=0.8)
        ax1.set_ylabel('desplacament [mm]')
        ax1.legend(loc='upper right')
        ax1.grid(alpha=0.3)
        ax2.plot(t, F_mitja + F, color='C2')
        ax2.axhline(F_mitja, color='k', ls='--', lw=0.8)
        ax2.set_xlabel('t [s]')
        ax2.set_ylabel('F per suport [N]')
        ax2.set_xlim(0, 4)
        ax2.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig('DiAM_07e1_SOL1_transitori.pdf')
