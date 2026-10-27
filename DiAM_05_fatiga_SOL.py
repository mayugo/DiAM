#%% Suport amb càrrega rotativa i càrrega constant: anàlisi a fatiga (solució completa)
# Autor: J.A. Mayugo, 2026

import numpy as np
import matplotlib.pyplot as plt
from f_fatigue_modified_Se import f_fatigue_modified_Se

plt.rcParams['svg.fonttype'] = 'none'

#%% Dades del problema
d   = 25.0          # mm, diàmetre secció petita
D   = 31.0          # mm, diàmetre secció gran
r   = 3.0           # mm, radi d'acord
y   = 76 + 7 - 3    # mm, distància entre la força i la secció d'estudi
F   = 1500.0        # N, força rotativa
Fc  = 2000.0        # N, força constant addicional
Sut = 750.0         # MPa, resistència última
Sy  = 550.0         # MPa, límit elàstic


#%% Funcions: factor de concentració teòric i sensibilitat a l'entalla
def K_t_rodo_flexio(D_d, r_d):
    """K_t d'un rodó sotmès a flexió amb canvi de diàmetre de D a d (D_d = D/d, r_d = r/d)."""
    return (0.632 + 0.377*D_d**(-4.4)
            + r_d**(-0.5)*np.sqrt((-0.14 - 0.363*D_d**2 + 0.503*D_d**4)
                                  /(1 - 2.39*D_d**2 + 3.368*D_d**4)))

def q_sensibilitat_entalla(r, Sut):
    """Sensibilitat a l'entalla q d'un acer a flexió (r en mm, Sut en MPa)."""
    Sut_kpsi = Sut/6.9
    r_in = r/25.4
    sqr_a = 0.246 - 3.08e-3*Sut_kpsi + 1.51e-5*Sut_kpsi**2 - 2.67e-8*Sut_kpsi**3
    return 1/(1 + sqr_a/np.sqrt(r_in))

def coeficients_seguretat(sigma_m, sigma_a, Se):
    """Coeficients de seguretat per a una tensió mitjana sigma_m i una amplitud sigma_a (uniaxial).
    Estàtic: fluència de von Mises amb la tensió màxima en valor absolut.
    Fatiga: Goodman modificat si sigma_m > 0; si sigma_m <= 0 no s'aprofita la compressió: Se/sigma_a."""
    sigma_max, sigma_min = sigma_m + sigma_a, sigma_m - sigma_a
    SF_e = Sy/max(abs(sigma_max), abs(sigma_min))
    SF_f = 1/(sigma_a/Se + sigma_m/Sut) if sigma_m > 0 else Se/sigma_a
    return SF_e, SF_f, sigma_max, sigma_min


#%% Tensions nominals (flexió)
M_a, M_m = F*y, Fc*y                                    # N·mm
sigma_0a = 32*M_a/(np.pi*d**3)                          # MPa, tensió nominal alterna
sigma_0m = 32*M_m/(np.pi*d**3)                          # MPa, tensió nominal mitjana


#%% Apartats 2) i 3): tensions i coeficients de seguretat
K_t = K_t_rodo_flexio(D/d, r/d)
q   = q_sensibilitat_entalla(r, Sut)
K_f = 1 + q*(K_t - 1)
Se, ka, kb, kc = f_fatigue_modified_Se('ground', d, 'bending', Sut, 'steel')
sigma_a = K_f*sigma_0a                                  # MPa
sigma_m = K_f*sigma_0m                                  # MPa (es comprova que no hi ha fluència local)
print("\nApartat 2)")
print(f"  sigma_0a = {sigma_0a:.1f} MPa, sigma_0m = {sigma_0m:.1f} MPa")
print(f"  K_t = {K_t:.3f}, q = {q:.3f}, K_f = {K_f:.3f}; Se = {Se:.1f} MPa "
      f"(ka = {ka:.3f}, kb = {kb:.3f}, kc = {kc:.1f})")
print(f"  comprovació K_f*(sigma_0a + sigma_0m) = {K_f*(sigma_0a + sigma_0m):.0f} MPa < Sy = {Sy:.0f} MPa")

print("\n  cas                      sigma_m  sigma_a  sigma_max  sigma_min   SF_e    SF_f")
casos = [("només rotativa", 0.0), ("rotativa + constant, tracció", sigma_m),
         ("rotativa + constant, compressió", -sigma_m)]
res = {}
for nom, s_m in casos:
    SF_e, SF_f, s_max, s_min = coeficients_seguretat(s_m, sigma_a, Se)
    res[nom] = (SF_e, SF_f)
    print(f"  {nom:<32} {s_m:7.1f}  {sigma_a:7.1f}  {s_max:9.1f}  {s_min:9.1f}  {SF_e:5.2f}  {SF_f:5.2f}")


#%% Apartat 3): sensibilitat a la fiabilitat (kd = 0.814 per a una fiabilitat del 99 %)
kd = 0.814
SF_f99 = 1/(sigma_a/(kd*Se) + sigma_m/Sut)
print(f"\nApartat 3) amb kd = {kd}: SF a fatiga (tracció) = {SF_f99:.2f}")


#%% Representació gràfica: diagrama de Haigh amb les línies de càrrega
# Línia de càrrega: recta de l'origen (O) al punt de treball (A); es prolonga fins al límit (B).
# SF = OB/OA: fins al límit de fatiga dona SF_f i fins a la fluència dona SF_e.
from matplotlib.lines import Line2D
sm = np.linspace(-Sy, Sut, 600)
fatiga = np.where(sm >= 0, Se*(1 - sm/Sut), Se)                      # Goodman modificat
sm_l = np.linspace(-Sy, Sy, 600)
langer = Sy - np.abs(sm_l)                                           # fluència (Langer)
plt.rcParams.update({'font.size': 12})
plt.figure(figsize=(7.0, 4.7))
OR = '#E8710A'   # taronja
BL = '#6FA3DC'   # blau pastel
plt.plot(sm, np.maximum(fatiga, 0), '-', color='0.35', linewidth=1.8)
plt.plot(sm_l, langer, '--', color='0.35', linewidth=1.2)
casos_g = [('Només rotativa', 0.0, 'o', BL),
           ('Rotativa + constant (tracció)', sigma_m, 's', OR),
           ('Rotativa + constant (compressió)', -sigma_m, '^', OR)]
handles = [Line2D([], [], color='0.35', linestyle='-', linewidth=1.8, label='Límit de fatiga'),
           Line2D([], [], color='0.35', linestyle='--', linewidth=1.2, label='Límit de fluència')]
for nom, s_m, mk, c in casos_g:
    SF_e, SF_f, _, _ = coeficients_seguretat(s_m, sigma_a, Se)
    plt.plot([0, SF_f*s_m], [0, SF_f*sigma_a], '-', color=c, linewidth=0.9)
    plt.plot([0, SF_e*s_m], [0, SF_e*sigma_a], '--', color=c, linewidth=0.9)
    plt.plot(s_m, sigma_a, mk, color=c, markersize=6)
    # desplaçament (punts) i alineació de les etiquetes segons el cas
    pos = {('f', 0): (6, 6, 'left'), ('e', 0): (0, 10, 'center'),
           ('f', 1): (12, -17, 'left'), ('e', 1): (6, 4, 'left'),
           ('f', 2): (0, 10, 'center'), ('e', 2): (-7, -3, 'right')}
    idx = [0.0, sigma_m, -sigma_m].index(s_m)
    for n, etq in ((SF_f, 'f'), (SF_e, 'e')):
        dx, dy, ha = pos[(etq, idx)]
        plt.plot(n*s_m, n*sigma_a, mk, color=c, markerfacecolor='none', markersize=7)
        plt.annotate(f'$SF_{etq}$ = {n:.2f}', (n*s_m, n*sigma_a), textcoords='offset points',
                     xytext=(dx, dy), ha=ha, va='center', fontsize=9, color=c,
                     bbox=dict(fc='white', ec='none', pad=0.6, alpha=0.85) if (etq, idx) == ('f', 1) else None)
    handles.append(Line2D([], [], color=c, marker=mk, linestyle='', markersize=6, label=nom))
ax = plt.gca()
ax.set_xlim(-Sy*1.05, Sut*1.05); ax.set_ylim(0, 1.3*Sy)
ax.set_aspect('equal', adjustable='box')   # mateixa escala a x i y (com 'axis equal'), sense canviar els límits
from matplotlib.ticker import MultipleLocator
ax.xaxis.set_major_locator(MultipleLocator(200)); ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_major_locator(MultipleLocator(200)); ax.yaxis.set_minor_locator(MultipleLocator(100))
plt.xlabel('$\\sigma_m$ [MPa]'); plt.ylabel('$\\sigma_a$ [MPa]')
plt.legend(handles=handles, loc='upper right', fontsize=9, framealpha=0.9)
ax.grid(True, which='major', alpha=0.45)
ax.grid(True, which='minor', alpha=0.2)
plt.savefig('DiAM_fatiga_Haigh_Fig1.pdf', bbox_inches='tight')
plt.savefig('DiAM_fatiga_Haigh_Fig1.svg', bbox_inches='tight')
plt.show()