#%% Suport amb càrrega rotativa: anàlisi a fatiga (solució parcial)
# Autor: J.A. Mayugo, 2026
# Cas base: només càrrega rotativa. Falta afegir la càrrega constant addicional i
# analitzar els punts més crítics a tracció i a compressió.

import numpy as np
from f_fatigue_modified_Se import f_fatigue_modified_Se

#%% Dades del problema
d   = 25.0          # mm, diàmetre secció petita
D   = 31.0          # mm, diàmetre secció gran
r   = 3.0           # mm, radi d'acord
y   = 76 + 7 - 3    # mm, distància entre la força i la secció d'estudi
F   = 1500.0        # N, força rotativa
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


#%% Tensió alterna al punt més carregat de la perifèria
M       = F*y                           # N·mm, moment flector a la secció d'estudi
sigma_0 = 32*M/(np.pi*d**3)             # MPa, tensió nominal
K_t = K_t_rodo_flexio(D/d, r/d)
q   = q_sensibilitat_entalla(r, Sut)
K_f = 1 + q*(K_t - 1)
sigma_a = K_f*sigma_0                   # MPa, amplitud de la tensió alterna

#%% Resistència: límit de fatiga modificat
Se, ka, kb, kc = f_fatigue_modified_Se('ground', d, 'bending', Sut, 'steel')

#%% Coeficients de seguretat (només tensió alterna)
SF_f = Se/sigma_a                       # a fatiga
SF_e = Sy/sigma_a                       # estàtic

print(f"M = {M:.0f} N·mm, sigma_0 = {sigma_0:.1f} MPa")
print(f"K_t = {K_t:.3f}, q = {q:.3f}, K_f = {K_f:.3f}, sigma_a = {sigma_a:.1f} MPa")
print(f"Se = {Se:.1f} MPa (ka = {ka:.3f}, kb = {kb:.3f}, kc = {kc:.1f})")
print(f"SF a fatiga = {SF_f:.2f}, SF estàtic = {SF_e:.2f}")
