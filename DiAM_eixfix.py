# Suport amb eix fix
# J.A. Mayugo 2015
import numpy as np

def K_t_rodo_flexio(D_d, r_d, i_flag=1):
# Calcula K_t d'un rodo sotmes a flexio amb canvi de seccio D -> d.
# D = diametre gran; d= diametre petit; r=radi de l'entalla
# D_d : Relacio D/d; r_d : Relacio r/D
# i_flag : =1 -> K_t per tensio normal,  !=1 -> K_t per von Mises
    if i_flag == 1:
        K_t = ( 0.632 + 0.377 * D_d**(-4.4) + r_d**(-0.5)
            * ( (-0.14 - 0.363 * D_d**2 + 0.503 * D_d**4)
                / (1 - 2.39 * D_d**2 + 3.368 * D_d**4) )**0.5 )
    else:
        K_t = ( 0.622 + 0.38 * D_d**(-4.3) + r_d**(-0.5)
            * ( (-0.322 - 0.277 * D_d**2 + 0.599 * D_d**4)
                / (1 - 2.55 * D_d**2 + 5.27 * D_d**4) )**0.5 )
    return K_t

# Dades problema
l = 60    # mm , distancia canvi seccio eix-politja
R = 100   # mm , radi politja
D = 25    # mm , diametre D eix politja
d = 20    # mm , diametre d eix politja
r = 4     # mm , R4, radi acord entre D i d

F1 = 100  # N , força corretja 1
F2 = 400  # N , força corretja 2

# Seccio A, punt a
# Calcul de sol·licitacions
M = (l - r) * (F1 + F2)   # N·mm
T = R * (F2 - F1)         # N·mm
V = F1 + F2               # N

# Calcul de factors de concentracio teorics K_t
K_t = K_t_rodo_flexio(D / d, r / d, 1)

# Càlcul de la tensio normal
sigma = (32 * M) / (np.pi * d**3) * K_t   # MPa
print(f"sigma = {sigma:.2f} MPa")