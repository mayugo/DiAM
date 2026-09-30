def K_t_rodo_torsio(D_d, r_d):
# Calcula K_t d'un rodó sotmès a torsió amb canvi de secció D -> d.
# D = diàmetre gran; d= diàmetre petit; r=radi de l'entalla
# D_d : Relació D/d; r_d : Relació r/D
    K_t = ( 0.78 + 0.2 * D_d**(-10) + r_d**(-0.46)
        * ( (0.002 - 0.125 * D_d**2 + 0.123 * D_d**4) / (1 - 2.75 * D_d**2 + 2.55 * D_d**4) )**0.5 )
    return K_t
