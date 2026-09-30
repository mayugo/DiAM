def K_t_rodo_flexio(D_d, r_d, i_flag=1):
# Calcula K_t d'un rodó sotmès a flexió amb canvi de secció D -> d.
# D = diàmetre gran; d= diàmetre petit; r=radi de l'entalla
# D_d : Relació D/d; r_d : Relació r/D
# i_flag : =1 -> K_t per tensió normal,  !=1 -> K_t per von Mises
    if i_flag == 1:
        K_t = ( 0.632 + 0.377 * D_d**(-4.4) + r_d**(-0.5)
            * ( (-0.14 - 0.363 * D_d**2 + 0.503 * D_d**4)
                / (1 - 2.39 * D_d**2 + 3.368 * D_d**4) )**0.5 )
    else:
        K_t = ( 0.622 + 0.38 * D_d**(-4.3) + r_d**(-0.5)
            * ( (-0.322 - 0.277 * D_d**2 + 0.599 * D_d**4)
                / (1 - 2.55 * D_d**2 + 5.27 * D_d**4) )**0.5 )
    return K_t
