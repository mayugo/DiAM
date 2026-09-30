def K_t_rodo_axial(D_d, r_d, i_flag=1):
# Calcula K_t d'un rodó sotmès a tracció amb canvi de secció D -> d.
# D = diàmetre gran; d= diàmetre petit; r=radi de l'entalla
# D_d : Relació D/d; r_d : Relació r/D
# i_flag : =1 -> K_t per tensió normal,  !=1 -> K_t per von Mises
    if i_flag == 1:
        K_t = (  0.493 + 0.480 * D_d ** (-2.43) + r_d ** (-0.48)
            *( (3.43 - 3.41 * D_d ** 2 - 0.0232 * D_d ** 4) / (1 - 8.85 * D_d ** 2 - 0.078 * D_d ** 4) )**0.5 )
    else:
        K_t = (  0.496 + 0.472 * D_d ** (-2.85) + r_d ** (-0.48)
            *( (2.921 - 2.945 * D_d ** 2 - 0.0217 * D_d ** 4) / (1 - 9.59 * D_d ** 2 - 0.053 * D_d ** 4) )**0.5  )
    return K_t