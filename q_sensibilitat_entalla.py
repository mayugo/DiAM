import numpy as np

def q_sensibilitat_entalla(r, Sut, iflag):
    # Sensibilitat a l'entalla (q) per elements d'acer a càrrega alterna a flexio, forca axial, o tallant
    #   r: radi de l'entalla, r input en [mm], formula en [in]
    #                   1 in = 25.4 mm
    #   Sut: resitencia ultima, Sut input en [MPa), formula en [kpsi]
    #                   1 kpsi = 6.9 MPa
    #   iflag = 0 flexio/axial, iflag <> 0 tallant

    Sut_kpsi = Sut / 6.9       # [kpsi]
    r_in = r / 25.4            # [in]

    if iflag == 0:
        sqr_a = ( 0.246 - 3.08e-3 * Sut_kpsi + 1.51e-5 * Sut_kpsi**2 - 2.67e-8 * Sut_kpsi**3 )
    else:
        sqr_a = ( 0.190 - 2.51e-3 * Sut_kpsi + 1.35e-5 * Sut_kpsi**2 - 2.67e-8 * Sut_kpsi**3 )

    q = 1 / (1 + sqr_a / np.sqrt(r_in))
    return q