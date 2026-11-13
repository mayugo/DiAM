import numpy as np

def f_fatiga_limitada(sigma_m, sigma_a, S_ut, S_0, N_0, S_e, N_e):
    """Calcula la vida esperada per a un estat de carrega equivalent segons la
    funcio de Basquin: Seqv = a N^b
      sigma_m = tensio equivalent mitjana (nomes per sigma_m positives)
      sigma_a = tensio equivalent alterna
      S_ut    = resistencia a traccio
      S_0     = resistencia a N_0, per R=-1
      N_0     = cicles
      S_e     = limit de fatiga a N_e, per R=-1
      N_e     = cicles
    Retorna (N, mode): nombre de cicles i descripcio del mode de falla
    """
    if sigma_m < 0:                   # comprova que tensio mitjana es positiva
        return 0, 'error: tensio mitjana negativa'

    sigma_max = sigma_m + sigma_a
    # A = sigma_a/sigma_m = (1-R)/(1+R)
    A = np.inf if sigma_m == 0 else sigma_a/sigma_m
    if np.isinf(A):
        Sfe = S_e
        Sf0 = S_0
    else:
        Sfe = (1 + A)*S_ut/(1 + A*S_ut/S_e)
        Sf0 = (1 + A)*S_ut/(1 + A*S_ut/S_0)

    if sigma_max > S_ut:
        N = 0
        mode = 'falla estatica'
    elif sigma_max > Sf0:
        b_1 = (np.log10(S_ut) - np.log10(Sf0))/(np.log10(1) - np.log10(N_0))
        a_1 = 10**(np.log10(S_ut) - b_1*np.log10(1))
        N = float((sigma_max/a_1)**(1/b_1))
        mode = 'baix #cicles: vida limitada'
    elif sigma_max > Sfe:
        b_2 = (np.log10(Sf0) - np.log10(Sfe))/(np.log10(N_0) - np.log10(N_e))
        a_2 = 10**(np.log10(Sf0) - b_2*np.log10(N_0))
        N = float((sigma_max/a_2)**(1/b_2))
        mode = 'alt #cicles: vida limitada'
    else:
        N = np.inf
        mode = 'alt #cicles: vida ilimitada'
    return N, mode