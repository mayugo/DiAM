# Biga corba vs biga recta a flexió (solució analítica)
import numpy as np
import matplotlib.pyplot as plt

# Check if Python is called outside optiSLang
if not 'OSL_REGULAR_EXECUTION' in locals():
    OSL_REGULAR_EXECUTION = False

# values to default if not existing yet
if not OSL_REGULAR_EXECUTION: # test run mode

## Dades
    F = 4500        # N, forces verticals
    d = 40          # mm, distància entre forces
    h = 40          # mm, cantell biga
    b = 20          # mm, gruix biga
    ro = 70         # mm, radi exterior

M = F*d         # N·mm, moment flector

## Propietats de la secció
A  = h*b                # mm^2
I  = b*h**3/12          # mm^4
ri = ro - h             # mm, radi interior
rc = (ro + ri)/2        # mm, radi del centroide
rn = h/np.log(ro/ri)    # mm, radi de la línia neutra
e  = rc - rn            # mm, excentricitat

## Tensions (path: 0 = fibra exterior, h = fibra interior)
path = np.linspace(0, h, 200)

y_recta = path - h/2            # des del centroide
y_corba = path - (h/2 + e)      # des de la línia neutra

sigma_recta = M*y_recta/I
sigma_corba = M*y_corba/(A*e*(rn - y_corba))

## Posició de les fibres neutres i tensions als extrems
path_n_recta = h/2              # mm des de la fibra exterior
path_n_corba = h/2 + e

ext = [0, -1]                   # índexs: fibra exterior, fibra interior
s_rec = sigma_recta[ext]
s_cor = sigma_corba[ext]
dif_pct = (np.abs(s_cor) - np.abs(s_rec))/np.abs(s_rec)*100   # > 0: corba més tensionada

sigma_max = s_cor[1]

if not OSL_REGULAR_EXECUTION: # test run mode

    print(f'Desplaçament fibra neutra (e) = {e:.2f} mm')
    print(f'Fibra exterior: recta {s_rec[0]:.1f} MPa | corba {s_cor[0]:.1f} MPa | {dif_pct[0]:+.1f} %')
    print(f'Fibra interior: recta {s_rec[1]:.1f} MPa | corba {s_cor[1]:.1f} MPa | {dif_pct[1]:+.1f} %')
    
    ## Gràfic
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.grid(True)
    ax.tick_params(labelsize=14)
    
    c1, c2 = 'C0', 'C1'
    ax.plot(sigma_recta, path, color=c1, lw=1.8, label='Biga recta')
    ax.plot(sigma_corba, path, color=c2, lw=1.8, label='Biga corba')
    ax.plot([0, 0], [0, h], 'k')
    
    xm = 2.3*np.max(np.abs(np.concatenate([sigma_recta, sigma_corba])))
    ax.set_xlim(-xm, xm)
    ax.set_ylim(-3, h + 3)
    
    # fibres neutres
    ax.plot([-xm, xm], [path_n_recta]*2, '--', color=c1)
    ax.plot([-xm, xm], [path_n_corba]*2, '--', color=c2)
    ax.text(-0.98*xm, path_n_recta, f'F. neutra recta\n(y = {path_n_recta:.1f} mm)',
            color=c1, va='top', fontsize=14)
    ax.text(-0.98*xm, path_n_corba, f'F. neutra corba\n(y = {path_n_corba:.1f} mm, e = {e:.2f} mm)',
            color=c2, va='bottom', fontsize=14)
    
    # extrems
    ax.plot(s_rec, path[ext], 'o', color=c1)
    ax.plot(s_cor, path[ext], 'o', color=c2)
    
    def txt(i):
        return (f'recta: {s_rec[i]:.1f} MPa\n'
                f'corba: {s_cor[i]:.1f} MPa\n'
                f'$\\Delta$ = {dif_pct[i]:+.1f} %')
    
    ax.text(min(s_rec[0], s_cor[0]) - 0.03*xm, path[0], txt(0), ha='right', va='center', fontsize=14)
    ax.text(max(s_rec[1], s_cor[1]) + 0.03*xm, path[-1], txt(1), ha='left', va='center', fontsize=14)
    
    ax.set_xlabel(r'$\sigma_L$ [MPa]', fontsize=18)
    ax.set_ylabel(r'$y$ [mm]  (0 = fibra exterior)', fontsize=16)
    ax.legend(loc='upper left', fontsize=16)
    
    plt.tight_layout()
    plt.show()