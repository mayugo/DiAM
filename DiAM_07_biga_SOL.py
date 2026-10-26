#%% Disseny d'una biga en voladís de rigidesa donada
# J.A. Mayugo, 2026
#
# Problema amb una sola variable de disseny: el cantell h d'una biga en voladís de secció
# rectangular, de longitud L i amplada b conegudes, perquè tingui exactament la rigidesa k_req.
# La solució analítica s'obté aïllant h. Aquí es busca també per prova i error:
#   Apartat 2): mètode de Montecarlo (cerca aleatòria)
#   Apartat 3): resolució de l'equació: bisecció, Newton i Brent (funció de biblioteca)
#   Apartat 4): comprovació del coeficient de seguretat i cas sense solució

import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


#%% Funcions emprades
# (en Python s'han de definir abans de cridar-les; en MatLab van al final)

def rigidesa(h, L, b, E):
    """Rigidesa a l'extrem lliure d'una biga en voladís de secció rectangular.
    Accepta escalars o vectors de NumPy (h).
    Entrada: h, L, b [mm], E [MPa]
    Sortida: k [N/mm]
    """
    I = b*h**3/12                               # moment d'inèrcia de la secció
    return 3*E*I/L**3                           # k = 3 E I / L^3


def resistencia(h, L, b, F, S_y):
    """Tensió màxima i coeficient de seguretat a fluència de la biga en voladís.
    Entrada: h, L, b [mm], F [N], S_y [MPa]
    Sortida: sigma [MPa], SF [-]
    """
    I     = b*h**3/12                           # moment d'inèrcia de la secció
    M     = F*L                                 # moment flector a l'encastament
    sigma = M*(h/2)/I                           # tensió màxima (fibra extrema)
    return sigma, S_y/sigma                     # coeficient de seguretat SF = S_y / sigma


def derivada_rigidesa(h, L, b, E):
    """Derivada de la rigidesa respecte del cantell, dk/dh = 3 E b h^2 / (4 L^3) [N/mm^2]."""
    return 3*E*b*h**2/(4*L**3)


def error_disseny(h, L, b, E, k_req):
    """Funció error: quadrat de l'error relatiu de la rigidesa, f = r^2."""
    r = (rigidesa(h, L, b, E) - k_req)/k_req
    return r**2


def montecarlo(n, h_max, L, b, E, k_req, rng):
    """Mètode de Montecarlo (Algorisme 1).
    Genera n valors de h uniformes a (0, h_max] i retorna el de menor error f.
    """
    h = h_max*(1 - rng.random(n))                       # uniforme a (0, h_max]
    f = error_disseny(h, L, b, E, k_req)
    i_opt = np.argmin(f)
    return dict(h=h, f=f, h_opt=h[i_opt], f_opt=f[i_opt])


def biseccio(h_lo, h_hi, tol, L, b, E, k_req):
    """Mètode de la bisecció (Algorisme 2).
    Busca h tal que k(h) = k_req a l'interval [h_lo, h_hi], on k(h) - k_req canvia de signe.
    S'atura quan l'interval és més petit que tol [mm].
    Retorna la solució, el nombre d'avaluacions i l'historial d'iteracions.
    """
    g = lambda h: rigidesa(h, L, b, E) - k_req               # equació g(h) = 0
    g_lo, g_hi = g(h_lo), g(h_hi)
    n_aval = 2
    assert g_lo*g_hi < 0, "L'interval inicial no conté l'arrel"
    historial = []
    while (h_hi - h_lo) > tol:
        h_m = (h_lo + h_hi)/2                           # punt mig
        g_m = g(h_m); n_aval += 1
        historial.append((h_m, g_m + k_req))
        if g_lo*g_m > 0:
            h_lo, g_lo = h_m, g_m                       # l'arrel és a la meitat superior
        else:
            h_hi, g_hi = h_m, g_m                       # l'arrel és a la meitat inferior
    return dict(h_opt=(h_lo + h_hi)/2, n_aval=n_aval, n_iter=len(historial),
                historial=np.array(historial))


def newton(h0, tol, L, b, E, k_req):
    """Mètode de Newton (Algorisme 3).
    Busca h tal que k(h) = k_req a partir del punt inicial h0.
    S'atura quan el pas és més petit que tol [mm].
    Retorna la solució, el nombre d'iteracions i l'historial d'iteracions.
    """
    h = h0
    historial = [(h, rigidesa(h, L, b, E))]
    for it in range(1, 101):
        g  = rigidesa(h, L, b, E) - k_req               # equació g(h) = 0
        dg = derivada_rigidesa(h, L, b, E)              # derivada g'(h)
        h_nou = h - g/dg                                # pas de Newton
        historial.append((h_nou, rigidesa(h_nou, L, b, E)))
        pas, h = abs(h_nou - h), h_nou
        if pas < tol:
            break
    return dict(h_opt=h, n_iter=it, n_aval=it, historial=np.array(historial))


def brent(h_lo, h_hi, tol, L, b, E, k_req):
    """Mètode de Brent amb la funció brentq de scipy (no es programa pas a pas).
    Retorna la solució, el nombre d'iteracions i el d'avaluacions de g.
    """
    g = lambda h: rigidesa(h, L, b, E) - k_req
    h, info = brentq(g, h_lo, h_hi, xtol=tol, full_output=True)
    return dict(h_opt=h, n_iter=info.iterations, n_aval=info.function_calls)


#%% Dades del problema
# Geometria [mm]
L     = 200.0       # longitud de la biga en voladís
b     = 20.0        # amplada de la secció
# Material (acer) [MPa]
E     = 200e3       # mòdul d'elasticitat
S_y   = 380.0       # límit elàstic
# Càrrega [N] i requisits
F     = 200.0       # força a l'extrem lliure
k_req = 125.0       # N/mm, rigidesa que ha de tenir la biga
SF_min = 2.0        # coeficient de seguretat mínim a fluència
# Espai de cerca [mm]
h_min = 1.0         # cantell mínim considerat
h_max = 30.0        # cantell màxim considerat
# Mètode de Montecarlo
n_MC  = [10**2, 10**3, 10**4]                   # nombre de realitzacions
rng   = np.random.default_rng(2030)             # llavor fixa per reproduir els resultats
# Mètode de la bisecció
tol   = 1e-2        # mm, mida de l'interval per aturar la cerca


#%% Apartat 1): solució analítica (per comparar)
h_an = (4*k_req*L**3/(E*b))**(1/3)                  # cantell que dona la rigidesa k_req
k_an = rigidesa(h_an, L, b, E)
sigma_an, SF_an = resistencia(h_an, L, b, F, S_y)
h_sigma = np.sqrt(6*F*L*SF_min/(b*S_y))             # cantell mínim per al coeficient de seguretat
def zona_no_factible(h_sigma, etiqueta=True):
    """Ombreja la zona de cantells que no compleixen SF >= SF_min (h < h_sigma)."""
    plt.axvspan(0, h_sigma, color='tab:red', alpha=0.12, linewidth=0,
                label='Zona no factible ($SF < SF_{\\min}$)' if etiqueta else None)
    plt.axvline(h_sigma, color='tab:red', linestyle='-', linewidth=0.9, alpha=0.8)

print(f"Solució analítica: h = {h_an:.3f} mm, k = {k_an:.2f} N/mm, "
      f"sigma = {sigma_an:.1f} MPa, SF = {SF_an:.2f}")
print(f"Cantell mínim per SF >= {SF_min}: h = {h_sigma:.3f} mm "
      f"({'governa la rigidesa' if h_an >= h_sigma else 'governa la resistència'})\n")

#%% Representació gràfica: rigidesa en funció del cantell
h_ = np.linspace(h_min, h_max, 400)
k_h = rigidesa(h_, L, b, E)
plt.figure(figsize=(4.0, 3))
zona_no_factible(h_sigma)
plt.plot(h_, k_h, '-', color='0.35', linewidth=1.4, label='$k(h)$')
plt.axhline(k_req, linestyle='--', color='0.35', linewidth=1.0, label='$k$ requerida')
plt.plot(h_an, k_req, 'b+', markersize=14, markeredgewidth=1.8, label='Solució analítica')
plt.xlim(0, h_max); plt.ylim(0, 500)
plt.xlabel('$h$ [mm]'); plt.ylabel('$k$ [N/mm]')
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=9)
plt.grid(True, alpha=0.4)
plt.savefig('DiAM_07_biga_analitica_Fig0.pdf', bbox_inches='tight')
plt.show()


#%% Apartat 2): mètode de Montecarlo
print("Apartat 2): mètode de Montecarlo")
print(f"{'n':>9} {'h [mm]':>9} {'k [N/mm]':>9} {'error h':>9} {'error k':>9} {'temps [ms]':>11}")
montecarlo(1000, h_max, L, b, E, k_req, rng)            # escalfament (no es compta)
for n in n_MC:
    t0 = time.perf_counter()
    R = montecarlo(n, h_max, L, b, E, k_req, rng)
    dt = (time.perf_counter() - t0)*1e3                     # ms
    k_R = rigidesa(R['h_opt'], L, b, E)
    err_h = (R['h_opt'] - h_an)/h_an*100
    err_k = (k_R - k_req)/k_req*100
    print(f"{n:9d} {R['h_opt']:9.4f} {k_R:9.3f} {err_h:8.3f}% {err_k:8.3f}% {dt:11.1f}")
R_vis = montecarlo(1000, h_max, L, b, E, k_req, rng)        # per al gràfic


#%% Representació gràfica: error de cada valor aleatori de h
plt.figure(figsize=(4, 3))
zona_no_factible(h_sigma)
fact = R_vis['h'] >= h_sigma
plt.semilogy(R_vis['h'][fact], R_vis['f'][fact], '.', markersize=3, color='tab:green', label='Dissenys factibles')
plt.semilogy(R_vis['h'][~fact], R_vis['f'][~fact], '.', markersize=3, color='tab:red', label='Dissenys no factibles')
plt.axvline(h_an, linestyle='--', color='k', linewidth=1.0)
plt.plot(R_vis['h_opt'], R_vis['f_opt'], 'x', color='k', markersize=10, markeredgewidth=2.0, label='Montecarlo')
#  plt.plot(h_an, 1e-9, 'b+', markersize=14, markeredgewidth=1.8, label='Solució analítica')
plt.xlim(0, h_max); plt.ylim(1e-9, 1e2)
plt.xlabel('$h$ [mm]'); plt.ylabel('Error $f$ [-]')
leg = plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=9)
for hd in leg.legend_handles[1:3]:                           # punts de la legenda més visibles
    hd.set_markersize(8)
plt.grid(True, alpha=0.4)
plt.savefig('DiAM_07_biga_MC_Fig1.pdf', bbox_inches='tight')
plt.show()


#%% Apartat 3): mètode de la bisecció
R_bis = biseccio(h_min, h_max, tol, L, b, E, k_req)         # (la primera crida és d'escalfament)
t0 = time.perf_counter()
R_bis = biseccio(h_min, h_max, tol, L, b, E, k_req)
dt_bis = (time.perf_counter() - t0)*1e3                     # ms
k_bis = rigidesa(R_bis['h_opt'], L, b, E)
err_h_bis = (R_bis['h_opt'] - h_an)/h_an*100
print("\nApartat 3): mètode de la bisecció")
print(f"  interval inicial: [{h_min:.1f}, {h_max:.1f}] mm, tolerància = {tol:.0e} mm")
print(f"  solució: h = {R_bis['h_opt']:.6f} mm, k = {k_bis:.6f} N/mm (error h {err_h_bis:.1e}%)")
print(f"  {R_bis['n_iter']} iteracions, {R_bis['n_aval']} avaluacions del model, {dt_bis:.2f} ms")
print("\n  Iteracions:")
for it, (hh, kk) in enumerate(R_bis['historial'], start=1):
    print(f"  {it:3d}  h = {hh:9.4f} mm  k = {kk:8.3f} N/mm  error k = {(kk - k_req)/k_req*100:7.2f}%")


#%% Representació gràfica: primeres iteracions de la bisecció
plt.figure(figsize=(4, 3))
zona_no_factible(h_sigma)
plt.plot(h_, k_h, '-', color='0.35', linewidth=1.4, label='$k(h)$')
plt.axhline(k_req, linestyle='--', color='0.35', linewidth=1.0, label='$k$ requerida')
Hb = R_bis['historial'][:6]
plt.plot(Hb[:, 0], Hb[:, 1], 'o', color='tab:red', markersize=5, label='Iteracions')
for it, (hh, kk) in enumerate(Hb, start=1):
    plt.annotate(str(it), (hh, kk), textcoords='offset points', xytext=(5, -12), fontsize=8, color='tab:red')
plt.plot(h_an, k_req, 'b+', markersize=14, markeredgewidth=1.8, label='Solució analítica')
plt.xlim(0, h_max); plt.ylim(0, 500)
plt.xlabel('$h$ [mm]'); plt.ylabel('$k$ [N/mm]')
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=9)
plt.grid(True, alpha=0.4)
plt.savefig('DiAM_07_biga_biseccio_Fig2.pdf', bbox_inches='tight')
plt.show()


#%% Apartat 3): mètode de Newton
h0_N = h_max                                            # punt inicial
R_new = newton(h0_N, tol, L, b, E, k_req)               # (la primera crida és d'escalfament)
t0 = time.perf_counter()
R_new = newton(h0_N, tol, L, b, E, k_req)
dt_new = (time.perf_counter() - t0)*1e3                 # ms
k_new = rigidesa(R_new['h_opt'], L, b, E)
err_h_new = (R_new['h_opt'] - h_an)/h_an*100
print("\nApartat 3): mètode de Newton")
print(f"  punt inicial: h0 = {h0_N:.1f} mm, tolerància = {tol:.0e} mm")
print(f"  solució: h = {R_new['h_opt']:.6f} mm, k = {k_new:.6f} N/mm (error h {err_h_new:.1e}%)")
print(f"  {R_new['n_iter']} iteracions, {R_new['n_aval']} avaluacions de g i {R_new['n_aval']} de g', {dt_new:.3f} ms")
print("\n  Iteracions:")
for it, (hh, kk) in enumerate(R_new['historial'][:8]):
    print(f"  {it:3d}  h = {hh:9.4f} mm  k = {kk:8.3f} N/mm  error k = {(kk - k_req)/k_req*100:7.2f}%")
R_new2 = newton(h_min, tol, L, b, E, k_req)             # punt inicial llunyà: h0 = h_min
print(f"  amb h0 = {h_min:.1f} mm: {R_new2['n_iter']} iteracions "
      f"(primer pas fins a h = {R_new2['historial'][1, 0]:.1f} mm)")


#%% Representació gràfica: primeres iteracions de Newton (tangents)
plt.figure(figsize=(4, 3))
zona_no_factible(h_sigma)
plt.plot(h_, k_h, '-', color='0.35', linewidth=1.4, label='$k(h)$')
plt.axhline(k_req, linestyle='--', color='0.35', linewidth=1.0, label='$k$ requerida')
Hn = R_new['historial'][:5]
for j in range(3):                                      # tangent de les tres primeres iteracions
    h_j, k_j = Hn[j]
    dk = derivada_rigidesa(h_j, L, b, E)
    h_tan = np.array([Hn[j + 1, 0], h_j])
    plt.plot(h_tan, k_j + dk*(h_tan - h_j), '-', color='tab:blue', linewidth=0.9,
             label='Tangent' if j == 0 else None)
    plt.plot([Hn[j + 1, 0]]*2, [k_req, Hn[j + 1, 1]], ':', color='tab:blue', linewidth=0.9)
plt.plot(Hn[:, 0], Hn[:, 1], 'o', color='tab:red', markersize=5, label='Iteracions')
for j, (hh, kk) in enumerate(Hn[:5]):
    plt.annotate(str(j), (hh, kk), textcoords='offset points', xytext=(-10, 3), fontsize=8, color='tab:red')
plt.plot(h_an, k_req, 'b+', markersize=14, markeredgewidth=1.8, label='Solució analítica')
plt.xlim(5, h_max); plt.ylim(0, 1200)
plt.xlabel('$h$ [mm]'); plt.ylabel('$k$ [N/mm]')
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=9)
plt.grid(True, alpha=0.4)
plt.savefig('DiAM_07_biga_newton_Fig3.pdf', bbox_inches='tight')
plt.show()


#%% Apartat 3): mètode de Brent (brentq)
R_bre = brent(h_min, h_max, tol, L, b, E, k_req)        # (la primera crida és d'escalfament)
t0 = time.perf_counter()
R_bre = brent(h_min, h_max, tol, L, b, E, k_req)
dt_bre = (time.perf_counter() - t0)*1e3                 # ms
k_bre = rigidesa(R_bre['h_opt'], L, b, E)
err_h_bre = (R_bre['h_opt'] - h_an)/h_an*100
print("\nApartat 3): mètode de Brent (brentq)")
print(f"  solució: h = {R_bre['h_opt']:.6f} mm, k = {k_bre:.6f} N/mm (error h {err_h_bre:.1e}%)")
print(f"  {R_bre['n_iter']} iteracions, {R_bre['n_aval']} avaluacions de g, {dt_bre:.3f} ms")


#%% Apartat 4): comprovació del coeficient de seguretat; casos amb i sense solució
print("\nApartat 4): casos amb i sense solució")
print(f"{'cas':>4} {'F [N]':>7} {'h_k [mm]':>9} {'h_sigma [mm]':>13} {'SF(h_k)':>8} {'k(h_sigma)':>11}  resultat")
for cas, F_cas in [('A', 200.0), ('B', 400.0)]:
    h_k = h_an                                          # el cantell que dona k_req no depèn de F
    h_s = np.sqrt(6*F_cas*L*SF_min/(b*S_y))             # cantell mínim per al coeficient de seguretat
    _, SF_k = resistencia(h_k, L, b, F_cas, S_y)
    k_s = rigidesa(h_s, L, b, E)
    ok = SF_k >= SF_min
    print(f"{cas:>4} {F_cas:7.0f} {h_k:9.3f} {h_s:13.3f} {SF_k:8.2f} {k_s:11.1f}  "
          f"{'té solució' if ok else 'no té solució'}")


#%% Representació gràfica: comparació dels casos A i B
casos = [('A', 200.0, 'tab:red', 0.12), ('B', 400.0, 'tab:red', 0.12)]
h_sA = np.sqrt(6*200.0*L*SF_min/(b*S_y))
h_sB = np.sqrt(6*400.0*L*SF_min/(b*S_y))
hh = np.linspace(h_min, h_max, 400)
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.6))
# (a) rigidesa: la solució és la mateixa; canvia la zona factible
ax[0].axvspan(0, h_sA, color='tab:red', alpha=0.12, linewidth=0, label='Zona no factible, cas A')
ax[0].axvspan(0, h_sB, color='tab:red', alpha=0.12, linewidth=0, label='Zona no factible, cas B')
ax[0].axvline(h_sA, color='tab:red', linestyle='-', linewidth=0.9)
ax[0].axvline(h_sB, color='tab:red', linestyle='--', linewidth=0.9)
ax[0].plot(hh, rigidesa(hh, L, b, E), '-', color='0.35', linewidth=1.4, label='$k(h)$')
ax[0].axhline(k_req, linestyle='--', color='0.35', linewidth=1.0, label='$k$ requerida')
ax[0].plot(h_an, k_req, 'b+', markersize=14, markeredgewidth=1.8, label='Solució analítica')
ax[0].set_xlim(0, 20); ax[0].set_ylim(0, 500)
ax[0].set_xlabel('$h$ [mm]'); ax[0].set_ylabel('$k$ [N/mm]')
ax[0].grid(True, alpha=0.4); ax[0].set_title('(a) Rigidesa', fontsize=10)
# (b) coeficient de seguretat per a les dues forces
for cas, F_cas, _, _ in casos:
    SFh = S_y/(6*F_cas*L/(b*hh**2))
    ax[1].plot(hh, SFh, '-' if cas == 'A' else '--', color='0.35', linewidth=1.4, label=f'$SF(h)$, cas {cas} ($F = {F_cas:.0f}$ N)')
ax[1].axhline(SF_min, color='tab:red', linewidth=1.0, label='$SF_{\\min}$')
ax[1].axvline(h_an, linestyle='--', color='k', linewidth=1.0)
for cas, F_cas, _, _ in casos:
    _, SFk = resistencia(h_an, L, b, F_cas, S_y)
    ax[1].plot(h_an, SFk, 'b+', markersize=14, markeredgewidth=1.8, label='Solució analítica' if cas == 'A' else None)
    ax[1].annotate(cas, (h_an, SFk), textcoords='offset points', xytext=(-14, 4), fontsize=9)
ax[1].set_xlim(0, 20); ax[1].set_ylim(0, 8)
ax[1].set_xlabel('$h$ [mm]'); ax[1].set_ylabel('$SF$ [-]')
ax[1].grid(True, alpha=0.4); ax[1].set_title('(b) Coeficient de seguretat', fontsize=10)
h0_, l0_ = ax[0].get_legend_handles_labels()
h1_, l1_ = ax[1].get_legend_handles_labels()
vist = {}
for hd, lb in zip(h0_ + h1_, l0_ + l1_):
    vist.setdefault(lb, hd)                                  # sense etiquetes repetides
fig.tight_layout()
fig.legend(vist.values(), vist.keys(), loc='upper left', bbox_to_anchor=(1.0, 0.95), fontsize=9)
fig.savefig('DiAM_07_biga_casos_Fig4.pdf', bbox_inches='tight')
plt.show()
