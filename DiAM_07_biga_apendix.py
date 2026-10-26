#%% Apèndix: mètodes numèrics per resoldre una equació f(x) = 0
# J.A. Mayugo, 2026
#   Figura 4: idea geomètrica de la bisecció, Newton, secant i Brent (funció genèrica)
#   Figura 5: convergència dels quatre mètodes en l'equació de la rigidesa de la biga

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

plt.rcParams['svg.fonttype'] = 'none'
vermell, blau, gris = 'tab:red', 'tab:blue', '0.35'


#%% Funció genèrica per a les quatre representacions
f  = lambda x: x**3 - 5                 # funció de l'exemple, arrel x* = 5^(1/3)
df = lambda x: 3*x**2                   # derivada
x_arrel = 5**(1/3)
x_ = np.linspace(0.6, 3.2, 300)
a0, b0 = 1.0, 3.0                       # interval inicial / punts inicials

fig, ax = plt.subplots(3, 2, figsize=(8, 9.6), sharex=True)
ax = ax.ravel()
for A in ax:
    A.plot(x_, f(x_), '-', color=gris, linewidth=1.4)
    A.axhline(0, color='k', linewidth=0.8)
    A.plot(x_arrel, 0, 'b+', markersize=12, markeredgewidth=1.8)
    A.grid(True, alpha=0.4)
    A.set_xlim(0.6, 3.2)
for A in ax[4:]:
    A.set_xlabel('$x$')
for A in ax[::2]:
    A.set_ylabel('$f(x)$')

# (a) bisecció: punts mitjos i intervals successius
a, b = a0, b0
A = ax[0]; A.set_ylim(-14, 28)
for j in range(4):
    A.hlines(-6.5 - 1.8*j, a, b, color=vermell, linewidth=3)
    A.text(0.65, -6.5 - 1.8*j, f'{j}', va='center', fontsize=8, color=vermell)
    if j < 3:
        xm = (a + b)/2
        A.plot(xm, f(xm), 'o', color=vermell, markersize=5)
        A.annotate(str(j + 1), (xm, f(xm)), textcoords='offset points', xytext=(5, 3), fontsize=8, color=vermell)
        if f(a)*f(xm) > 0: a = xm
        else: b = xm
A.set_title('(a) Bisecció', fontsize=10)

# (b) regula falsi: recta pels extrems de l'interval; el tall amb l'eix substitueix l'extrem del mateix signe
A = ax[1]; A.set_ylim(-6, 28)
a, b = a0, b0
A.plot([a, b], [f(a), f(b)], 'o', color=vermell, markersize=5)
A.annotate('a', (a, f(a)), textcoords='offset points', xytext=(5, -10), fontsize=8, color=vermell)
A.annotate('b', (b, f(b)), textcoords='offset points', xytext=(5, 3), fontsize=8, color=vermell)
estils = ['-', '--', ':']
for j in range(3):
    xn = b - f(b)*(b - a)/(f(b) - f(a))
    A.plot([a, b], [f(a), f(b)], estils[j], color=blau, linewidth=1.1)
    A.plot(xn, 0, 's', color=blau, markersize=4)
    A.plot([xn, xn], [0, f(xn)], ':', color=gris, linewidth=0.8)
    A.plot(xn, f(xn), 'o', color=vermell, markersize=5)
    A.annotate(str(j + 1), (xn, f(xn)), textcoords='offset points', xytext=(5, -10), fontsize=8, color=vermell)
    if f(xn)*f(a) > 0: a = xn
    else: b = xn
A.set_title('(b) Regula falsi', fontsize=10)

# (c) Illinois: com la regula falsi, però quan el mateix extrem es manté dues iteracions seguides
# es divideix per dos el valor de f d'aquest extrem (punt buit: valor de f modificat)
A = ax[2]; A.set_ylim(-6, 28)
a, b = a0, b0
fa, fb = f(a), f(b)
costat = 0
A.plot([a, b], [fa, fb], 'o', color=vermell, markersize=5)
A.annotate('a', (a, fa), textcoords='offset points', xytext=(5, -10), fontsize=8, color=vermell)
A.annotate('b', (b, fb), textcoords='offset points', xytext=(5, 3), fontsize=8, color=vermell)
estils = ['-', '--', ':', '-.']
for j in range(4):
    xn = (a*fb - b*fa)/(fb - fa)
    A.plot([a, b], [fa, fb], estils[j], color=blau, linewidth=1.1)
    A.plot(xn, 0, 's', color=blau, markersize=4)
    A.plot([xn, xn], [0, f(xn)], ':', color=gris, linewidth=0.8)
    A.plot(xn, f(xn), 'o', color=vermell, markersize=5)
    A.annotate(str(j + 1), (xn, f(xn)), textcoords='offset points', xytext=(5, -10), fontsize=8, color=vermell)
    fn = f(xn)
    if fn*fb > 0:                                   # el nou punt substitueix b
        b, fb = xn, fn
        if costat == -1:
            fa = fa/2
            A.plot(a, fa, 'o', mfc='none', color=vermell, markersize=5)
        costat = -1
    else:                                           # el nou punt substitueix a
        a, fa = xn, fn
        if costat == 1:
            fb = fb/2
            A.plot(b, fb, 'o', mfc='none', color=vermell, markersize=5)
        costat = 1
A.set_title('(c) Illinois', fontsize=10)

# (d) Newton: tangents
A = ax[3]; A.set_ylim(-6, 28)
x = b0
for j in range(3):
    xn = x - f(x)/df(x)
    A.plot([x, xn], [f(x), 0], '-', color=blau, linewidth=0.9)
    A.plot([xn, xn], [0, f(xn)], ':', color=blau, linewidth=0.9)
    A.plot(x, f(x), 'o', color=vermell, markersize=5)
    A.annotate(str(j), (x, f(x)), textcoords='offset points', xytext=(5, 3), fontsize=8, color=vermell)
    x = xn
A.set_title('(d) Newton', fontsize=10)

# (e) secant: es parteix de dos punts qualssevol (no cal canvi de signe).
# Cada recta passa pels dos darrers punts; el seu tall amb l'eix dona el punt següent.
A = ax[4]; A.set_ylim(-6, 28)
x0, x1 = 3.0, 2.5
P = [x0, x1]
estils = ['-', '--', ':']
for j in range(3):
    xa, xb = P[-2], P[-1]
    xn = xb - f(xb)*(xb - xa)/(f(xb) - f(xa))
    xs = np.array([xn, max(xa, xb)]) if xn < min(xa, xb) else np.array([min(xa, xb), xn])
    xs = np.array([min(xa, xb, xn), max(xa, xb, xn)])
    pend = (f(xb) - f(xa))/(xb - xa)
    A.plot(xs, f(xb) + pend*(xs - xb), estils[j], color=blau, linewidth=1.1)
    A.plot(xn, 0, 's', color=blau, markersize=4)
    A.plot([xn, xn], [0, f(xn)], ':', color=gris, linewidth=0.8)
    P.append(xn)
for j, xj in enumerate(P):
    A.plot(xj, f(xj), 'o', color=vermell, markersize=5)
    A.annotate(str(j), (xj, f(xj)), textcoords='offset points', xytext=(5, -10) if j > 1 else (5, 3), fontsize=8, color=vermell)
A.set_title('(e) Secant', fontsize=10)

# (f) Brent: interpolació inversa quadràtica enfront del punt mig
A = ax[5]; A.set_ylim(-6, 28)
xs3 = np.array([a0, 2.0, b0]); ys3 = f(xs3)
coef = np.polyfit(ys3, xs3, 2)                           # x com a funció quadràtica de y
yy = np.linspace(-5, 24, 200)
A.plot(np.polyval(coef, yy), yy, '--', color=blau, linewidth=1.0)
A.plot(xs3, ys3, 'o', color=vermell, markersize=5)
x_iqi = np.polyval(coef, 0.0)
A.plot(x_iqi, 0, 's', color=blau, markersize=6)
A.annotate('interpolació', (x_iqi, 0), textcoords='offset points', xytext=(-8, 14), fontsize=8, color=blau)
A.plot((a0 + b0)/2, 0, 'd', color=vermell, markersize=6)
A.annotate('punt mig', ((a0 + b0)/2, 0), textcoords='offset points', xytext=(2, -14), fontsize=8, color=vermell)
A.set_title('(f) Brent', fontsize=10)

fig.tight_layout()
fig.savefig('DiAM_07_biga_apendix_Fig5.pdf', bbox_inches='tight')
fig.savefig('DiAM_07_biga_apendix_Fig5.svg', bbox_inches='tight')
plt.show()


#%% Convergència dels mètodes en l'equació de la rigidesa de la biga
L, b, E, k_req, h_min, h_max = 200.0, 20.0, 200e3, 125.0, 1.0, 30.0
h_an = 10.0
g  = lambda h: E*b*h**3/(4*L**3) - k_req
dg = lambda h: 3*E*b*h**2/(4*L**3)
nmax = 30
err = lambda h: max(abs(h - h_an), 1e-14)

def millor(punts):                                       # millor aproximació fins al moment
    e = [err(h) for h in punts]
    return np.minimum.accumulate(e)

# bisecció: g(h_min), g(h_max) i punts mitjos
lo, hi = h_min, h_max
punts_bis = [lo, hi]
for _ in range(nmax):
    m = (lo + hi)/2; punts_bis.append(m)
    if g(lo)*g(m) > 0: lo = m
    else: hi = m
# regula falsi amb [h_min, h_max]
lo, hi = h_min, h_max
punts_rf = [lo, hi]
for _ in range(nmax):
    xn = hi - g(hi)*(hi - lo)/(g(hi) - g(lo)); punts_rf.append(xn)
    if g(xn)*g(lo) > 0: lo = xn
    else: hi = xn
# Illinois amb [h_min, h_max]
lo, hi = h_min, h_max
glo, ghi = g(lo), g(hi)
costat = 0
punts_ill = [lo, hi]
for _ in range(nmax):
    xn = (lo*ghi - hi*glo)/(ghi - glo); gn = g(xn); punts_ill.append(xn)
    if gn*ghi > 0:
        hi, ghi = xn, gn
        if costat == -1: glo = glo/2
        costat = -1
    else:
        lo, glo = xn, gn
        if costat == 1: ghi = ghi/2
        costat = 1
# Newton des de h_max
h = h_max; punts_new = [h]
for _ in range(nmax):
    h = h - g(h)/dg(h); punts_new.append(h)
# secant amb h_min i h_max
hp, h = h_min, h_max; punts_sec = [hp, h]
for _ in range(nmax):
    if g(h) == g(hp): break
    hn = h - g(h)*(h - hp)/(g(h) - g(hp)); punts_sec.append(hn); hp, h = h, hn
# Brent (brentq): punts avaluats
punts_bre = []
def g_reg(h):
    punts_bre.append(h); return g(h)
brentq(g_reg, h_min, h_max, xtol=1e-12)

plt.figure(figsize=(6.2, 4.2))
for punts, etiq, estil, color in [(punts_bis, 'Bisecció', '-', 'tab:orange'),
                                  (punts_rf, 'Regula falsi', '-', 'tab:purple'),
                                  (punts_ill, 'Illinois', '-', 'tab:brown'),
                                  (punts_new, 'Newton', '-', 'tab:red'),
                                  (punts_sec, 'Secant', '-', 'tab:green'),
                                  (punts_bre, 'Brent', '-', 'tab:blue')]:
    e = millor(punts)[:nmax]
    plt.semilogy(np.arange(1, len(e) + 1), e, estil, marker='o', markersize=3, color=color, label=etiq)
plt.xlim(0, nmax); plt.ylim(1e-15, 1e2)
plt.xlabel('Avaluacions de $g$ [-]'); plt.ylabel('Error de la millor aproximació [mm]')
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=9)
plt.grid(True, alpha=0.4)
plt.savefig('DiAM_07_biga_apendix_Fig6.pdf', bbox_inches='tight')
plt.show()
