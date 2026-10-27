"""
Analisi a fatiga de la manovella (apartat 2)
Cicle de carrega: la forca F es desplaca des de d = 20 mm fins a d = 60 mm i torna.
Per a cada punt estudiat (seccio 1 punts a i b, seccio 2 punts a i b):
  1) tensions principals (I, II, III) a cada distancia d (instant de carrega)
  2) components mitjana i alterna de cada tensio principal al llarg del cicle
  3) tensions equivalents (von Mises) mitjana i alterna
  4) factor de seguretat a fatiga (Goodman modificada)
  5) comprovacio de la fluencia (von Mises) a l'instant mes desfavorable
J.A. Mayugo (2026)
"""
import numpy as np
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------------
# Funcions
# ----------------------------------------------------------------------------
def section_properties(b, h, t, eta_a, eta_b, alpha=1.1):
    """Inertia a flexio Ix [mm^4], constant de torsio It [mm^4] i centre de
    tallant c [mm] de la seccio en U (amplada b, cantell h, gruix t)."""
    Ix = b*h**3/12 - (b - t)*(h - 2*t)**3/12
    It = alpha*(2*eta_a*b*t**3 + eta_b*(h - 2*t)*t**3)
    c = (b - t/2)**2*(h - t)**2*t/4/Ix
    return Ix, It, c


def loads(F, d, s, c, t):
    """Tallant V [N], moment flector M [N mm] i moment torsor T [N mm].
    La distancia d pot ser un vector (un valor per cada instant del cicle)."""
    V = F
    M = F*s
    T = F*(d + t/2 - c)
    return V, M, T


def nominal_stresses(V, M, T, b, h, t, Ix, It, point):
    """Tensio normal per flexio sig i tallant total tau [MPa] al punt 'a'
    (fibra extrema) o 'b' (eix neutre)."""
    tau_T = T*t/It
    if point == 'a':
        sig = M*(h/2)/Ix*np.ones_like(tau_T)
        Me = (b - t)*t*(h - t)/2
    else:
        sig = np.zeros_like(tau_T)
        Me = b*t*(h - t)/2 + (h/2 - t)*t*(h/2 - t)/2
    tau_V = V*Me/t/Ix
    return sig, tau_V + tau_T


def principal_stresses(sig, tau):
    """Tensions principals ordenades I >= II >= III (estat pla: una val 0).
    Retorna una matriu 3 x n (una columna per cada instant del cicle)."""
    R = np.sqrt((sig/2)**2 + tau**2)
    S = np.vstack((sig/2 + R, np.zeros_like(sig), sig/2 - R))
    return -np.sort(-S, axis=0)


def cycle_components(S):
    """Components mitjana i alterna de cada tensio principal (fila de S)
    al llarg del cicle."""
    Smax = S.max(axis=1)
    Smin = S.min(axis=1)
    return (Smax + Smin)/2, (Smax - Smin)/2


def von_mises(S):
    """Tensio equivalent de von Mises a partir de les tensions principals."""
    return np.sqrt(2)/2*np.sqrt((S[0] - S[1])**2 + (S[1] - S[2])**2
                                + (S[2] - S[0])**2)


def equivalent_stresses(S_m, S_a):
    """Tensions equivalents mitjana (amb signe de la suma de tensions
    principals) i alterna."""
    signe = 1.0 if S_m.sum() >= 0 else -1.0   # si la suma es 0, signe positiu
    return signe*von_mises(S_m), von_mises(S_a)


def goodman(sig_a, sig_m, Se, Sut):
    """Factor de seguretat a fatiga, regla de Goodman modificada."""
    if sig_m > 0:
        return 1/(sig_a/Se + sig_m/Sut)
    return Se/sig_a


# ----------------------------------------------------------------------------
# Dades del problema
# ----------------------------------------------------------------------------
F = 150.0                           # N, forca aplicada
d = np.arange(20.0, 60.0 + 1, 1.0)  # mm, distancies del cicle (cada 1 mm)
Sut = 150.0                         # MPa, resistencia ultima
Sy = 70.0                           # MPa, limit elastic
Se = 45.0                           # MPa, limit de fatiga (ja corregit)

t = 4.0                             # mm, gruix (seccions 1 i 2)
sections = {                        # s: cota; b: amplada; eta: coef. torsio
    'Seccio 1': {'s': 180.0, 'b': 25.0, 'eta': (0.300, 0.306)},
    'Seccio 2': {'s': 15.0,  'b': 20.0, 'eta': (0.291, 0.306)},
}
points = ('a', 'b')

# ----------------------------------------------------------------------------
# Calcul
# ----------------------------------------------------------------------------
results = []
for name, sec in sections.items():
    h = 30 + (40 - 30)/200*sec['s']
    Ix, It, c = section_properties(sec['b'], h, t, *sec['eta'])
    for point in points:
        V, M, T = loads(F, d, sec['s'], c, t)
        sig, tau = nominal_stresses(V, M, T, sec['b'], h, t, Ix, It, point)
        S = principal_stresses(sig, tau)               # 3 x n
        S_m, S_a = cycle_components(S)                 # 3 components
        sig_m, sig_a = equivalent_stresses(S_m, S_a)
        seqv = von_mises(S)                            # a cada instant
        i_crit = np.argmax(seqv)                       # instant mes desfavorable
        results.append({
            'label': f"{name}, punt {point}", 'S': S, 'S_m': S_m, 'S_a': S_a,
            'sig_m': sig_m, 'sig_a': sig_a,
            'SFf': goodman(sig_a, sig_m, Se, Sut),
            'SFy': Sy/seqv[i_crit], 'd_crit': d[i_crit]})

# ----------------------------------------------------------------------------
# Resultats
# ----------------------------------------------------------------------------
for r in results:
    print(f"{r['label']}: sigma_m = {r['sig_m']:5.1f} MPa, "
          f"sigma_a = {r['sig_a']:5.1f} MPa, "
          f"SFf = {r['SFf']:4.2f}, SFy = {r['SFy']:4.2f} "
          f"(instant mes desfavorable: d = {r['d_crit']:.0f} mm)")
    print("    components mitjanes   (I, II, III) = "
          + ", ".join(f"{v:6.1f}" for v in r['S_m']) + " MPa")
    print("    components alternes   (I, II, III) = "
          + ", ".join(f"{v:6.1f}" for v in r['S_a']) + " MPa")

# ----------------------------------------------------------------------------
# Representacio del cicle complet (anada 20 -> 60 mm i tornada 60 -> 20 mm)
# ----------------------------------------------------------------------------
d_cicle = np.concatenate((d, d[-2::-1]))     # la tornada repeteix els valors
x = np.arange(d_cicle.size)                  # instant del cicle
x_cota = x.max() + 4                         # posicio de la cota (fora del requadre)
tol = 0.05                                   # MPa, tensio considerada nul.la
colors = ['C0', 'C1', 'C2']
roman = ['I', 'II', 'III']
plt.rcParams.update({'font.size': 13})

fig, axs = plt.subplots(4, 1, figsize=(10, 14), sharex=True)
for ax, r in zip(axs, results):
    S_cicle = np.concatenate((r['S'], r['S'][:, -2::-1]), axis=1)
    for k in range(3):
        sm, sa = r['S_m'][k], r['S_a'][k]
        ax.plot(x, S_cicle[k], '-', color=colors[k], lw=2,
                label=rf'$\sigma_{{\mathrm{{{roman[k]}}}}}$')
        ax.hlines(sm, x.min(), x.max(), color=colors[k], ls='--', lw=1.2)
        if abs(sm) < tol and abs(sa) < tol:                  # tensio nul.la
            txt = (rf'$\sigma_{{\mathrm{{{roman[k]}}},m}}=0$' + '\n'
                   rf'$\sigma_{{\mathrm{{{roman[k]}}},a}}=0$')
        else:
            ext = S_cicle[k, np.argmax(np.abs(S_cicle[k] - sm))]   # extrem
            ax.plot([x.max(), x_cota], [ext, ext], ':', color=colors[k],
                    lw=1.2, clip_on=False)
            ax.annotate('', xy=(x_cota, ext), xytext=(x_cota, sm),
                        arrowprops=dict(arrowstyle='<->', color=colors[k],
                                        lw=1.5), annotation_clip=False)
            txt = (rf'$\sigma_{{\mathrm{{{roman[k]}}},m}}={sm:.1f}$' + '\n'
                   rf'$\sigma_{{\mathrm{{{roman[k]}}},a}}={sa:.1f}$')
        ax.text(x_cota + 1.5, sm, txt, color=colors[k], fontsize=13,
                va='center', ha='left')
    ax.axhline(0, color='k', lw=0.6)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylabel('Tensio (MPa)')
    ax.set_title(r['label'], fontsize=14)
    ax.grid(True, alpha=0.3)
ticks = np.arange(0, x.max() + 1, 10)
axs[-1].set_xticks(ticks)
axs[-1].set_xticklabels([f'{v:.0f}' for v in d_cicle[ticks]])
axs[-1].set_xlabel('Instant del cicle, d (mm): anada 20 -> 60 mm i tornada 60 -> 20 mm')
fig.tight_layout(rect=(0, 0, 0.84, 1))
fig.savefig('DiAM_05b_manovella_SOL2_cicle.pdf')
fig.savefig('DiAM_05b_manovella_SOL2_cicle.svg')

# ----------------------------------------------------------------------------
# Esquerra: estatic (tensions principals maximes i limit de von Mises)
# Dreta: fatiga (diagrama de vida infinita, tipus Haigh)
# ----------------------------------------------------------------------------
fig2, (axp, axh) = plt.subplots(1, 2, figsize=(13, 6))
pcol = ['C3', 'C4', 'C5', 'C6']
marks = ['o', 's', '^', 'D']
msz = [13, 10, 13, 10]    # mides diferents perque els punts propers es vegin tots
names_p = ['1a', '1b', '2a', '2b']

# a) Estatic: tensions principals no nul.les (A = I, B = III) a l'instant mes
#    desfavorable (d = 60 mm) i limit de von Mises
#    (estat pla: sigma_A^2 - sigma_A sigma_B + sigma_B^2 = Sy^2)
sA, sB = np.meshgrid(np.linspace(-5, 90, 400), np.linspace(-90, 5, 400))
SF_lines = [1, 1.5, 2, 3]                       # linies de factor de seguretat
bbox = dict(fc='w', ec='none', pad=0.5)
for SF in SF_lines:                          # linies de SF constant (blau fi)
    axp.contour(sA, sB, sA**2 - sA*sB + sB**2 - (Sy/SF)**2, levels=[0],
                colors='tab:blue', linewidths=0.8)
    axp.text(0.5, -Sy/SF, f'SF = {SF:g}', color='tab:blue', fontsize=11,
             va='center', ha='left', bbox=bbox)
axp.contour(sA, sB, sA**2 - sA*sB + sB**2 - Sy**2, levels=[0],
            colors='darkorange', linewidths=5)    # limit de von Mises
for r, pc, mk, nm, ms in zip(results, pcol, marks, names_p, msz):
    i_c = np.argmax(von_mises(r['S']))                   # instant mes desfavorable
    axp.plot(r['S'][0, i_c], r['S'][2, i_c], mk, color=pc, ms=ms-1,
             label=f'Punt {nm}')
axp.axhline(0, color='k', lw=0.6)
axp.axvline(0, color='k', lw=0.6)
axp.set_xlim(-5, 90)
axp.set_ylim(-90, 5)
axp.set_aspect('equal')
axp.set_xlabel(r'$\sigma_A=\sigma_{\mathrm{I}}$ (MPa)')
axp.set_ylabel(r'$\sigma_B=\sigma_{\mathrm{III}}$ (MPa)')
axp.grid(True, alpha=0.3)
axp.legend(loc='lower right', fontsize=11)

# b) Fatiga: diagrama de vida infinita (Haigh) amb Goodman modificada (plataforma
#    sigma_a = Se per a sigma_m < 0)
for SF in SF_lines:                          # linies de SF constant (blau fi)
    axh.plot([-20, 0, Sut/SF], [Se/SF, Se/SF, 0], '-', color='tab:blue', lw=0.8)
    axh.text(-10, Se/SF, f'SF = {SF:g}', color='tab:blue', fontsize=11,
             ha='center', va='bottom')
axh.plot([-20, 0, Sut], [Se, Se, 0], '-', color='darkorange', lw=5,
         label='Goodman modificada')
for r, pc, mk, nm, ms in zip(results, pcol, marks, names_p, msz):
    axh.plot(r['sig_m'], r['sig_a'], mk, color=pc, ms=ms, mfc='none', mew=2,
             label=f'Punt {nm}')
axh.set_xlim(-20, 160)
axh.set_ylim(0, 80)
axh.set_xlabel(r'$\sigma_m$ (MPa)')
axh.set_ylabel(r'$\sigma_a$ (MPa)')
axh.grid(True, alpha=0.3)
axh.legend(loc='upper right', fontsize=11)
fig2.tight_layout()
fig2.savefig('DiAM_05b_manovella_SOL2_haigh.pdf')
fig2.savefig('DiAM_05b_manovella_SOL2_haigh.svg')
plt.show()
