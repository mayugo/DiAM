# %%% Disseny de una cel·la càrrega a tracció
# %%% J.A. Mayugo 8/01/2022

import numpy as np
import matplotlib.pyplot as plt


# %% Dades del material
# Valors típics (no de certificat). Sy: límit elàstic 0,2 %
materials = {
    # 'inox_austenitic': {'nom': 'Inox austenític (304/316, recuit)',
    #                     'E': 190e3, 'nu': 0.29, 'Sy': 240},
    '17-4PH_H1150':    {'nom': '17-4 PH H1150',
                        'E': 197e3, 'nu': 0.27, 'Sy': 724},
    '17-4PH_H1025':    {'nom': '17-4 PH H1025',
                        'E': 197e3, 'nu': 0.27, 'Sy': 1000},
    '17-4PH_H900':     {'nom': '17-4 PH H900',
                        'E': 197e3, 'nu': 0.27, 'Sy': 1170},
}

materials_a_provar = list(materials)     # tots; o p.ex. ['17-4PH_H1025']


# %% Dades problema
F = 10e3            # N, força nominal / rang

C_d = 2e-3          # V/V, sensibilitat del sensor desitjada

pont = 'complet'    # 'mig'     : 2 galgues longitudinals
                    # 'complet' : 2 longitudinals + 2 transversals (Poisson)

K = 2               # -  , factor de galga
Vo_Vi_sens = 2e-6   # V/V, sensibilitat de l'equip (+-)

SF_d = 3            # -  , SF estàtic mínim

inc = 0.1           # mm, arrodoniment de la geometria
verbose = False     # True per veure les iteracions de la bisecció


# %% Funcions emprades

def factor_pont(nu, pont):
    #% Relació eps_wheastone / eps_x segons la configuració del pont
    if pont == 'mig':
        return 2                    # ε1 + ε3 (2 longitudinals)
    elif pont == 'complet':
        return 2*(1 + nu)           # ε1 - ε2 + ε3 - ε4 = 2εx - 2εy
    raise ValueError(f"pont desconegut: {pont}")


def C_trial(d, F, E, nu, Sy, SF_d, Vo_Vi_sens, K, pont):
    #% Càlcul de la sensibilitat del disseny
    #% Input:    d (mm, diàmetre), F (N, força), E i nu (prop. elàstiques),
    #%           Sy i SF_d (prop. resistència), Vo_Vi_sens i K (extensiometria),
    #%           pont ('mig' o 'complet')
    #% Output:   C (V/V, sensibilitat del sensor)
    #%           sens_F_e (N, sensibilitat en força del disseny)
    #%           SF_e (factor seguretat)
    #%           factible (bool)
    #% Accepta d escalar o vector (numpy)

    A = np.pi*d**2/4                # % àrea
    sig = F/A                       # % tensió longitudinal (uniaxial)
    eps_x = sig/E                   # % longitudinal
    eps_wheastone = factor_pont(nu, pont)*eps_x
    C = K/4*eps_wheastone           # % V/V, sensibilitat del sensor
    sens_F_e = Vo_Vi_sens/C*F       # % N, força mínima detectable

    SF_e = Sy/sig                   # factor de seguretat
    factible = SF_e >= SF_d
    return C, sens_F_e, SF_e, factible


def biseccio_C(C_d, F, E, nu, Sy, SF_d, Vo_Vi_sens, K, pont,
               d_min=1, d_max=1000, tol=1e-3, max_iter=100, verbose=False):
    #% Diàmetre que dona C = C_d (bisecció; C decreix amb d)
    C, cont = 0, 0
    while np.abs(C_d - C) > tol*C_d and cont < max_iter:
        cont += 1
        d = (d_min + d_max)/2
        C, sF, SF, ok = C_trial(d, F, E, nu, Sy, SF_d, Vo_Vi_sens, K, pont)
        if C >= C_d:                # C massa gran → d massa petit
            d_min = d
        else:
            d_max = d
        if verbose:
            print(f'   # iter. {cont}: d = {d:.3f} mm | C = {C*1e3:.3f} mV/V | '
                  f'ΔF = {sF:.1f} N | SF = {SF:.2f}')
    return d, cont


def txt_factible(ok):
    return 'disseny factible' if ok else 'disseny NO factible'


# %% Càlcul per a cada material

print(f'Cel·la de tracció: F = {F/1e3:.0f} kN, C_d = {C_d*1e3:.2f} mV/V, '
      f'pont {pont}, SF_d = {SF_d}\n')

resultats = {}
for clau in materials_a_provar:
    mat = materials[clau]
    E, nu, Sy = mat['E'], mat['nu'], mat['Sy']
    args = (F, E, nu, Sy, SF_d, Vo_Vi_sens, K, pont)

    # Deformació necessària per a C_d
    eps_req = 4/K*C_d/factor_pont(nu, pont)

    # (1) Diàmetre que dona C_d, arrodonit a 'inc' per sota (C ≥ C_d)
    d_C, cont = biseccio_C(C_d, *args, verbose=verbose)
    d_C = np.floor(d_C/inc)*inc
    C, sF, SF, ok = C_trial(d_C, *args)

    # (2) Límit factible: diàmetre mínim que compleix SF_d
    d_SF = np.sqrt(4*F*SF_d/(np.pi*Sy))
    C_max, sF_max, _, _ = C_trial(d_SF, *args)

    resultats[clau] = dict(nom=mat['nom'], E=E, nu=nu, Sy=Sy, eps_req=eps_req,
                           d_C=d_C, C=C, sF=sF, SF=SF, ok=ok, cont=cont,
                           d_SF=d_SF, C_max=C_max, sF_max=sF_max)

    print(f"── {mat['nom']}  (E = {E/1e3:.0f} GPa, ν = {nu}, Sy = {Sy} MPa)")
    print(f'   Deformació necessària per a C_d: {eps_req*1e6:.0f} µε '
          f'(admissible amb SF_d: {Sy/E/SF_d*1e6:.0f} µε)')
    print(f'   Disseny per C_d  : d = {d_C:.1f} mm | C = {C*1e3:.3f} mV/V | '
          f'ΔF = {sF:.1f} N | {np.floor(F/sF):.0f} mesures | '
          f'SF = {SF:.2f} ({txt_factible(ok)}) [{cont} iter.]')
    if not ok:
        print('   ' + '!'*57)
        print(f'   ATENCIÓ: σ = {F/(np.pi*d_C**2/4):.0f} MPa > Sy/SF_d = '
              f'{Sy/SF_d:.0f} MPa')
        print('   ' + '!'*57)
    print(f'   Límit factible   : d_min = {d_SF:.2f} mm | C_max = '
          f'{C_max*1e3:.3f} mV/V | ΔF = {sF_max:.1f} N | '
          f'{np.floor(F/sF_max):.0f} mesures\n')


# %% Taula resum

print(f"{'Material':<34}{'d_C (mm)':>9}{'SF':>7}{'':>4}"
      f"{'d_SF (mm)':>10}{'C_max':>8}{'ΔF_max':>8}")
print('-'*80)
for r in resultats.values():
    print(f"{r['nom']:<34}{r['d_C']:>9.1f}{r['SF']:>7.2f}"
          f"{'✓' if r['ok'] else '✗':>4}{r['d_SF']:>10.2f}"
          f"{r['C_max']*1e3:>8.3f}{r['sF_max']:>8.1f}")
print('(C_max en mV/V, ΔF_max en N)')


# %% Gràfic 1: espai de disseny en funció del diàmetre (tots els materials)
#    Traç continu = zona factible (SF ≥ SF_d); traç discontinu = NO factible

colors = plt.cm.viridis(np.linspace(0.05, 0.85, len(resultats)))
d_v = np.linspace(4, 8, 500)                             # mm

fig, axs = plt.subplots(3, 1, figsize=(9, 11), sharex=True)

for (clau, r), col in zip(resultats.items(), colors):
    C_v, sF_v, SF_v, ok_v = C_trial(d_v, F, r['E'], r['nu'], r['Sy'], SF_d,
                                    Vo_Vi_sens, K, pont)
    for ax, y in zip(axs, (C_v*1e3, sF_v, SF_v)):
        ax.plot(d_v, np.where(ok_v, y, np.nan), color=col, lw=2.2)
        ax.plot(d_v, np.where(~ok_v, y, np.nan), color=col, lw=1.2, ls='--')
        ax.axvline(r['d_SF'], color=col, lw=0.8, ls=':')

    # Punt de disseny per C_d (cercle) i límit factible (quadrat)
    mk = dict(color=col, ms=8, mec='k', mew=0.6)
    axs[0].plot(r['d_C'], r['C']*1e3, 'o', **mk)
    axs[0].plot(r['d_SF'], r['C_max']*1e3, 's', **mk,
                label=f"{r['nom']}  (d_SF = {r['d_SF']:.1f} mm)")
    axs[1].plot(r['d_C'], r['sF'], 'o', **mk)
    axs[1].plot(r['d_SF'], r['sF_max'], 's', **mk)
    axs[2].plot(r['d_C'], r['SF'], 'o', **mk)
    axs[2].plot(r['d_SF'], SF_d, 's', **mk)

axs[0].axhline(C_d*1e3, color='gray', ls='--', lw=1)
axs[0].text(d_v[-1], C_d*1e3, '$C_d$ ', va='bottom', ha='right', color='gray')
axs[0].set_ylabel('C (mV/V)')
axs[0].set_title(f'Espai de disseny — cel·la de tracció {F/1e3:.0f} kN, '
                 f'pont {pont}\n● disseny per $C_d$    ■ límit factible '
                 f'($SF = SF_d$)    — factible   - - NO factible', fontsize=10)
axs[0].legend(loc='upper right', fontsize=8.5)

axs[1].set_ylabel(r'$\Delta F$ (N)')

axs[2].axhline(SF_d, color='gray', ls='--', lw=1)
axs[2].text(d_v[-1], SF_d, '$SF_d$ ', va='bottom', ha='right', color='gray')
axs[2].axhline(1, color='#c62828', lw=0.8)
axs[2].text(d_v[-1], 1, 'plastificació (SF=1) ', va='bottom', ha='right',
            color='#c62828', fontsize=9)
axs[2].set_ylabel('SF (-)')
axs[2].set_xlabel('Diàmetre d (mm)')

for ax in axs:
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('espai_disseny_materials.png', dpi=150)


# %% Gràfic 2: comparativa de materials

noms = [r['nom'] for r in resultats.values()]
x = np.arange(len(noms))

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))

# (a) SF obtingut dissenyant per C_d
SF_vals = [r['SF'] for r in resultats.values()]
bars = ax1.bar(x, SF_vals, color=['#2e7d32' if r['ok'] else '#c62828'
                                  for r in resultats.values()])
ax1.axhline(SF_d, color='k', ls='--', lw=1)
ax1.text(-0.55, SF_d, ' $SF_d$', va='bottom', ha='left')
ax1.axhline(1, color='#c62828', lw=0.8)
ax1.bar_label(bars, fmt='%.2f', padding=2)
ax1.set_xticks(x, noms, rotation=15, ha='right', fontsize=9)
ax1.set_ylabel('SF (-)')
ax1.set_title(f'SF dissenyant per $C_d$ = {C_d*1e3:.1f} mV/V')

# (b) C màxima que permet el SF exigit
Cmax_vals = [r['C_max']*1e3 for r in resultats.values()]
bars = ax2.bar(x, Cmax_vals, color=colors)
ax2.axhline(C_d*1e3, color='k', ls='--', lw=1)
ax2.text(-0.55, C_d*1e3, ' $C_d$', va='bottom', ha='left')
ax2.bar_label(bars, labels=[f"{c:.2f}\n({r['sF_max']:.0f} N)"
                            for c, r in zip(Cmax_vals, resultats.values())],
              padding=2, fontsize=9)
ax2.set_xticks(x, noms, rotation=15, ha='right', fontsize=9)
ax2.set_ylabel('$C_{max}$ (mV/V)')
ax2.set_ylim(0, max(Cmax_vals)*1.25)
ax2.set_title(f'Sensibilitat màxima amb $SF$ = {SF_d}')

for ax in (ax1, ax2):
    ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('comparativa_materials.png', dpi=150)
plt.show()