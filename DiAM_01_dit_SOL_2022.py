#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calcul de tensions: dit ortopédic
DiAM, 2023
"""

import math
import numpy as np
import matplotlib.pyplot as plt

def f_invariants2D(sigX, sigY, tauXY):
#f_invariants2D Aquesta funció calcula 4 invariants d'un estat de tensió
#               plana (bidimensional):
#               sig_A  : tensió principal màxima en el pla
#               sig_B  : tensió principal mínima en el pla
#               tau_max: tallant màxim
#               sig_VM : tensió equivalent de von Mises
    sigXY = np.array([[sigX, tauXY], [tauXY, sigY]]) #Estat de tensions 
    
    sig,v = np.linalg.eig(sigXY)
    sig = np.sort(sig)[::-1]

    sig_A = sig[0]  # Tensio principal A
    sig_B = sig[1]  # Tensio principal B
    tau_max= 0.5*(sig[0] -sig[1]) # Tallant maxim
    sig_VM = (((sig[0])**2)+((sig[1])**2)-(sig[0]*sig[1]))**0.5 # Tensio Von Mises
    return [sig_A,sig_B,tau_max,sig_VM] 

def f_cercle_Mohr2D(sigma_A,sigma_B,desc_punt,t):
# Representa el cercle de Mohr d'un estat de tensió plana (bidimensional)
# inputs:   sigma_A  : tensió principal màxima en el pla
#           sigma_B  : tensió principal mínima en el pla
#           desc_punt, etiqueta del punt representat
#           t, gruix en mm de la secció

    fig, ax = plt.subplots()
    theta= np.arange(0, 2*np.pi, 0.05) # [rad]

    # cercle sigmaI-sigmaII
    ax.plot((sigma_A-sigma_B)/2 * np.cos(theta) + (sigma_A+sigma_B)/2,
                       (sigma_A-sigma_B)/2 * np.sin(theta), 
                       '-b',label='cercle', linewidth=1.6)
    ax.plot([sigma_A, sigma_B],[0, 0],'ob',markersize=10)

    ax.set_title("Cercle de Mohr del punt " + desc_punt +
                 " amb gruix t="+str(round(t,1))+'mm')
    
    ax.set_aspect('equal')
    ax.grid(True, which='both')
    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.yaxis.tick_left()
    ax.invert_yaxis()
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    ax.xaxis.tick_bottom()
    
    ax.set_xlabel(r"$\sigma$ [MPa]", loc = "right")
    ax.set_ylabel(r"$\tau$ [MPa]")
    
    ax.plot((1), (0), ls="", marker=">", ms=10, color="k",
            transform=ax.get_yaxis_transform(), clip_on=False)
    ax.plot((0), (0), ls="", marker="v", ms=10, color="k",
            transform=ax.get_xaxis_transform(), clip_on=False)
    gap = ax.get_ylim()[0]/10
    plt.text(sigma_A,-gap,r'$\sigma_{A}$',fontsize=16)
    plt.text(sigma_B,-gap,r'$\sigma_{B}$',fontsize=16)

    plt.show()
    
def f_grafic_r(r_,desc_punt):
    fig = plt.figure()
    sig_A,  = plt.plot(t_, r_[:,0], '-',label=r'$\sigma_{A}$', linewidth=2.4)
    sig_B,  = plt.plot(t_, r_[:,1], '-',label=r'$\sigma_{B}$', linewidth=2.4)
    tau_max,= plt.plot(t_, r_[:,2], '-',label=r'$\tau_{max}$', linewidth=2.4)
    sig_eqv,= plt.plot(t_, r_[:,3], '-',label=r'$\sigma_{VM}$', linewidth=2.4)
    plt.xlabel(r'$t$ [mm]',fontsize=14)
    plt.ylabel(r'$\sigma$ [MPa]',fontsize=14)
    plt.title("Tensions principals del punt "+desc_punt+" en funció del gruix")
    plt.legend()
    plt.show()

# Dades problema
l = 53     # [mm] distància entre la forca (F2) i la secció més propera a B
b = 27     # [mm] distància entre la forca (F1) i la secció més propera a B
D =  8     # [mm] diàmetre exterior de la secció a analitzar
t_= np.arange(1, 3+0.1, 0.1) # [mm] gruix perfil circular de 1 a 3 mm

F1 = 190   # [N]  força aplicada F1
F2 = 30    # [N]  força aplicada F2

# Càlcul de les sol.licitacions sobre la secció C-C
N =   F1   # forca axial a la seccio
V_z = F2   # forca tallant en l'eix z a la seccio
M_z = b*F1 # moment creat en l'eix z
M_y = l*F2 # moment creat en l'eix y
T_x = b*F2 # torsor creat en l'eix x


r_a = []
r_b = []
r_c = []
for t in t_:
    # Càlcul de les propietats de la secció C-C (depenen del gruix t)
    d = D-2*t              # Diametre interior
    A = (math.pi/4)*(D**2-d**2)   # Area de la seccio
    I = (math.pi/64)*(D**4-d**4)  # Moment de segon ordre de la seccio circular
    J = (math.pi/32)*(D**4-d**4)  # Moment de segon ordre polar de la seccio circular
    r= D/2                 # Distancia maxima linia neutre / centroide
    
    # Estudi dels punts crítics 'a', 'b' i 'c'
    
    # Tensions en el punt 'a'
    sig_N = N/A           # Tensio normal de la forca axial en el punt a
    sig_M = M_y*r/I       # Tensio normal del moment flector en el punt a
    tau_T = T_x*r/J       # Tensio tallant del torsor en el punt a
    
    [sig_A,sig_B,tau_max,sig_VM]  = f_invariants2D(sig_N + sig_M, 0, tau_T)
    if (t==t_[0]) | (t==t_[-1]):   # calcular per t_min i t_max
        f_cercle_Mohr2D(sig_A,sig_B,desc_punt="'a'",t=t)
    r_a.append([sig_A,sig_B,tau_max,sig_VM])

    # Tensions en el punt 'b'
    sig_N = N/A        # Tensio normal de la forca axial en el punt b
    tau_V = 2*V_z/A    # Tensio tallant de la forca tallant en el punt b
    sig_M = M_z*r/I    # Tensio normal del moment flector en el punt b
    tau_T = T_x*r/J    # Tensio tallant del torsor en el punt b
    
    [sig_A,sig_B,tau_max,sig_VM]  = f_invariants2D(sig_M+sig_N,0,tau_T+tau_V);
    if (t==t_[0]) | (t==t_[-1]):   # calcular per t_min i t_max
        f_cercle_Mohr2D(sig_A,sig_B,desc_punt="'b'",t=t)
    r_b.append([sig_A,sig_B,tau_max,sig_VM])

    # Tensions en el punt 'c'
    phi = math.atan2(M_y,M_z)      # rad, valor de l'angle quan les tensions normals son maximes
    
    sig_N = N/A                    # Tensio normal de la forca radial 
    sig_M= ((M_y**2+M_z**2)**.5)*r/I  # Tensio normal del moment flector
    tau_V= 0                       # Tensio tallant [o tambe aproximadament = (2*V_z/A)*sin(phi)]
    tau_T= T_x*r/J                 # Tensio tallant del torsor
    
    [sig_A,sig_B,tau_max,sig_VM]  = f_invariants2D(sig_M+sig_N,0,tau_T+tau_V);
    if (t==t_[0]) | (t==t_[-1]):   # calcular per t_min i t_max
        f_cercle_Mohr2D(sig_A,sig_B,desc_punt="'c'",t=t)
    r_c.append([sig_A,sig_B,tau_max,sig_VM])

r_a = np.asarray(r_a)
r_b = np.asarray(r_b)
r_c = np.asarray(r_c)

f_grafic_r(r_a,"'a'")
f_grafic_r(r_b,"'b'")
f_grafic_r(r_c,"'c'")