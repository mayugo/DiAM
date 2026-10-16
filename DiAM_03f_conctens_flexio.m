%% DiAM - Exercici 04: Semieix tubular a flexio amb forat transversal i
%  canvi de seccio.
%
%  Calcul analitic de:
%    1) les sol.licitacions a les seccions del forat i del canvi de diametre;
%    2) els factors de concentracio de tensions K_t, les tensions nominals i
%       maximes, la seccio critica i el coeficient de seguretat respecte a S_ut.
%
%  Unitats: N, mm, MPa (N/mm^2), N.mm.

clear; clc;

%% ========================================================================
%  DADES
%  ========================================================================
F = 20000.0;        % [N]   carrega transversal a l'extrem lliure
L = 210.0;          % [mm]  distancia de la cara de la placa a la forca

D   = 50.0;         % [mm]  diametre exterior del tram gran (on hi ha el forat)
d   = 40.0;         % [mm]  diametre exterior del tram petit
d_i = 20.0;         % [mm]  diametre del forat axial (interior)
a   = 10.0;         % [mm]  diametre del forat transversal
r   = 3.0;          % [mm]  radi d'acord del canvi de diametre

x_forat = 30.0;     % [mm]  posicio del forat (des de la cara de la placa)
x_canvi = 60.0;     % [mm]  posicio del canvi de diametre

S_ut = 1800.0;      % [MPa] resistencia a traccio del 100Cr6 templat

Kt_canvi_grafica = 1.78;   % K_t llegit a la grafica A-15-9 de Shigley (per comparar)


%% ========================================================================
%  1) SOL.LICITACIONS
%  ========================================================================
moment_flector = @(x) F * (L - x);   % voladis amb carrega F a x = L

V        = F;
M_forat  = moment_flector(x_forat);
M_canvi  = moment_flector(x_canvi);
M_encast = moment_flector(0.0);

% Tensio tangencial maxima (fibra neutra) a la seccio petita, ordre de magnitud
tau_max = 2.0 * V / area_tub(d, d_i);


%% ========================================================================
%  2) FORAT TRANSVERSAL
%  ========================================================================
a_D      = a / D;
dD_forat = d_i / D;
[A_forat, Kt_forat] = K_t_forat_flexio(a_D, dD_forat);

Z_net          = A_forat * modul_resistent_tub(D, d_i);  % Z_net = pi A (D^4-d^4)/(32 D)
sigma0_forat   = M_forat / Z_net;
sigmamax_forat = Kt_forat * sigma0_forat;
SF_forat       = S_ut / sigmamax_forat;


%% ========================================================================
%  2) CANVI DE DIAMETRE
%  ========================================================================
D_d = D / d;
r_d = r / d;
Kt_canvi = K_t_rodo_flexio(D_d, r_d, 1);

Z_canvi        = modul_resistent_tub(d, d_i);           % seccio de diametre petit
sigma0_canvi   = M_canvi / Z_canvi;
sigmamax_canvi = Kt_canvi * sigma0_canvi;
SF_canvi       = S_ut / sigmamax_canvi;

% Mateix calcul amb el K_t llegit a la grafica
sigmamax_canvi_graf = Kt_canvi_grafica * sigma0_canvi;
SF_canvi_graf       = S_ut / sigmamax_canvi_graf;


%% ========================================================================
%  RESULTATS
%  ========================================================================
linia = repmat('=', 1, 66);
fprintf('%s\n', linia);
fprintf('1) SOL.LICITACIONS\n');
fprintf('%s\n', linia);
fprintf('Esforc tallant V                     = %8.1f kN\n', V/1000);
fprintf('Moment a l''encastament  (x =   0 mm) = %8.1f N.m\n', M_encast/1000);
fprintf('Moment al forat         (x = %3.0f mm) = %8.1f N.m\n', x_forat, M_forat/1000);
fprintf('Moment al canvi seccio  (x = %3.0f mm) = %8.1f N.m\n', x_canvi, M_canvi/1000);
fprintf('Tensio tangencial maxima (tallant)   = %8.1f MPa\n', tau_max);

fprintf('\n%s\n', linia);
fprintf('2) CONCENTRACIO DE TENSIONS\n');
fprintf('%s\n', linia);
fprintf('Forat:  a/D = %.3f,  d/D = %.3f\n', a_D, dD_forat);
fprintf('        A   = %.3f,  K_t = %.3f\n', A_forat, Kt_forat);
fprintf('        Z_net = %.1f mm^3\n', Z_net);
fprintf('Canvi:  D/d = %.3f,  r/d = %.3f\n', D_d, r_d);
fprintf('        K_t = %.3f\n', Kt_canvi);
fprintf('        Z   = %.1f mm^3\n', Z_canvi);

fprintf('\n%16s%10s%8s%10s%12s%8s\n', '', 'M [N.m]', 'K_t', 'sigma_0', 'sigma_max', 'SF');
fprintf('%16s%10s%8s%10s%12s%8s\n', '', '', '', '[MPa]', '[MPa]', '');
fprintf('%s\n', repmat('-', 1, 64));
fprintf('%-16s%10.1f%8.3f%10.1f%12.1f%8.2f\n', 'Forat', ...
        M_forat/1000, Kt_forat, sigma0_forat, sigmamax_forat, SF_forat);
fprintf('%-16s%10.1f%8.3f%10.1f%12.1f%8.2f\n', 'Canvi seccio', ...
        M_canvi/1000, Kt_canvi, sigma0_canvi, sigmamax_canvi, SF_canvi);

if sigmamax_forat >= sigmamax_canvi
    critica = 'el forat';           smax_c = sigmamax_forat; smax_a = sigmamax_canvi;
else
    critica = 'el canvi de seccio'; smax_c = sigmamax_canvi; smax_a = sigmamax_forat;
end
dif = 100.0 * (smax_c - smax_a) / smax_a;
fprintf('\nSeccio critica: %s (sigma_max %.1f %% superior a l''altra discontinuitat)\n', ...
        critica, dif);
fprintf('Coeficient de seguretat minim: SF = %.2f\n', S_ut / smax_c);

fprintf('\n%s\n', repmat('-', 1, 66));
fprintf('Comparacio amb el K_t llegit a la grafica del canvi de seccio\n');
fprintf('%s\n', repmat('-', 1, 66));
fprintf('K_t funcio  = %.3f  ->  sigma_max = %6.1f MPa,  SF = %.2f\n', ...
        Kt_canvi, sigmamax_canvi, SF_canvi);
fprintf('K_t grafica = %.3f  ->  sigma_max = %6.1f MPa,  SF = %.2f\n', ...
        Kt_canvi_grafica, sigmamax_canvi_graf, SF_canvi_graf);
fprintf('Diferencia en K_t: %.1f %%\n', 100*(Kt_canvi/Kt_canvi_grafica - 1));
if sigmamax_forat >= sigmamax_canvi_graf
    critica_graf = 'el forat';
else
    critica_graf = 'el canvi de seccio';
end
fprintf('Amb el K_t de la grafica, la seccio critica seria %s.\n', critica_graf);


%% ========================================================================
%  FUNCIONS LOCALS
%  ========================================================================

function [A, K_t] = K_t_forat_flexio(a_D, d_D)
%K_T_FORAT_FLEXIO  Retorna A i K_t per a la flexio d'una barra rodona o un
%   tub amb un forat transversal, interpolant linealment a la taula 3.1.
%
%   a_D : relacio a/D (diametre del forat / diametre exterior)
%   d_D : relacio d/D (diametre interior / diametre exterior); 0 si es massis
%
%   Primer s'interpola cada columna (d/D fix) en funcio de a/D i despres
%   s'interpola el resultat en funcio de d/D.

    % --- Dades de la taula 3.1 ---
    t_a_D = [0.050 0.075 0.10 0.125 0.15 0.175 0.20 0.225 0.25 0.275 0.30]';
    t_d_D = [0.0 0.6 0.9];

    % Files: a/D; columnes: d/D = 0, 0.6, 0.9
    t_A = [0.88 0.91 0.92
           0.86 0.88 0.89
           0.83 0.85 0.86
           0.80 0.82 0.82
           0.76 0.79 0.79
           0.72 0.75 0.76
           0.68 0.72 0.73
           0.65 0.68 0.69
           0.61 0.64 0.67
           0.58 0.61 0.66
           0.54 0.58 0.64];
    t_Kt = [2.42 2.55 2.63
            2.35 2.43 2.55
            2.27 2.36 2.49
            2.20 2.32 2.41
            2.15 2.29 2.39
            2.10 2.26 2.38
            2.07 2.23 2.39
            2.04 2.21 2.40
            2.00 2.18 2.42
            1.97 2.16 2.48
            1.94 2.14 2.52];

    % --- Comprovacio del rang ---
    if a_D < t_a_D(1) || a_D > t_a_D(end)
        error('a/D = %.3f fora del rang de la taula 3.1', a_D);
    end
    if d_D < t_d_D(1) || d_D > t_d_D(end)
        error('d/D = %.3f fora del rang de la taula 3.1', d_D);
    end

    % --- Interpolacio lineal doble ---
    interpola = @(taula) interp1(t_d_D, interp1(t_a_D, taula, a_D), d_D);

    A   = interpola(t_A);
    K_t = interpola(t_Kt);
end


function [K_t] = K_t_rodo_flexio(D_d, r_d, i_flag)
%Calcula la K_t d'un rodo sotmes a flexio amb un canvi de seccio de D a d:
%   D = diametre gran; d= diametre petit; r=radi de l'entalla
%   D_d = D/d ; r_d = r/d
%   i_flag = 1, si es vol K_t per tensio normal
%   i_flag <> 1, si es vol K_t per von Mises

    if i_flag == 1
        K_t = 0.632+0.377*D_d^(-4.4)+r_d^(-0.5)*...
          sqrt((-0.14-0.363*D_d^2+0.503*D_d^4)/(1-2.39*D_d^2+3.368*D_d^4));
    else
        K_t = 0.622+0.38*D_d^(-4.3)+r_d^(-0.5)*...
          sqrt((-0.322-0.277*D_d^2+0.599*D_d^4)/(1-2.55*D_d^2+5.27*D_d^4));
    end
end


function Z = modul_resistent_tub(D_ext, D_int)
%Modul resistent a flexio d'un tub: Z = pi (D^4 - d^4) / (32 D)
    Z = pi * (D_ext^4 - D_int^4) / (32.0 * D_ext);
end


function A = area_tub(D_ext, D_int)
    A = pi / 4.0 * (D_ext^2 - D_int^2);
end
