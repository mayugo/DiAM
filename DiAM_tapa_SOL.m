%% Disseny paramètric d'un suport de tapa de dipòsit
% Solució proposada. J.A. Mayugo (revisió 2026)
%
% Falla estàtica d'un material fràgil (fosa gris EN-GJL-300).
%   Apartat 2): coeficient de seguretat de les seccions A i B del disseny inicial
%   Apartat 3): redisseny per a diferents distàncies D (mètode de la bisecció)

clear; close all;

%% Dades del problema
% Material [MPa]
S_ut = 300;     % resistència a tracció
S_uc = 780;     % resistència a compressió
% Càrrega [N] i excentricitat [mm]
F  = 10e3;      % força aplicada al centre del dipòsit
de = 20;        % descentrament lateral màxim de la força
% Geometria [mm]
D0   = 180;     % distància de l'eix del dipòsit al suport (disseny inicial)
r_o0 = 120;     % radi exterior de la corba (variable de disseny)
a0   = 80;      % cantell de la secció en T (variable de disseny)
b    = 70;      % amplada de l'ala de la T (constant)
t    = 8;       % gruix de l'ala i de l'ànima (constant)
% Redisseny (apartat 3)
SF_d  = 2.2;                    % coeficient de seguretat mínim exigit
D3    = [240, 300, 360, 420];   % mm, noves distàncies D
pas   = 5;                      % mm, increment de les variables de disseny
k_min = -10;                    % increment mínim (k_min*pas = -50 mm)
k_max = 16;                     % increment màxim (k_max*pas = 80 mm)

%% Apartat 1)
% El criteri emprat a l'apèndix és el de Coulomb-Mohr modificat. És adequat
% per a una fosa gris, que és un material fràgil amb S_uc > S_ut.

%% Apartat 2): disseny inicial
[SF0, R0] = SF_suport(D0, r_o0, a0, b, t, de, F, S_ut, S_uc);

fprintf('Apartat 2): D = %.0f mm, r_o = %.0f mm, a = %.0f mm\n', D0, r_o0, a0);
fprintf('  A = %.0f mm^2, x_c = %.2f mm, r_n = %.2f mm, e = %.2f mm\n', ...
    R0.A, R0.x_c, R0.r_n, R0.e);
fprintf('  I = %.3g mm^4, I'' = %.3g mm^4, I_t = %.0f mm^4, M_e = %.0f mm^3\n', ...
    R0.I, R0.I_p, R0.I_t, R0.M_e);
fprintf('  Secció A, punt a: sigma = %5.1f, tau = %5.1f, sigma_A = %6.1f, sigma_B = %6.1f MPa -> SF = %.2f\n', ...
    R0.sAa, R0.tAa, R0.s1_Aa, R0.s2_Aa, R0.SF_Aa);
fprintf('  Secció A, punt b: sigma = %5.1f, tau = %5.1f, sigma_A = %6.1f, sigma_B = %6.1f MPa -> SF = %.2f\n', ...
    R0.sAb, R0.tAb, R0.s1_Ab, R0.s2_Ab, R0.SF_Ab);
fprintf('  Secció B, punt a: sigma = %5.1f MPa -> SF = %.2f\n', R0.sBa, R0.SF_Ba);
fprintf('  SF de les seccions A i B = %.2f\n\n', SF0);

%% Apartat 3): redisseny per a diferents D (bisecció)
% Les variables de disseny augmenten simultàniament k*pas: a = a0 + k*pas,
% r_o = r_o0 + k*pas (el radi interior r_i = r_o - a no canvia).
k3 = zeros(size(D3)); n3 = k3; SF3 = k3;
fprintf('Apartat 3): mètode de la bisecció\n');
for i = 1:length(D3)
    [k3(i), n3(i)] = biseccio(D3(i), r_o0, a0, b, t, de, F, S_ut, S_uc, ...
                              SF_d, pas, k_min, k_max);
    SF3(i) = SF_suport(D3(i), r_o0 + k3(i)*pas, a0 + k3(i)*pas, b, t, de, F, S_ut, S_uc);
    fprintf('  D = %.0f mm: %d iteracions, r_o = %.0f mm, a = %.0f mm, SF = %.2f\n', ...
        D3(i), n3(i), r_o0 + k3(i)*pas, a0 + k3(i)*pas, SF3(i));
end

% Comprovació amb el mètode incremental (augmentar k d'un en un)
fprintf('\nComprovació: mètode incremental\n');
for i = 1:length(D3)
    k = k_min;
    SF_trial = SF_suport(D3(i), r_o0 + k*pas, a0 + k*pas, b, t, de, F, S_ut, S_uc);
    cont = 0;
    while SF_trial < SF_d
        cont = cont + 1;
        k = k + 1;
        SF_trial = SF_suport(D3(i), r_o0 + k*pas, a0 + k*pas, b, t, de, F, S_ut, S_uc);
    end
    fprintf('  D = %.0f mm: %d iteracions, r_o = %.0f mm, a = %.0f mm, SF = %.2f\n', ...
        D3(i), cont, r_o0 + k*pas, a0 + k*pas, SF_trial);
end

%% Representació gràfica de la solució
figure
hold on, grid on, set(gca, 'fontsize', 20)
plot([D0, D3], [a0, a0 + k3*pas], 'x', 'MarkerSize', 12, 'LineWidth', 2.8);
xlabel('$D$ [mm]', 'Interpreter', 'latex', 'fontsize', 24)
ylabel('$a$ [mm]', 'Interpreter', 'latex', 'fontsize', 24)
print('Fig_DiAM_biga_T', '-dpdf');

%% Funcions emprades
% (en MatLab les funcions locals van al final de l'script; en Python, al principi)

function [SF, R] = SF_suport(D, r_o, a, b, t, de, F, S_ut, S_uc)
% Coeficient de seguretat estàtic de les seccions A i B del suport.
% Entrada: D, r_o, a, b, t, de [mm], F [N], S_ut, S_uc [MPa]
% Sortida: SF, coeficient de seguretat mínim [-]
%          R,  estructura amb resultats intermedis

% Propietats de la secció en T de la zona corba
r_i = r_o - a;                                          % radi interior
A   = b*t + (a - t)*t;                                  % àrea
x_c = (b*t*t/2 + (a - t)*t*(t + (a - t)/2))/A;          % centroide des de r_i
r_c = r_i + x_c;                                        % radi del centroide
r_n = A/(b*log((r_i + t)/r_i) + t*log(r_o/(r_i + t)));  % radi de la línia neutra
e   = r_c - r_n;
c_i = r_n - r_i;
c_o = r_o - r_n;
I   = b*t^3/12 + b*t*(x_c - t/2)^2 + t*(a - t)^3/12 + (a - t)*t*(t + (a - t)/2 - x_c)^2;
I_p = t*b^3/12 + (a - t)*t^3/12;                        % respecte l'eix de simetria
alpha = 1.1;                                            % factor per a seccions en T
[~, eta3_1] = calc_eta(b/t);                            % ala
[~, eta3_2] = calc_eta((a - t)/t);                      % ànima
I_t = alpha*(eta3_1*b*t^3 + eta3_2*(a - t)*t^3);        % rigidesa a torsió
M_e = (b - t)/2*t*(x_c - t/2);                          % moment estàtic de mitja ala
[eta1_b, ~] = calc_eta(a/t);                            % relació tau(T) punt b / punt a

% Secció A: sol·licitacions
M_A = (D - r_i)*F;      % N·mm, moment flector
T_A = de*F;             % N·mm, moment torsor
V_A = F;                % N,    esforç tallant

% Secció A, punt a: flexió (biga corba) + torsió + tallant
sAa  = M_A*c_i/(A*e*r_i);
tAa  = T_A*t/I_t + V_A*M_e/(t*I);
[s1_Aa, s2_Aa] = s_principals(sAa, 0, tAa);
SF_Aa = criteri_MC_mod(S_ut, S_uc, s1_Aa, s2_Aa);

% Secció A, punt b: flexió (biga corba) + torsió
sAb  = -M_A*c_o/(A*e*r_o);
tAb  = eta1_b*T_A*t/I_t;
[s1_Ab, s2_Ab] = s_principals(sAb, 0, tAb);
SF_Ab = criteri_MC_mod(S_ut, S_uc, s1_Ab, s2_Ab);

% Secció B: sol·licitacions
M_B  = (D + x_c)*F;     % N·mm, moment flector
Mp_B = de*F;            % N·mm, moment flector degut al descentrament
N_B  = F;               % N,    esforç axial

% Secció B, punt a: flexió (biga corba) + flexió (biga recta) + axial
sBa   = M_B*c_i/(A*e*r_i) + Mp_B*(b/2)/I_p + N_B/A;
SF_Ba = criteri_MC_mod(S_ut, S_uc, sBa, 0);

% Coeficient de seguretat mínim
SF = min([SF_Aa, SF_Ab, SF_Ba]);

R = struct('A', A, 'x_c', x_c, 'r_n', r_n, 'e', e, 'I', I, 'I_p', I_p, ...
           'I_t', I_t, 'M_e', M_e, 'sAa', sAa, 'tAa', tAa, 's1_Aa', s1_Aa, ...
           's2_Aa', s2_Aa, 'SF_Aa', SF_Aa, 'sAb', sAb, 'tAb', tAb, ...
           's1_Ab', s1_Ab, 's2_Ab', s2_Ab, 'SF_Ab', SF_Ab, 'sBa', sBa, 'SF_Ba', SF_Ba);
end

function [k_sol, cont] = biseccio(D, r_o0, a0, b, t, de, F, S_ut, S_uc, SF_d, pas, k_min, k_max)
% Mètode de la bisecció sobre l'increment enter k (múltiples de pas).
% Hipòtesi: SF(k_min) < SF_d i SF(k_max) >= SF_d.
% Retorna el k mínim que compleix SF >= SF_d i el nombre d'iteracions.
k_lo = k_min; k_hi = k_max; cont = 0;
while (k_hi - k_lo) > 1
    cont = cont + 1;
    k = floor((k_lo + k_hi)/2);                 % punt mig
    SF = SF_suport(D, r_o0 + k*pas, a0 + k*pas, b, t, de, F, S_ut, S_uc);
    if SF >= SF_d
        k_hi = k;                               % compleix: es pot reduir
    else
        k_lo = k;                               % no compleix: cal augmentar
    end
end
k_sol = k_hi;
end

function [eta1, eta3] = calc_eta(hb)
% Coeficients de torsió de seccions rectangulars en funció de la relació h/b
h_b  = [1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 100];
eta1 = [1.0, 0.858, 0.796, 0.775, 0.753, 0.743, 0.743, 0.743, 0.743, 0.743, 0.743, 0.743, 0.743];
eta3 = [0.141, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.299, 0.303, 0.307, 0.310, 0.313, 0.333];
eta1 = interp1(h_b, eta1, hb);
eta3 = interp1(h_b, eta3, hb);
end

function [sigmaA, sigmaB] = s_principals(s_x, s_y, s_xy)
% Tensions principals d'un estat pla de tensions (s_x, s_y, s_xy)
s_prin = eig([s_x s_xy; s_xy s_y]);
sigmaA = max(s_prin);
sigmaB = min(s_prin);
end

function [SF, mode] = criteri_MC_mod(S_ut, S_uc, sigmaA, sigmaB)
% Coeficient de seguretat amb el criteri de Coulomb-Mohr modificat
%   sigmaA = tensió principal major, sigmaB = tensió principal menor
if sigmaA < sigmaB                              % ordena sigmaA >= sigmaB
    [sigmaA, sigmaB] = deal(sigmaB, sigmaA);
end
if sigmaA >= 0 && sigmaB >= 0
    SF = S_ut/sigmaA;               mode = 'traccio';
elseif sigmaA <= 0 && sigmaB <= 0
    SF = S_uc/abs(sigmaB);          mode = 'compressio';
elseif abs(sigmaB) <= abs(sigmaA)
    SF = S_ut/sigmaA;               mode = 'traccio';
else
    SF = (S_ut*S_uc)/(S_uc*sigmaA - S_ut*(sigmaA + sigmaB));
    mode = 'traccio-compressio';
end
end
