%% Suport del contenidor: analisi a fatiga (apartat 2)
% Autor: J.A. Mayugo, 2026
clear, close all

%% Dades del problema
m_max = 1000;  m_min = 250;     % kg, massa del contenidor ple i buit
d_e = 100;  d_i = 80;           % mm, diametres exterior i interior de l'eix
D = 120;                        % mm, diametre del reforc central
r = 5;                          % mm, radi d'acord
l1 = 1100;  l2 = 500;           % mm, distancies contenidor i cilindre a l'eix de rotacio
alpha = 30*pi/180;             % angle del piston respecte la vertical
l_AB = 1200;  l_CD = 120;       % mm, longitud de l'eix i del reforc central
g = 9.81;                       % m/s2
Sut = 500;  Sy = 370;           % MPa, Acer 1
Sep = 0.5*Sut;                  % MPa, limit de fatiga estimat

%% Forces internes i propietats de la seccio
R  = m_min/m_max;                                 % reversio de la carrega
Fp = m_max*g*l1/(l2*cos(alpha));                  % N
V  = Fp/2;                                        % N
L  = (l_AB - l_CD)/2 - r;                         % mm
M  = L*Fp/2;                                      % N mm
T  = l1*m_max/2*g;                                % N mm
I  = pi*(d_e^4 - d_i^4)/64;                       % mm4
J  = pi*(d_e^4 - d_i^4)/32;                       % mm4
fprintf('Forces internes i seccio\n');
fprintf('  Fp = %.1f N, V = %.1f N, M = %.4g N mm, T = %.4g N mm, I = %.4g mm4, J = %.4g mm4\n', Fp, V, M, T, I, J);

%% Tensions nominals i coeficient de seguretat estatic (von Mises)
sig_max = M*(d_e/2)/I;  tau_max = T*(d_e/2)/J;    % MPa
sig_min = R*sig_max;    tau_min = R*tau_max;
s_p = tensions_principals(sig_max, tau_max);
sig_eqv_max = von_mises(s_p);
SF_e = Sy/sig_eqv_max;
fprintf('\nFalla estatica\n');
fprintf('  sigma_M,max = %.1f MPa, tau_T,max = %.1f MPa, sigma_eqv,max = %.1f MPa\n', sig_max, tau_max, sig_eqv_max);
fprintf('  SF_e = %.2f\n', SF_e);

%% Concentracio de tensions (nomes a les components alternes)
K_t  = K_t_rodo_flexio(D/d_e, r/d_e);
K_ts = K_t_rodo_torsio(D/d_e, r/d_e);
q    = q_sensibilitat_entalla(r, Sut, false);
q_s  = q_sensibilitat_entalla(r, Sut, true);
K_f  = 1 + q*(K_t - 1);
K_fs = 1 + q_s*(K_ts - 1);

%% Components mitjana i alterna, i criteri de Sines
sig_m_M = (sig_max + sig_min)/2;           tau_m_T = (tau_max + tau_min)/2;
sig_a_M = K_f*(sig_max - sig_min)/2;       tau_a_T = K_fs*(tau_max - tau_min)/2;
Im = tensions_principals(sig_m_M, tau_m_T);       % [I II III] mitjanes
Ia = tensions_principals(sig_a_M, tau_a_T);       % [I II III] alternes
sigma_m = sum(Im);                                % Sines: suma de les tensions principals mitjanes
sigma_a = von_mises(Ia);                          % Sines: von Mises de les components alternes
fprintf('\nEstat de tensio ciclica\n');
fprintf('  q = %.3f, K_t = %.2f, K_f = %.2f, q_s = %.3f, K_ts = %.2f, K_fs = %.2f\n', q, K_t, K_f, q_s, K_ts, K_fs);
fprintf('  (sigma_M)m = %.1f, (tau_T)m = %.1f, (sigma_M)a = %.1f, (tau_T)a = %.1f MPa\n', sig_m_M, tau_m_T, sig_a_M, tau_a_T);
fprintf('  sigma_I,m = %.1f, sigma_II,m = %.1f, sigma_III,m = %.1f MPa\n', Im);
fprintf('  sigma_I,a = %.1f, sigma_II,a = %.1f, sigma_III,a = %.1f MPa\n', Ia);
fprintf('  sigma_m (Sines) = %.1f MPa, sigma_a (von Mises) = %.1f MPa\n', sigma_m, sigma_a);

%% Coeficient de seguretat a fatiga (Goodman modificat)
[Se, ka, kb, kc] = Se_modificat(Sut, Sep, d_e);
SF_f = 1/(sigma_m/Sut + sigma_a/Se);
fprintf('\nFalla a fatiga\n');
fprintf('  ka = %.2f, kb = %.2f, kc = %.2f, Se'' = %.1f MPa, Se = %.1f MPa\n', ka, kb, kc, Sep, Se);
fprintf('  SF_f = %.2f\n', SF_f);

%% Diagrama de Haigh
% Linia de carrega: de l'origen (O) al punt de treball (A), prolongada fins a la recta de Goodman (B).
% SF_f = OB/OA.
OR = [232 113 10]/255;     % taronja
gris = [0.35 0.35 0.35];
figure('Units', 'inches', 'Position', [1 1 7.5 3.8]); hold on
h_g = plot([0 Sut], [Se 0], '-', 'Color', gris, 'LineWidth', 1.8);                  % Goodman modificat
h_c = plot([0 SF_f*sigma_m], [0 SF_f*sigma_a], '-', 'Color', OR, 'LineWidth', 1.1);   % linia de carrega O-B
plot(sigma_m, sigma_a, 'o', 'Color', OR, 'MarkerFaceColor', OR, 'MarkerSize', 6);   % A
plot(SF_f*sigma_m, SF_f*sigma_a, 'o', 'Color', OR, 'MarkerSize', 7);                % B
text(sigma_m - 22, sigma_a + 12, 'A', 'Color', OR, 'FontSize', 12);
text(SF_f*sigma_m + 12, SF_f*sigma_a + 8, 'B', 'Color', OR, 'FontSize', 12);
text(SF_f*sigma_m + 25, SF_f*sigma_a + 30, sprintf('SF_f = OB/OA = %.2f', SF_f), ...
     'Color', OR, 'FontSize', 12, 'BackgroundColor', 'w', 'Margin', 1);
text(8, Se + 12, 'S_e', 'Color', gris, 'FontSize', 14);
text(Sut - 22, 18, 'S_{ut}', 'Color', gris, 'FontSize', 14);
ax = gca;
axis equal
xlim([0, 1.08*Sut]);  ylim([0, 1.5*Se]);
ax.XTick = 0:100:500;  ax.YTick = 0:100:300;
ax.XAxis.MinorTickValues = 0:50:550;  ax.YAxis.MinorTickValues = 0:50:300;
ax.XMinorTick = 'on';  ax.YMinorTick = 'on';
ax.XGrid = 'on';  ax.YGrid = 'on';  ax.XMinorGrid = 'on';  ax.YMinorGrid = 'on';
ax.GridAlpha = 0.45;  ax.MinorGridAlpha = 0.2;  ax.Box = 'on';
ax.FontSize = 12;
xlabel('\sigma_m (Sines) [MPa]');  ylabel('\sigma_a [MPa]');
legend([h_g h_c], {'Goodman modificat', 'Linia de carrega'}, 'Location', 'northeast', 'FontSize', 12);
exportgraphics(gcf, 'DiAM_05c_contenidor_Haigh_Fig1.pdf');

%% Funcions
function K_t = K_t_rodo_flexio(D_d, r_d)
% K_t d'un rodo a flexio amb canvi de diametre (D_d = D/d, r_d = r/d)
K_t = 0.632 + 0.377*D_d^(-4.4) + r_d^(-0.5)*...
      sqrt((-0.14 - 0.363*D_d^2 + 0.503*D_d^4)/(1 - 2.39*D_d^2 + 3.368*D_d^4));
end

function K_ts = K_t_rodo_torsio(D_d, r_d)
% K_ts d'un rodo a torsio amb canvi de diametre (D_d = D/d, r_d = r/d)
K_ts = 0.78 + 0.2*D_d^(-10) + r_d^(-0.46)*...
       sqrt((0.002 - 0.125*D_d^2 + 0.123*D_d^4)/(1 - 2.75*D_d^2 + 2.55*D_d^4));
end

function q = q_sensibilitat_entalla(r, Sut, cisalla)
% Sensibilitat a l'entalla d'un acer (r en mm, Sut en MPa): q (traccio/flexio) o q_s (tallant)
S = Sut/6.9;
r_in = r/25.4;
if cisalla
    sqr_a = 0.190 - 2.51e-3*S + 1.35e-5*S^2 - 2.67e-8*S^3;
else
    sqr_a = 0.246 - 3.08e-3*S + 1.51e-5*S^2 - 2.67e-8*S^3;
end
q = 1/(1 + sqr_a/sqrt(r_in));
end

function s = tensions_principals(sx, tau)
% Tensions principals [I II III] (I >= II >= III) d'un estat pla (sx, 0, tau)
c = sx/2;  rad = hypot(sx/2, tau);
A = c + rad;  B = c - rad;
if B > 0
    s = [A, B, 0];
elseif A > 0
    s = [A, 0, B];
else
    s = [0, A, B];
end
end

function s_vm = von_mises(s)
% Tensio equivalent de von Mises a partir de les tensions principals s = [I II III]
s_vm = sqrt(2)/2*sqrt((s(1) - s(2))^2 + (s(2) - s(3))^2 + (s(3) - s(1))^2);
end

function [Se, ka, kb, kc] = Se_modificat(Sut, Sep, d)
% Limit de fatiga modificat Se = Sep*ka*kb*kc (mecanitzat, flexio i torsio)
ka = 4.51*Sut^(-0.265);
if d < 2.79
    kb = (2.79/7.62)^(-0.107);
elseif d <= 51
    kb = (d/7.62)^(-0.107);
else
    kb = 1.51*d^(-0.157);
end
kc = 1;
Se = Sep*ka*kb*kc;
end
