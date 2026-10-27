%% Suport del contenidor: coeficients de seguretat per a diferents acers i radis (apartat 3)
% Autor: J.A. Mayugo, 2026
clear, close all

%% Dades del problema
P.m_max = 1000;  P.m_min = 250;     % kg, massa del contenidor ple i buit
P.d_e = 100;  P.d_i = 80;           % mm, diametres exterior i interior de l'eix
P.D = 120;                          % mm, diametre del reforc central
P.l1 = 1100;  P.l2 = 500;           % mm, distancies contenidor i cilindre a l'eix de rotacio
P.alpha = 30*pi/180;               % angle del piston respecte la vertical
P.l_AB = 1200;  P.l_CD = 120;       % mm, longitud de l'eix i del reforc central
P.g = 9.81;                         % m/s2

% Acers: [Sut Sy Se']; Se' = NaN indica que s'estima com 0.5*Sut
acers = [500 370 NaN;
         550 450 305;
         600 470 NaN;
         650 500 360];
radis = [3.5 5 7.5];                % mm

%% Calcul per a tots els acers i radis
n_a = size(acers, 1);  n_r = numel(radis);
SF_e = zeros(n_a, n_r);
SF_f = zeros(n_a, n_r);
for i = 1:n_a
    Sut = acers(i, 1);  Sy = acers(i, 2);  Sep = acers(i, 3);
    if isnan(Sep), Sep = 0.5*Sut; end
    for j = 1:n_r
        [SF_e(i, j), SF_f(i, j)] = coeficients_seguretat(Sut, Sy, Sep, radis(j), P);
    end
end

fprintf('Acer   r [mm]    SF_e    SF_f\n');
for i = 1:n_a
    for j = 1:n_r
        if j == 2, etq = sprintf('%d', i); else, etq = ''; end
        fprintf('%4s   %6.2f   %6.2f  %6.2f\n', etq, radis(j), SF_e(i, j), SF_f(i, j));
    end
end

%% Representacio grafica
colors = [0 114 189; 217 83 25; 237 177 32; 126 47 142]/255;
noms = arrayfun(@(x) ['Acer ' num2str(x)], 1:n_a, 'UniformOutput', false);
figure('Units', 'inches', 'Position', [1 1 11.2 4.2]);
for k = 1:2
    subplot(1, 2, k)
    if k == 1, SF = SF_e; etq = 'SF_e'; lim = [1 4]; else, SF = SF_f; etq = 'SF_f'; lim = [1 2.5]; end
    b = bar(SF', 'grouped');
    for i = 1:n_a, b(i).FaceColor = colors(i, :); end
    set(gca, 'XTickLabel', arrayfun(@(x) num2str(x), radis, 'UniformOutput', false), 'FontSize', 14);
    xlabel('radi r [mm]');  ylabel(etq);
    ylim(lim);  grid on
    if k == 2, legend(noms, 'Location', 'northwest', 'FontSize', 11); end
end
exportgraphics(gcf, 'DiAM_05c_contenidor_Fig2.pdf');

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

function Se = Se_modificat(Sut, Sep, d)
% Limit de fatiga modificat Se = Sep*ka*kb*kc (mecanitzat, flexio i torsio)
ka = 4.51*Sut^(-0.265);
if d < 2.79
    kb = (2.79/7.62)^(-0.107);
elseif d <= 51
    kb = (d/7.62)^(-0.107);
else
    kb = 1.51*d^(-0.157);
end
Se = Sep*ka*kb*1;
end

function [SF_e, SF_f] = coeficients_seguretat(Sut, Sy, Sep, r, P)
% Coeficients de seguretat estatic (von Mises) i a fatiga (Goodman modificat amb el criteri de Sines)
R = P.m_min/P.m_max;
Fp = P.m_max*P.g*P.l1/(P.l2*cos(P.alpha));
M = ((P.l_AB - P.l_CD)/2 - r)*Fp/2;
T = P.l1*P.m_max/2*P.g;
I = pi*(P.d_e^4 - P.d_i^4)/64;
J = pi*(P.d_e^4 - P.d_i^4)/32;
sig_max = M*(P.d_e/2)/I;  tau_max = T*(P.d_e/2)/J;
sig_min = R*sig_max;      tau_min = R*tau_max;
% estatic
SF_e = Sy/von_mises(tensions_principals(sig_max, tau_max));
% concentracio de tensions (nomes a les components alternes)
K_f  = 1 + q_sensibilitat_entalla(r, Sut, false)*(K_t_rodo_flexio(P.D/P.d_e, r/P.d_e) - 1);
K_fs = 1 + q_sensibilitat_entalla(r, Sut, true)*(K_t_rodo_torsio(P.D/P.d_e, r/P.d_e) - 1);
% criteri de Sines
sigma_m = sum(tensions_principals((sig_max + sig_min)/2, (tau_max + tau_min)/2));
sigma_a = von_mises(tensions_principals(K_f*(sig_max - sig_min)/2, K_fs*(tau_max - tau_min)/2));
Se = Se_modificat(Sut, Sep, P.d_e);
SF_f = 1/(sigma_m/Sut + sigma_a/Se);
end
