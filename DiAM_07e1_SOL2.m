clear; clc

% ---------------------------------------------------------------
% ESPECIFICACIONS (mm, N, MPa, kg, rad/s)
% ---------------------------------------------------------------
E = 210e3;                      % modul d'elasticitat [MPa]
Sy = 820;                       % limit elastic [MPa]
Sut = 950;                      % resistencia a traccio [MPa]
zeta = 0.02;                    % esmorteiment relatiu
X_lim = 5e-3;                   % amplitud desitjada de x [m]
g = 9.81;
SF_min = 1.8;                   % factor de seguretat minim (estatic i fatiga)
tol_k = 0.05;                   % tolerancia relativa de la rigidesa
L_min = 100; L_max = 150; dL = 5;   % longitud de flexio [mm]

% configuracions: omega [rad/s], y0 [m], m [kg]
w_tab = [60 120 45 30];
y0_tab = [1.5 0.9 1.4 1.1]*1e-3;
m_tab = [1.6 1.2 1.8 2.0];

% cataleg de seccions: gruix t i amplada b [mm]
cataleg = [
    0.25 15; 0.25 20; 0.25 25; 0.3 20; 0.3 30; 0.3 40; 0.3 50; 0.4 20;
    0.5 6; 0.5 8; 0.5 10; 0.5 12; 0.5 15; 0.5 20; 0.5 25; 0.5 28;
    0.5 30; 0.5 40; 0.5 50; 0.5 60; 0.5 70; 0.5 80; 0.5 100; 0.6 40;
    0.6 50; 0.8 10; 0.8 20; 0.8 25; 0.8 30; 0.8 40; 0.8 50; 1 10;
    1 12; 1 15; 1 20; 1 25; 1 30; 1 35; 1 40; 1 50;
    1 60; 1 70; 1 80; 1 100; 1.2 20; 1.2 30; 1.5 10; 1.5 12;
    1.5 15; 1.5 20; 1.5 25; 1.5 30; 1.5 35; 1.5 40; 1.5 45; 1.5 50;
    1.5 60; 1.5 70; 1.5 80; 1.5 100; 2 10; 2 12; 2 15; 2 20;
    2 25; 2 30; 2 35; 2 40; 2 45; 2 50; 2 60; 2 70;
    2 80; 2 100; 2.5 20; 2.5 25; 2.5 30; 2.5 50; 3 15; 3 20;
    3 30; 3 40; 3 50; 3 60; 3 70; 3 80; 3 100;
];

% ---------------------------------------------------------------
% CALCULS
% ---------------------------------------------------------------
t = cataleg(:, 1);                  % columna
b = cataleg(:, 2);                  % columna
L = L_min:dL:L_max;                 % fila

for cfg = 1:4
    [k_req, F_m, F_a] = part1(w_tab(cfg), y0_tab(cfg), m_tab(cfg), zeta, X_lim, g);

    k = rigidesa(b, t, L, E);       % matriu seccions x longituds
    err = (k - k_req)/k_req;
    SFe = SF_estatic(F_m, F_a, b, t, L, Sy);
    SFf = SF_fatiga(F_m, F_a, b, t, L, Sut);
    SF = min(SFe, SFf);

    idx = find(abs(err) <= tol_k & SF >= SF_min);
    [i, j] = ind2sub(size(k), idx);
    dis = [t(i), b(i), L(j)', k(idx), 100*err(idx), SFe(idx), SFf(idx), SF(idx)];
    [~, ordre] = sort(dis(:, 8), 'descend');
    dis = dis(ordre, :);

    fprintf('\nConfiguracio %d: k_req = %.2f N/mm, F_m = %.2f N, F_a = %.2f N -> %d dissenys\n', ...
        cfg, k_req, F_m, F_a, size(dis, 1));
    fprintf('   t x b [mm]   L [mm]  k [N/mm]  err k [%%]  SF_est  SF_fat  SF\n');
    for n = 1:size(dis, 1)
        fprintf('   %g x %g   %7d  %8.2f  %9.2f  %6.2f  %6.2f  %5.2f\n', ...
            dis(n, 1), dis(n, 2), dis(n, 3), dis(n, 4), dis(n, 5), dis(n, 6), dis(n, 7), dis(n, 8));
    end
end

% ---------------------------------------------------------------
% PART 1: rigidesa necessaria i forces (resum del SOL1)
% ---------------------------------------------------------------
function T = transmitancia(Omega, zeta)
T = sqrt(1 + (2*zeta*Omega).^2)./sqrt((1 - Omega.^2).^2 + (2*zeta*Omega).^2);
end

function [k_req, F_m, F_a] = part1(omega, y0, m, zeta, X_lim, g)
% k_req [N/mm] total, F_m i F_a [N] per suport
Omega = fzero(@(O) transmitancia(O, zeta) - X_lim/y0, [1 100]);
k_req = m*(omega/Omega)^2/1e3;
F_m = m*g/2;
F_a = m*omega^2*X_lim/2;
end

% ---------------------------------------------------------------
% PART 2: formulacio del disseny (b, t columnes; L fila)
% ---------------------------------------------------------------
function k = rigidesa(b, t, L, E)
% Rigidesa total dels dos suports [N/mm] (biga en voladiu)
I = b.*t.^3/12;
k = 2*3*E*I./L.^3;
end

function sig = tensio(F, b, t, L)
% Tensio maxima de flexio a l'encastament [MPa]
sig = 6*F*L./(b.*t.^2);
end

function SF = SF_estatic(F_m, F_a, b, t, L, Sy)
SF = Sy./tensio(F_m + F_a, b, t, L);
end

function Se = limit_fatiga(b, t, Sut)
% Limit d'enduriment corregit Se [MPa]
ka = 1.58*Sut^(-0.085);             % superficie rectificada
de = 0.808*sqrt(b.*t);              % diametre equivalent
kb = ones(size(de));
mig = de >= 2.79 & de <= 51;
gran = de > 51;
kb(mig) = (de(mig)/7.62).^(-0.107);
kb(gran) = 1.51*de(gran).^(-0.157);
kc = 1;                             % flexio
Se = 0.5*Sut*ka*kb*kc;
end

function SF = SF_fatiga(F_m, F_a, b, t, L, Sut)
% Goodman modificada: sigma_m > 0 (el pes dona F_m > 0)
sig_a = tensio(F_a, b, t, L);
sig_m = tensio(F_m, b, t, L);
Se = limit_fatiga(b, t, Sut);
SF = 1./(sig_a./Se + sig_m/Sut);
end
