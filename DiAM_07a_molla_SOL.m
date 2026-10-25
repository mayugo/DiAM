%% Disseny d'un ressort de compressió per a una vàlvula d'admissió
% J.A. Mayugo, 2026
%
% Problema determinat: cal trobar la rigidesa k i la longitud lliure L_0 d'una molla que
% faci exactament les forces F1 i F2 a les dues posicions de funcionament.
% La solució analítica es calcula a l'apartat 1. Aquí es busca la mateixa solució per prova i error:
%   Apartat 2): mètode de Montecarlo (cerca aleatòria)
%   Apartat 3): optimització (mínims quadrats no lineals, lsqnonlin)
% Cal la Optimization Toolbox (lsqnonlin).

clear; close all;
global N_AVAL HISTORIAL           % comptador d'avaluacions i historial (només per a lsqnonlin)

blau = [0.122 0.467 0.706];       % color 'tab:blue'
vermell = [0.839 0.153 0.157];    % color 'tab:red'
gris = [0.35 0.35 0.35];

%% Dades del problema
% Posicions de funcionament [mm]
L1     = 28.0;      % longitud de la molla amb la vàlvula tancada
y12    = 4.0;       % recorregut de la vàlvula
L2     = L1 - y12;  % longitud de la molla amb la vàlvula oberta
% Forces [N]
F1     = 20.0;      % força amb la vàlvula tancada
F2     = 40.0;      % força amb la vàlvula oberta
% Espai de cerca
k_max  = 10.0;      % N/mm, rigidesa màxima considerada
L0_max = 40.0;      % mm, longitud lliure màxima considerada
% Mètode de Montecarlo
n_MC   = [1e3, 1e4, 1e5, 1e6];  % nombre de realitzacions
rng(2030);                      % llavor fixa per reproduir els resultats

%% Apartat 1): solució analítica (per comparar)
k_an  = (F2 - F1)/y12;          % rigidesa
y1_an = F1/k_an;                % deformació amb la vàlvula tancada
y2_an = F2/k_an;                % deformació amb la vàlvula oberta
L0_an = L1 + y1_an;             % longitud lliure
fprintf('Solució analítica: k = %.3f N/mm, y1 = %.2f mm, y2 = %.2f mm, L0 = %.3f mm\n\n', ...
    k_an, y1_an, y2_an, L0_an);

%% Representació gràfica: les dues corbes de les forces i la solució del sistema
k_    = linspace(0.05, k_max, 400);
L0_F1 = L1 + F1./k_;            % corba F_1 = F1
L0_F2 = L2 + F2./k_;            % corba F_2 = F2
figure('Position', [100 100 400 360])
hold on, grid on, set(gca, 'fontsize', 12)
plot(k_, L0_F1, '--', 'Color', gris, 'LineWidth', 1.4);
plot(k_, L0_F2, '-.', 'Color', gris, 'LineWidth', 1.4);
plot(k_an, L0_an, '+', 'Color', 'b', 'MarkerSize', 14, 'LineWidth', 1.8);
xlim([0 k_max]); ylim([L1 L0_max]);
xlabel('$k$ [N/mm]', 'Interpreter', 'latex')
ylabel('$L_0$ [mm]', 'Interpreter', 'latex')
legend({'Corba $F_{(1)} = F_1$', 'Corba $F_{(2)} = F_2$', 'Solució del sistema'}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 9)
print('DiAM_07a_molla_analitica_Fig0', '-dpdf');

%% Apartat 2): mètode de Montecarlo
fprintf('Apartat 2): mètode de Montecarlo\n');
fprintf('%9s %9s %9s %8s %8s %8s %9s %11s\n', 'n', 'k [N/mm]', 'L0 [mm]', 'F_1 [N]', ...
    'F_2 [N]', 'error k', 'error L0', 'temps [ms]');
montecarlo(1000, k_max, L0_max, L1, L2, F1, F2);        % escalfament (no es compta)
for i = 1:length(n_MC)
    n = n_MC(i);
    tic;
    R = montecarlo(n, k_max, L0_max, L1, L2, F1, F2);
    dt = toc*1e3;                                       % ms
    [F_1, F_2] = molla(R.k_opt, R.L0_opt, L1, L2);
    err_k  = (R.k_opt - k_an)/k_an*100;
    err_L0 = (R.L0_opt - L0_an)/L0_an*100;
    fprintf('%9d %9.3f %9.3f %8.3f %8.3f %7.2f%% %8.3f%% %11.1f\n', ...
        n, R.k_opt, R.L0_opt, F_1, F_2, err_k, err_L0, dt);
end
R_vis = montecarlo(10000, k_max, L0_max, L1, L2, F1, F2);   % per al gràfic

%% Representació gràfica: error de cada disseny aleatori
figure('Position', [100 100 600 400])
hold on, grid on, set(gca, 'fontsize', 12)
scatter(R_vis.k, R_vis.L0, 4, R_vis.f, 'filled');
try, colormap(flipud(viridis)); catch, colormap(flipud(parula)); end
set(gca, 'ColorScale', 'log'); caxis([1e-3 10]);
cb = colorbar; ylabel(cb, 'Error $f$ [-]', 'Interpreter', 'latex')
h1 = plot(k_, L0_F1, '--', 'Color', 'k', 'LineWidth', 1.4);
h2 = plot(k_, L0_F2, '-.', 'Color', 'k', 'LineWidth', 1.4);
h3 = plot(R_vis.k_opt, R_vis.L0_opt, 'x', 'Color', vermell, 'MarkerSize', 10, 'LineWidth', 2.5);
h4 = plot(k_an, L0_an, '+', 'Color', 'b', 'MarkerSize', 14, 'LineWidth', 1.8);
xlim([0 k_max]); ylim([L1 L0_max]);
xlabel('$k$ [N/mm]', 'Interpreter', 'latex')
ylabel('$L_0$ [mm]', 'Interpreter', 'latex')
legend([h1 h2 h3 h4], {'Corba $F_{(1)} = F_1$', 'Corba $F_{(2)} = F_2$', 'Montecarlo', ...
    'Solució analítica'}, 'Interpreter', 'latex', 'Location', 'northeastoutside', 'FontSize', 9)
print('DiAM_07a_molla_MC_Fig1', '-dpdf');

%% Apartat 3): optimització per mínims quadrats no lineals
x0 = [k_max, L0_max];                                   % punt inicial (no cal que sigui factible)
optimitzacio(x0, k_max, L0_max, L1, L2, F1, F2);        % escalfament (no es compta)
R_opt = optimitzacio(x0, k_max, L0_max, L1, L2, F1, F2);
dt_opt = R_opt.t;                                       % ms (només la resolució)
err_k_opt  = (R_opt.k_opt - k_an)/k_an*100;
err_L0_opt = (R_opt.L0_opt - L0_an)/L0_an*100;
[F_1o, F_2o] = molla(R_opt.k_opt, R_opt.L0_opt, L1, L2);
fprintf('\nApartat 3): optimització per mínims quadrats\n');
fprintf('  punt inicial: k = %.3f N/mm, L0 = %.1f mm\n', x0(1), x0(2));
fprintf('  solució: k = %.6f N/mm, L0 = %.6f mm (errors %.1e%% i %.1e%%)\n', ...
    R_opt.k_opt, R_opt.L0_opt, err_k_opt, err_L0_opt);
fprintf('  %d iteracions, %d avaluacions del model, %.1f ms, convergència: %d (nfev = %d)\n', ...
    R_opt.n_iter, R_opt.n_aval, dt_opt, R_opt.exit, R_opt.n_f);
fprintf('  forces: F_1 = %.6f N, F_2 = %.6f N\n', F_1o, F_2o);
fprintf('\n  Iteracions:\n');
for it = 1:size(R_opt.historial, 1)
    kk = R_opt.historial(it, 1); LL = R_opt.historial(it, 2);
    fprintf('  %3d  k = %8.4f N/mm  L0 = %8.4f mm  f = %.2e\n', it - 1, kk, LL, ...
        error_disseny(kk, LL, L1, L2, F1, F2));
end

%% Representació gràfica: camí de l'optimitzador
[kk_, LL_] = meshgrid(linspace(0.05, k_max, 300), linspace(L1, L0_max, 300));
f_grid = error_disseny(kk_, LL_, L1, L2, F1, F2);
figure('Position', [100 100 500 400])
hold on, grid on, set(gca, 'fontsize', 12)
[C, hc] = contour(kk_, LL_, f_grid, [1e-2 1e-1 1 10], ':', 'LineColor', blue, 'LineWidth', 0.8);
clabel(C, hc, 'FontSize', 8, 'Color', blue);
h0 = plot(NaN, NaN, ':', 'Color', blue, 'LineWidth', 0.8);
h1 = plot(k_, L0_F1, '--', 'Color', gris, 'LineWidth', 1.4);
h2 = plot(k_, L0_F2, '-.', 'Color', gris, 'LineWidth', 1.4);
H = R_opt.historial;
h3 = plot(H(:, 1), H(:, 2), 'o-', 'Color', vermell, 'MarkerSize', 4, 'LineWidth', 1.2, ...
    'MarkerFaceColor', vermell);
h4 = plot(k_an, L0_an, '+', 'Color', 'b', 'MarkerSize', 14, 'LineWidth', 1.8);
xlim([0 k_max]); ylim([L1 L0_max]);
xlabel('$k$ [N/mm]', 'Interpreter', 'latex')
ylabel('$L_0$ [mm]', 'Interpreter', 'latex')
legend([h0 h1 h2 h3 h4], {'$f$ constant', 'Corba $F_{(1)} = F_1$', 'Corba $F_{(2)} = F_2$', ...
    'Iteracions', 'Solució analítica'}, 'Interpreter', 'latex', 'Location', 'northeastoutside', ...
    'FontSize', 9)
print('DiAM_07a_molla_LSQ_Fig2', '-dpdf');

%% Funcions emprades
% (en MatLab les funcions locals van al final de l'script; en Python, al principi)

function [F_1, F_2] = molla(k, L0, L1, L2)
% Forces de la molla a les dues posicions de funcionament (llei de Hooke).
% Accepta escalars o vectors (k i L0 de la mateixa mida).
% Entrada: k [N/mm], L0, L1, L2 [mm]
% Sortida: forces F_1 i F_2 [N]
F_1 = k.*(L0 - L1);                             % vàlvula tancada (longitud L1)
F_2 = k.*(L0 - L2);                             % vàlvula oberta (longitud L2)
end

function [r_1, r_2] = residus(k, L0, L1, L2, F1, F2)
% Errors relatius de les forces respecte de les desitjades: r_1 i r_2 [-]
[F_1, F_2] = molla(k, L0, L1, L2);
r_1 = (F_1 - F1)/F1;
r_2 = (F_2 - F2)/F2;
end

function f = error_disseny(k, L0, L1, L2, F1, F2)
% Funció objectiu: suma dels quadrats dels errors relatius, f = r_1^2 + r_2^2
[r_1, r_2] = residus(k, L0, L1, L2, F1, F2);
f = r_1.^2 + r_2.^2;
end

function R = montecarlo(n, k_max, L0_max, L1, L2, F1, F2)
% Mètode de Montecarlo (Algorisme 1).
% Genera n parells (k, L0) uniformes a (0, k_max] x (L1, L0_max] i retorna
% el de menor error f.
k  = k_max*(1 - rand(n, 1));                    % uniforme a (0, k_max]
L0 = L1 + (L0_max - L1)*(1 - rand(n, 1));       % uniforme a (L1, L0_max]
f  = error_disseny(k, L0, L1, L2, F1, F2);
[f_opt, i_opt] = min(f);
R = struct('k', k, 'L0', L0, 'f', f, 'k_opt', k(i_opt), 'L0_opt', L0(i_opt), 'f_opt', f_opt);
end

function R = optimitzacio(x0, k_max, L0_max, L1, L2, F1, F2)
% Optimització per mínims quadrats no lineals (Algorisme 2).
% Minimitza f = r_1^2 + r_2^2 amb 0 < k <= k_max i L1 <= L0 <= L0_max.
% Retorna el disseny òptim, el nombre d'avaluacions i l'historial d'iteracions.
global N_AVAL HISTORIAL
N_AVAL = 0;                                     % comptador d'avaluacions del model
HISTORIAL = x0(:)';                             % punts de cada iteració (el primer és x0)

limits_inf = [1e-3, L1];
limits_sup = [k_max, L0_max];
opcions = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective', ...
    'Display', 'off', 'OutputFcn', @desa_iteracio);
tic;
[x, resnorm, ~, exitflag, output] = lsqnonlin(@(x) vector_residus(x, L1, L2, F1, F2), ...
    x0, limits_inf, limits_sup, opcions);
dt = toc*1e3;                                   % ms
R = struct('k_opt', x(1), 'L0_opt', x(2), 'f_opt', resnorm, 'n_aval', N_AVAL, ...
    'n_iter', size(HISTORIAL, 1) - 1, 'n_f', output.funcCount, 't', dt, ...
    'exit', exitflag > 0, 'historial', HISTORIAL);
end

function r = vector_residus(x, L1, L2, F1, F2)
% Vector d'errors r = [r_1; r_2] (el que minimitza lsqnonlin)
global N_AVAL
N_AVAL = N_AVAL + 1;
[r_1, r_2] = residus(x(1), x(2), L1, L2, F1, F2);
r = [r_1; r_2];
end

function stop = desa_iteracio(x, optimValues, state)
% Desa el punt de cada iteració en què el disseny canvia
global HISTORIAL
stop = false;
if strcmp(state, 'iter') && norm(x(:)' - HISTORIAL(end, :)) > 1e-6
    HISTORIAL(end+1, :) = x(:)';
end
end
