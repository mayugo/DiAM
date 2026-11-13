clear; clc; close all

% ---------------------------------------------------------------
% DADES (SI: m, kg, N, rad/s)
% ---------------------------------------------------------------
w_tab = [60 120 45 30];                 % rad/s
y0_tab = [1.5 0.9 1.4 1.1];             % mm
m_tab = [1.6 1.2 1.8 2.0];              % kg
X_lim = 5e-3;                           % amplitud desitjada de x
zeta = 0.02;
g = 9.81;

% ---------------------------------------------------------------
% CALCULS I SIMULACIO PER A LES 4 CONFIGURACIONS
% ---------------------------------------------------------------
for cfg = 1:4
    omega = w_tab(cfg);
    y0 = y0_tab(cfg)*1e-3;
    m = m_tab(cfg);

    [k, Omega] = rigidesa_necessaria(omega, m, zeta, X_lim, y0);
    wn = sqrt(k/m);
    Tr = transmitancia(Omega, zeta);
    [F_mitja, F_alt, F_alt_comp] = forces_suport(omega, m, zeta, k, X_lim, y0, g);

    [t, x, y, F, X_ss, F_ss] = DiAM_07e1_simulacio(m, k, zeta, y0, omega);

    fprintf('\nConfiguracio %d: w = %g rad/s, y0 = %.2f mm, m = %g kg\n', cfg, omega, y0*1e3, m);
    fprintf('  Omega = %.3f   wn = %.1f rad/s   T = %.2f\n', Omega, wn, Tr);
    fprintf('  k total = %.2f N/mm   k per suport = %.2f N/mm\n', k/1e3, k/2e3);
    fprintf('  F_m = %.2f N   F_a = %.2f N (comprovacio %.2f N)\n', F_mitja, F_alt, F_alt_comp);
    fprintf('  F_max = %.2f N   F_min = %.2f N\n', F_mitja + F_alt, F_mitja - F_alt);
    fprintf('  Simulacio: X_ss = %.2f mm, max|x| transitori = %.2f mm\n', X_ss*1e3, max(abs(x))*1e3);
    fprintf('             F_ss = %.2f N, max|F| transitori = %.2f N, F total max/min = %.2f / %.2f N\n', ...
        F_ss, max(abs(F)), F_mitja + max(F), F_mitja + min(F));

    % figura nomes per a la configuracio 1
    if cfg == 1
        figure('Position', [100 100 700 500])
        subplot(2, 1, 1)
        plot(t, x*1e3, t, y*1e3, 'LineWidth', 1), hold on
        yline(X_lim*1e3, 'k--'); yline(-X_lim*1e3, 'k--');
        ylabel('desplacament [mm]')
        legend('x (motlle)', 'y (base)'), grid on
        subplot(2, 1, 2)
        plot(t, F_mitja + F, 'LineWidth', 1), hold on
        yline(F_mitja, 'k--');
        xlabel('t [s]'), ylabel('F per suport [N]')
        xlim([0 4]), grid on
        exportgraphics(gcf, 'DiAM_07e1_SOL1_transitori.pdf')
    end
end

% ---------------------------------------------------------------
% FORMULACIO (funcions locals)
% ---------------------------------------------------------------
function T = transmitancia(Omega, zeta)
% Transmitancia de desplacament X/y0 amb esmorteiment
T = sqrt(1 + (2*zeta*Omega).^2)./sqrt((1 - Omega.^2).^2 + (2*zeta*Omega).^2);
end

function [k, Omega] = rigidesa_necessaria(omega, m, zeta, X, y0)
% Rigidesa total k tal que la transmitancia sigui X/y0 amb Omega > 1
Omega = fzero(@(Omega) transmitancia(Omega, zeta) - X/y0, [1 100]);
wn = omega/Omega;
k = m*wn^2;
end

function [F_mitja, F_alt, F_alt_comp] = forces_suport(omega, m, zeta, k, X, y0, g)
% Forces sobre cada suport: mitjana (pes) i alterna (dinamica)
F_mitja = m*g/2;
F_alt = m*omega^2*X/2;                  % forca d'inercia del motlle / 2
Omega = omega/sqrt(k/m);
Z = y0*Omega^2/sqrt((1 - Omega^2)^2 + (2*zeta*Omega)^2);
c = 2*zeta*sqrt(k*m);
F_alt_comp = sqrt(k^2 + (omega*c)^2)*Z/2;
end
