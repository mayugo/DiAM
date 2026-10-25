%% Disseny de la secció transversal d'un suport (columna)
% J.A. Mayugo, 2021 (revisió 2026)
%
% Falla per inestabilitat. Disseny paramètric amb dues variables de disseny (a, b).
%   Apartat 1): mètode de Montecarlo
%   Apartat 2): optimització amb restriccions (fmincon, algorisme SQP)
%   Apartat 3): comprovació validesa de la fórmula d'Euler

clear; close all;

%% Dades del problema
% Geometria [mm]
L     = 500;        % longitud del suport
a_max = 80;         % valor màxim de les variables de disseny a i b
% Càrrega [N]
P     = 20e3;       % força de compressió centrada
% Material: alumini 2014-T6 [MPa]
E     = 70e3;       % mòdul d'elasticitat
Sy    = 212;        % límit de fluència
% Especificacions de disseny
SF    = 2.5;        % coeficient de seguretat (resistència i vinclament)
% Condicions de contorn: longitud efectiva Le = K*L
K_xy  = 0.7;        % pla x-y: encastat-articulat
K_xz  = 2.0;        % pla x-z: encastat-lliure
% Mètode de Montecarlo
n_MC  = [1e3, 1e4, 1e5, 1e6];   % nombre de realitzacions
rng(2030, 'twister');           % llavor fixa per reproduir els resultats

%% Solució analítica (apèndix) per comparar
A_res = SF*P/Sy;                                        % àrea mínima per resistència
ratio = K_xy/K_xz;                                      % a/b del disseny més eficient
b_an  = (12*(K_xz*L)^2*SF*P/(pi^2*E*ratio))^0.25;
a_an  = ratio*b_an;
A_an  = a_an*b_an;
fprintf('Solució analítica: a = %.2f mm, b = %.2f mm, A = %.1f mm^2\n', a_an, b_an, A_an);
fprintf('  (àrea mínima per resistència: %.1f mm^2)\n\n', A_res);

%% Apartat 1): mètode de Montecarlo
fprintf('Apartat 1): mètode de Montecarlo\n');
fprintf('%9s %8s %8s %9s %8s %10s %11s\n', 'n', 'a [mm]', 'b [mm]', 'A [mm2]', ...
    'error A', 'factibles', 'temps [ms]');
montecarlo(1000, a_max, L, P, E, Sy, K_xy, K_xz, SF);   % escalfament (no es compta)
res = zeros(length(n_MC), 6);                           % [a, b, A, error, factibles, temps]
for k = 1:length(n_MC)
    n = n_MC(k);
    t0 = tic;
    R = montecarlo(n, a_max, L, P, E, Sy, K_xy, K_xz, SF);
    dt = toc(t0)*1e3;                                   % ms
    err = (R.A_opt - A_an)/A_an*100;
    res(k, :) = [R.a_opt, R.b_opt, R.A_opt, err, mean(R.factible)*100, dt];
    fprintf('%9d %8.2f %8.2f %9.1f %7.2f%% %9.1f%% %11.1f\n', n, res(k, :));
end
R_vis = montecarlo(1e4, a_max, L, P, E, Sy, K_xy, K_xz, SF);   % per al gràfic

%% Representació gràfica: dissenys factibles i no factibles
c_nofact = [1.000 0.498 0.055];     % taronja
c_fact   = [0.173 0.627 0.173];     % verd
c_MC     = [0.839 0.153 0.157];     % vermell
c_lim    = 0.65*[1 1 1];            % gris
f = R_vis.factible;

figure('Position', [100 100 600 550]);
hold on, box on
plot(R_vis.a(~f), R_vis.b(~f), '.', 'MarkerSize', 3, 'Color', c_nofact, 'HandleVisibility', 'off');
plot(R_vis.a(f),  R_vis.b(f),  '.', 'MarkerSize', 3, 'Color', c_fact,   'HandleVisibility', 'off');
plot(NaN, NaN, '.', 'MarkerSize', 15, 'Color', c_nofact, 'DisplayName', 'No factibles');  % llegenda
plot(NaN, NaN, '.', 'MarkerSize', 15, 'Color', c_fact,   'DisplayName', 'Factibles');
plot(R_vis.a_opt, R_vis.b_opt, 'x', 'Color', c_MC, 'MarkerSize', 10, 'LineWidth', 2.5, ...
    'DisplayName', 'Montecarlo');

% Frontera analítica de vinclament en cada pla i de resistència
a_    = linspace(1, a_max, 400);
b_xy  = SF*P*12*(K_xy*L)^2./(pi^2*E*a_.^3);             % Pcr_xy = SF*P
b_xz  = (SF*P*12*(K_xz*L)^2./(pi^2*E*a_)).^(1/3);       % Pcr_xz = SF*P
b_res = A_res./a_;                                      % A = A_res
plot(a_, b_xy,  '--', 'Color', c_lim, 'LineWidth', 1.4, 'DisplayName', 'Límit vinclament pla $x$-$y$');
plot(a_, b_xz,  '-.', 'Color', c_lim, 'LineWidth', 1.4, 'DisplayName', 'Límit vinclament pla $x$-$z$');
plot(a_, b_res, '-',  'Color', c_lim, 'LineWidth', 1.0, 'DisplayName', 'Límit resistència');
plot(a_an, b_an, 'k+', 'MarkerSize', 14, 'LineWidth', 1.8, 'DisplayName', 'Analítica');
plot(R_vis.a_opt, R_vis.b_opt, 'x', 'Color', c_MC, 'MarkerSize', 10, 'LineWidth', 2.5, ...
    'HandleVisibility', 'off');

xlim([0 a_max]); ylim([0 a_max]);
xlabel('$a$ [mm]', 'Interpreter', 'latex');
ylabel('$b$ [mm]', 'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.4);
exportgraphics(gcf, 'Fig1_DiAM_columna_MC.pdf', 'ContentType', 'vector');

%% Apartat 2): optimització amb restriccions (fmincon, SQP)
x0 = [a_max, a_max];                                        % punt inicial factible
optimitzacio(x0, a_max, L, P, E, Sy, K_xy, K_xz, SF);       % escalfament (no es compta)
t0 = tic;
R_opt = optimitzacio(x0, a_max, L, P, E, Sy, K_xy, K_xz, SF);
dt_opt = toc(t0)*1e3;                                       % ms
err_opt = (R_opt.A_opt - A_an)/A_an*100;
fprintf('\nApartat 2): optimització amb restriccions (SQP)\n');
fprintf('  punt inicial: a = %.0f mm, b = %.0f mm\n', x0);
fprintf('  solució: a = %.3f mm, b = %.3f mm, A = %.2f mm^2 (error %.4f%%)\n', ...
    R_opt.a_opt, R_opt.b_opt, R_opt.A_opt, abs(err_opt));
fprintf('  %d iteracions, %d avaluacions del model, %.1f ms, exitflag = %d\n', ...
    R_opt.n_iter, R_opt.n_aval, dt_opt, R_opt.exitflag);
fprintf('\n  Iteracions:\n');
for k = 1:size(R_opt.historial, 1)
    fprintf('  %3d  a = %7.3f mm  b = %7.3f mm  A = %8.2f mm^2\n', k-1, ...
        R_opt.historial(k, 1), R_opt.historial(k, 2), prod(R_opt.historial(k, :)));
end

%% Representació gràfica: camí de l'optimitzador
figure('Position', [100 100 840 550]);
hold on, box on
b_lim = max([b_xy; b_xz; b_res]);                           % zona factible
fill([a_, fliplr(a_)], [min(b_lim, a_max), a_max*ones(size(a_))], c_fact, ...
    'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'Zona factible');
plot(a_, b_xy,  '--', 'Color', c_lim, 'LineWidth', 1.4, 'DisplayName', 'Límit vinclament pla $x$-$y$');
plot(a_, b_xz,  '-.', 'Color', c_lim, 'LineWidth', 1.4, 'DisplayName', 'Límit vinclament pla $x$-$z$');
plot(a_, b_res, '-',  'Color', c_lim, 'LineWidth', 1.0, 'DisplayName', 'Límit resistència');
c_iso = [0.122 0.467 0.706];                                % blau
for Ak = [A_an, 1000, 2000, 4000]                           % corbes d'isoàrea
    plot(a_, Ak./a_, ':', 'Color', c_iso, 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
plot(NaN, NaN, ':', 'Color', c_iso, 'LineWidth', 0.8, 'DisplayName', 'Corbes $A = a\,b$ constant');
H = R_opt.historial;
plot(H(:, 1), H(:, 2), 'o-', 'Color', c_MC, 'MarkerSize', 4, 'MarkerFaceColor', c_MC, ...
    'LineWidth', 1.2, 'DisplayName', 'Iteracions SQP');
plot(a_an, b_an, 'k+', 'MarkerSize', 14, 'LineWidth', 1.8, 'DisplayName', 'Analítica');
xlim([0 a_max]); ylim([0 a_max]); axis square
xlabel('$a$ [mm]', 'Interpreter', 'latex');
ylabel('$b$ [mm]', 'Interpreter', 'latex');
legend('Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.4);
exportgraphics(gcf, 'Fig2_DiAM_columna_SLSQP.pdf', 'ContentType', 'vector');

%% Apartat 3): validesa de la fórmula d'Euler: tensió crítica en funció de l'esveltesa
% Columna llarga (Euler) si Le/r >= (Le/r)_1; columna intermèdia (Johnson) si no
lam_1  = sqrt(2*pi^2*E/Sy);                         % esveltesa límit Euler-Johnson
lam_   = linspace(1, 200, 800);
sig_E  = pi^2*E./lam_.^2;                           % Euler
sig_J  = Sy - (Sy*lam_/(2*pi)).^2/E;                % Johnson
a_d = R_opt.a_opt; b_d = R_opt.b_opt;               % disseny òptim (SQP)
lam_xy = K_xy*L*sqrt(12)/a_d;                       % esveltesa pla x-y (r_z = a/sqrt(12))
lam_xz = K_xz*L*sqrt(12)/b_d;                       % esveltesa pla x-z (r_y = b/sqrt(12))
sig_d  = SF*P/(a_d*b_d);                            % tensió crítica de disseny, Pcr/A
fprintf('\nApartat 3): validesa de la fórmula d''Euler\n');
fprintf('  esveltesa límit (Le/r)_1 = %.1f\n', lam_1);
fprintf('  esveltesa pla x-y = %.1f, pla x-z = %.1f\n', lam_xy, lam_xz);
fprintf('  tensió crítica de disseny = %.1f MPa (Sy/2 = %.0f MPa)\n', sig_d, Sy/2);
if min(lam_xy, lam_xz) >= lam_1
    fprintf('  -> columna llarga: Euler és aplicable\n');
else
    fprintf('  -> columna intermèdia: cal aplicar Johnson\n');
end

c_blau = [0.122 0.467 0.706];                       % blau
figure('Position', [100 100 650 500]);
hold on, box on
plot(lam_(lam_ >= lam_1), sig_E(lam_ >= lam_1), '-',  'Color', c_lim, 'LineWidth', 2.8, ...
    'DisplayName', 'F\''ormula d''Euler');
plot(lam_(lam_ <  lam_1), sig_E(lam_ <  lam_1), '--', 'Color', c_lim, 'LineWidth', 1.8, ...
    'HandleVisibility', 'off');
plot(lam_(lam_ <= lam_1), sig_J(lam_ <= lam_1), '-',  'Color', c_blau, 'LineWidth', 2.8, ...
    'DisplayName', 'F\''ormula de Johnson');
yline(Sy,   '--', 'Color', c_blau, 'LineWidth', 2.0, 'HandleVisibility', 'off');
yline(Sy/2, '--', 'Color', c_blau, 'LineWidth', 2.0, 'HandleVisibility', 'off');
xline(lam_1, '-.', 'Color', c_lim, 'LineWidth', 2.8, 'HandleVisibility', 'off');   % límit entre zones
plot(lam_xz, sig_d, 'o', 'Color', c_MC, 'MarkerFaceColor', c_MC, 'MarkerSize', 8, ...
    'DisplayName', 'Disseny, pla $x$-$z$');
plot(lam_xy, sig_d, 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Disseny, pla $x$-$y$');

text(4, Sy/2 + 8, sprintf('$S_y/2 = %.0f$ MPa', Sy/2), 'Interpreter', 'latex', 'Color', c_blau);
text(4, Sy + 8,   sprintf('$S_y = %.0f$ MPa', Sy),     'Interpreter', 'latex', 'Color', c_blau);
text(lam_1 + 2, 10, sprintf('$(L_e/r)_1 = %.1f$', lam_1), 'Interpreter', 'latex');
text(lam_xz + 8, sig_d - 12, {sprintf('$L_e/r = %.1f$', lam_xz), ...
    sprintf('$\\sigma_{cr} = %.1f$ MPa', sig_d)}, 'Interpreter', 'latex');
% Cotes de les dues zones (cal fixar els límits abans de dibuixar-les)
xlim([0 200]); ylim([0 300]);
y_cota = 1.22*Sy;
annotation_doble(0, lam_1, y_cota);
annotation_doble(lam_1, 200, y_cota);
text(lam_1/2, y_cota + 4, 'Columna interm\`edia', 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(lam_1/2, y_cota - 4, '(Johnson)', 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text((lam_1 + 200)/2, y_cota + 4, 'Columna llarga', 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text((lam_1 + 200)/2, y_cota - 4, '(Euler)', 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

xlabel('Esveltesa $L_e/r$ [--]', 'Interpreter', 'latex');
ylabel('Tensi\''o cr\''itica $\sigma_{cr} = P_{cr}/A$ [MPa]', 'Interpreter', 'latex');
legend('Location', 'east', 'Interpreter', 'latex', 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.4);
exportgraphics(gcf, 'Fig3_DiAM_columna_Euler_Johnson.pdf', 'ContentType', 'vector');

%% Funcions emprades
% (en MatLab les funcions locals van al final de l'script; en Python, al principi)

function [A, SF_res, SF_cr, Pcr] = suport(a, b, L, P, E, Sy, K_xy, K_xz)
% Coeficients de seguretat del suport per resistència i per vinclament.
% Accepta escalars o vectors (a i b de la mateixa mida).
% Entrada: a, b [mm], L [mm], P [N], E, Sy [MPa], K_xy, K_xz [-]
% Sortida: A [mm^2], SF_res, SF_cr [-] i Pcr [N]
A = a.*b;                                   % àrea de la secció
% Resistència a compressió
sigma  = P./A;                              % MPa, tensió de compressió
SF_res = Sy./sigma;                         % coeficient de seguretat per resistència
% Vinclament en el pla x-y (flexió al voltant de z, encastat-articulat)
I_z    = b.*a.^3/12;
Pcr_xy = pi^2*E*I_z/(K_xy*L)^2;
% Vinclament en el pla x-z (flexió al voltant de y, encastat-lliure)
I_y    = a.*b.^3/12;
Pcr_xz = pi^2*E*I_y/(K_xz*L)^2;
% Càrrega crítica del suport: la menor dels dos plans
Pcr    = min(Pcr_xy, Pcr_xz);
SF_cr  = Pcr/P;                             % coeficient de seguretat per vinclament
end

function R = montecarlo(n, a_max, L, P, E, Sy, K_xy, K_xz, SF)
% Mètode de Montecarlo (Algorisme 1).
% Genera n parells (a, b) uniformes a (0, a_max], selecciona els factibles
% i retorna el de secció mínima.
a = a_max*(1 - rand(n, 1));                 % uniforme a (0, a_max]
b = a_max*(1 - rand(n, 1));
[A, SF_res, SF_cr, Pcr] = suport(a, b, L, P, E, Sy, K_xy, K_xz);
factible = (SF_res >= SF) & (SF_cr >= SF);  % restriccions de disseny
i_fact = find(factible);
[~, j] = min(A(factible));
i_opt = i_fact(j);
R = struct('a', a, 'b', b, 'A', A, 'Pcr', Pcr, 'factible', factible, ...
           'a_opt', a(i_opt), 'b_opt', b(i_opt), 'A_opt', A(i_opt), 'Pcr_opt', Pcr(i_opt));
end

function R = optimitzacio(x0, a_max, L, P, E, Sy, K_xy, K_xz, SF)
% Optimització amb restriccions amb fmincon i l'algorisme SQP (Algorisme 2).
% min A = a*b  subjecte a  SF_res >= SF,  SF_cr >= SF  i  0 < a, b <= a_max
% Retorna el disseny òptim, el nombre d'avaluacions del model i
% l'historial de les iteracions.
n_aval = 0;                                 % comptador d'avaluacions del model
historial = [];                             % punts de cada iteració

    function A = objectiu(x)                % funció objectiu: àrea de la secció
        n_aval = n_aval + 1;
        A = suport(x(1), x(2), L, P, E, Sy, K_xy, K_xz);
    end

    function [c, ceq] = restriccions(x)     % restriccions: c(x) <= 0
        n_aval = n_aval + 1;
        a = x(1); b = x(2);
        [~, SF_res] = suport(a, b, L, P, E, Sy, K_xy, K_xz);
        Pcr_xy = pi^2*E*(b*a^3/12)/(K_xy*L)^2;
        Pcr_xz = pi^2*E*(a*b^3/12)/(K_xz*L)^2;
        c = [SF - Pcr_xy/P;                 % vinclament pla x-y
             SF - Pcr_xz/P;                 % vinclament pla x-z
             SF - SF_res];                  % resistència
        ceq = [];
    end

    function stop = desa(x, ~, ~)           % desa el punt de cada iteració
        historial = [historial; x(:)'];
        stop = false;
    end

opcions = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
                       'OptimalityTolerance', 1e-6, 'ConstraintTolerance', 1e-6, ...
                       'OutputFcn', @desa);
[x, A_opt, exitflag, output] = fmincon(@objectiu, x0, [], [], [], [], ...
    [1e-3, 1e-3], [a_max, a_max], @restriccions, opcions);
R = struct('a_opt', x(1), 'b_opt', x(2), 'A_opt', A_opt, 'n_aval', n_aval, ...
           'n_iter', output.iterations, 'exitflag', exitflag, 'historial', historial);
end

function annotation_doble(x1, x2, y)
% Dibuixa una cota (fletxa de doble cap) entre x1 i x2 a l'altura y,
% en coordenades de les dades de l'eix actual
plot([x1 x2], [y y], 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
dx = 0.02*diff(xlim); dy = 0.012*diff(ylim);
fill([x1, x1+dx, x1+dx], [y, y+dy, y-dy], 'k', 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([x2, x2-dx, x2-dx], [y, y+dy, y-dy], 'k', 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
