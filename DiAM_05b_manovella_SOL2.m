%% Analisi a fatiga de la manovella (apartat 2)
% Cicle de carrega: la forca F es desplaca des de d = 20 mm fins a d = 60 mm i torna.
% Per a cada punt estudiat (seccio 1 punts a i b, seccio 2 punts a i b):
%   1) tensions principals (I, II, III) a cada distancia d (instant de carrega)
%   2) components mitjana i alterna de cada tensio principal al llarg del cicle
%   3) tensions equivalents (von Mises) mitjana i alterna
%   4) factor de seguretat a fatiga (Goodman modificada)
%   5) comprovacio de la fluencia (von Mises) a l'instant mes desfavorable
% J.A. Mayugo (2026)
clear; close all; clc

%% Dades del problema
F   = 150;            % N, forca aplicada
d   = 20:1:60;        % mm, distancies del cicle (cada 1 mm)
Sut = 150;            % MPa, resistencia ultima
Sy  = 70;             % MPa, limit elastic
Se  = 45;             % MPa, limit de fatiga (ja corregit)
t   = 4;              % mm, gruix (seccions 1 i 2)

sec_name = {'Seccio 1', 'Seccio 2'};
s   = [180, 15];                     % mm, cota de cada seccio
b   = [25, 20];                      % mm, amplada de cada seccio
eta = [0.300, 0.306; 0.291, 0.306];  % coeficients de torsio de cada seccio
pts = 'ab';                          % punts analitzats

%% Calcul
n = 0;
for i = 1:numel(sec_name)
    h = 30 + (40-30)/200*s(i);
    [Ix, It, c] = section_properties(b(i), h, t, eta(i,1), eta(i,2));
    for ip = 1:numel(pts)
        n = n + 1;
        [V, M, T]  = loads(F, d, s(i), c, t);
        [sig, tau] = nominal_stresses(V, M, T, b(i), h, t, Ix, It, pts(ip));
        S = principal_stresses(sig, tau);                 % 3 x numel(d)
        [S_m, S_a] = cycle_components(S);                 % 3 components
        [sig_m, sig_a] = equivalent_stresses(S_m, S_a);
        [seqv_max, i_crit] = max(von_mises(S));           % instant mes desfavorable
        res(n).label = sprintf('%s, punt %c', sec_name{i}, pts(ip));
        res(n).S     = S;
        res(n).S_m   = S_m;
        res(n).S_a   = S_a;
        res(n).sig_m = sig_m;
        res(n).sig_a = sig_a;
        res(n).SFf   = goodman(sig_a, sig_m, Se, Sut);
        res(n).SFy   = Sy/seqv_max;
        res(n).d_crit = d(i_crit);
    end
end

%% Resultats
for k = 1:numel(res)
    fprintf('%s: sigma_m = %5.1f MPa, sigma_a = %5.1f MPa, SFf = %4.2f, SFy = %4.2f (instant mes desfavorable: d = %.0f mm)\n', ...
        res(k).label, res(k).sig_m, res(k).sig_a, res(k).SFf, res(k).SFy, res(k).d_crit);
    fprintf('    components mitjanes (I, II, III) = %6.1f, %6.1f, %6.1f MPa\n', res(k).S_m);
    fprintf('    components alternes (I, II, III) = %6.1f, %6.1f, %6.1f MPa\n', res(k).S_a);
end

%% Representacio del cicle complet (anada 20 -> 60 mm i tornada 60 -> 20 mm)
d_cicle = [d, fliplr(d(1:end-1))];        % la tornada repeteix els valors
x       = 0:numel(d_cicle)-1;             % instant del cicle
x_cota  = max(x) + 4;                     % posicio de la cota (fora del requadre)
tol     = 0.05;                           % MPa, tensio considerada nul.la
colors  = lines(3);
roman   = {'I', 'II', 'III'};
fig = figure('Position', [100 50 800 1000]);
for k = 1:numel(res)
    subplot(4, 1, k); hold on; grid on
    S_cicle = [res(k).S, fliplr(res(k).S(:, 1:end-1))];
    for j = 1:3
        sm = res(k).S_m(j);  sa = res(k).S_a(j);
        plot(x, S_cicle(j,:), '-', 'Color', colors(j,:), 'LineWidth', 2, ...
            'DisplayName', ['\sigma_{' roman{j} '}']);
        plot([min(x) max(x)], [sm sm], '--', 'Color', colors(j,:), ...
            'HandleVisibility', 'off');                       % component mitjana
        if abs(sm) < tol && abs(sa) < tol                     % tensio nul.la
            txt = {['\sigma_{' roman{j} ',m} = 0'], ['\sigma_{' roman{j} ',a} = 0']};
        else
            [~, ie] = max(abs(S_cicle(j,:) - sm));
            ext = S_cicle(j, ie);                             % extrem del cicle
            plot([max(x) x_cota], [ext ext], ':', 'Color', colors(j,:), ...
                'HandleVisibility', 'off', 'Clipping', 'off');
            plot([x_cota x_cota], [sm ext], '-', 'Color', colors(j,:), ...
                'LineWidth', 1.5, 'HandleVisibility', 'off', 'Clipping', 'off'); % cota alterna
            plot(x_cota, ext, '^', 'Color', colors(j,:), 'MarkerSize', 5, ...
                'MarkerFaceColor', colors(j,:), 'HandleVisibility', 'off', 'Clipping', 'off');
            txt = {sprintf('\\sigma_{%s,m} = %.1f', roman{j}, sm), ...
                   sprintf('\\sigma_{%s,a} = %.1f', roman{j}, sa)};
        end
        text(x_cota + 1.5, sm, txt, 'Color', colors(j,:), 'FontSize', 12, ...
            'VerticalAlignment', 'middle');
    end
    yline(0, 'k-', 'HandleVisibility', 'off');
    xlim([min(x) max(x)]);
    xticks(0:10:max(x)); xticklabels(arrayfun(@(v) sprintf('%.0f', v), d_cicle(1:10:end), 'UniformOutput', false));
    title(res(k).label); ylabel('Tensio (MPa)');
    set(gca, 'FontSize', 12);
    p = get(gca, 'Position'); p(3) = 0.68*p(3); set(gca, 'Position', p);   % espai per a les etiquetes
    if k == numel(res)
        xlabel('Instant del cicle, d (mm): anada 20 -> 60 mm i tornada 60 -> 20 mm');
    end
end
exportgraphics(fig, 'DiAM_05b_manovella_SOL2_cicle.pdf');   % R2020a o posterior

%% Esquerra: estatic (tensions principals maximes i limit de von Mises)
%% Dreta: fatiga (diagrama de vida infinita, tipus Haigh)
fig2 = figure('Position', [100 100 1300 600]);
pcol  = lines(7);  pcol = pcol(4:7, :);
marks = {'o', 's', '^', 'd'};
msz   = [13, 10, 13, 10];          % mides diferents perque els punts propers es vegin tots
names_p = {'1a', '1b', '2a', '2b'};

% a) Estatic: tensions principals no nul.les (A = I, B = III) a l'instant mes
%    desfavorable (d = 60 mm) i limit de von Mises
%    (estat pla: sigma_A^2 - sigma_A sigma_B + sigma_B^2 = Sy^2)
subplot(1, 2, 1); hold on; grid on
[sA, sB] = meshgrid(linspace(-5, 90, 400), linspace(-90, 5, 400));
SF_lines = [1 1.5 2 3];                          % linies de factor de seguretat
for SF = SF_lines                              % linies de SF constant (blau fi)
    contour(sA, sB, sA.^2 - sA.*sB + sB.^2 - (Sy/SF)^2, [0 0], 'LineColor', [0 0.45 0.74], 'LineWidth', 0.8);
    text(0.5, -Sy/SF, sprintf('SF = %g', SF), 'Color', [0 0.45 0.74], 'FontSize', 11, ...
        'BackgroundColor', 'w', 'Margin', 0.5);
end
orange = [1 0.55 0];
contour(sA, sB, sA.^2 - sA.*sB + sB.^2 - Sy^2, [0 0], 'LineColor', orange, 'LineWidth', 5);   % limit de von Mises
for k = 1:numel(res)
    [~, ic] = max(von_mises(res(k).S));                  % instant mes desfavorable
    hp(k) = plot(res(k).S(1,ic), res(k).S(3,ic), marks{k}, 'Color', pcol(k,:), ...
        'MarkerFaceColor', pcol(k,:), 'MarkerSize', msz(k)-1);
end
xline(0, 'k-'); yline(0, 'k-');
axis equal; xlim([-5 90]); ylim([-90 5]);
xlabel('\sigma_A = \sigma_I (MPa)'); ylabel('\sigma_B = \sigma_{III} (MPa)');
legend(hp, strcat('Punt ', names_p), 'Location', 'southeast');
set(gca, 'FontSize', 12);

% b) Fatiga: diagrama de vida infinita (Haigh) amb Goodman modificada (plataforma
%    sigma_a = Se per a sigma_m < 0)
subplot(1, 2, 2); hold on; grid on
for SF = SF_lines                              % linies de SF constant (blau fi)
    plot([-20 0 Sut/SF], [Se/SF Se/SF 0], '-', 'Color', [0 0.45 0.74], 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
    text(-10, Se/SF, sprintf('SF = %g', SF), 'Color', [0 0.45 0.74], 'FontSize', 11, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
plot([-20 0 Sut], [Se Se 0], '-', 'Color', orange, 'LineWidth', 5, 'DisplayName', 'Goodman modificada');
for k = 1:numel(res)
    plot(res(k).sig_m, res(k).sig_a, marks{k}, 'Color', pcol(k,:), ...
        'MarkerSize', msz(k), 'LineWidth', 2, 'DisplayName', ['Punt ' names_p{k}]);
end
xlim([-20 160]); ylim([0 80]);
xlabel('\sigma_m (MPa)'); ylabel('\sigma_a (MPa)');
legend('Location', 'northeast'); set(gca, 'FontSize', 12);
exportgraphics(fig2, 'DiAM_05b_manovella_SOL2_haigh.pdf');

%% Funcions
function [Ix, It, c] = section_properties(b, h, t, eta_a, eta_b)
% Inertia a flexio Ix [mm^4], constant de torsio It [mm^4] i centre de
% tallant c [mm] de la seccio en U (amplada b, cantell h, gruix t)
alpha = 1.1;
Ix = b*h^3/12 - (b-t)*(h-2*t)^3/12;
It = alpha*(2*eta_a*b*t^3 + eta_b*(h-2*t)*t^3);
c  = (b-t/2)^2*(h-t)^2*t/4/Ix;
end

function [V, M, T] = loads(F, d, s, c, t)
% Tallant V [N], moment flector M [N mm] i moment torsor T [N mm]
% La distancia d pot ser un vector (un valor per cada instant del cicle)
V = F;
M = F*s;
T = F*(d + t/2 - c);
end

function [sig, tau] = nominal_stresses(V, M, T, b, h, t, Ix, It, point)
% Tensio normal per flexio sig i tallant total tau [MPa] al punt 'a'
% (fibra extrema) o 'b' (eix neutre)
tau_T = T*t/It;
if point == 'a'
    sig = M*(h/2)/Ix*ones(size(tau_T));
    Me  = (b-t)*t*(h-t)/2;
else
    sig = zeros(size(tau_T));
    Me  = b*t*(h-t)/2 + (h/2-t)*t*(h/2-t)/2;
end
tau_V = V*Me/t/Ix;
tau   = tau_V + tau_T;
end

function S = principal_stresses(sig, tau)
% Tensions principals ordenades I >= II >= III (estat pla: una val 0)
% Retorna una matriu 3 x n (una columna per cada instant del cicle)
R = sqrt((sig/2).^2 + tau.^2);
S = [sig/2 + R; zeros(size(sig)); sig/2 - R];
S = sort(S, 1, 'descend');
end

function [S_m, S_a] = cycle_components(S)
% Components mitjana i alterna de cada tensio principal (fila de S) al
% llarg del cicle
Smax = max(S, [], 2);
Smin = min(S, [], 2);
S_m  = (Smax + Smin)/2;
S_a  = (Smax - Smin)/2;
end

function seqv = von_mises(S)
% Tensio equivalent de von Mises a partir de les tensions principals
seqv = sqrt(2)/2*sqrt((S(1,:)-S(2,:)).^2 + (S(2,:)-S(3,:)).^2 ...
    + (S(3,:)-S(1,:)).^2);
end

function [sig_m, sig_a] = equivalent_stresses(S_m, S_a)
% Tensions equivalents mitjana (amb signe de la suma de tensions principals)
% i alterna
if sum(S_m) >= 0, signe = 1; else, signe = -1; end   % si la suma es 0, signe positiu
sig_m = signe*von_mises(S_m);
sig_a = von_mises(S_a);
end

function SFf = goodman(sig_a, sig_m, Se, Sut)
% Factor de seguretat a fatiga, regla de Goodman modificada
if sig_m > 0
    SFf = 1/(sig_a/Se + sig_m/Sut);
else
    SFf = Se/sig_a;
end
end
