%% Implementació cas d'estudi biga corba (suport en U amb perfil en I)
% J.A. Mayugo, 2020 (revisió 2026)
%
% Falla estàtica. Disseny paramètric d'un element mecànic.
%   Apartats 1) i 2): tensions i factor de seguretat amb t = 4 mm i F = 3000 N
%   Apartat 3):       gruix t per a F = 2000 N amb el mateix SF (bisecció)
%   Apartat 4):       gruix t per a diverses forces F (bisecció)

clear; close all;

%% Dades del problema (Taules 1, 8 i 9)
% Geometria [mm]
d   = 80;       % distància de la força a la part interior U
r_i = 40;       % radi interior suport U
r_o = 64;       % radi exterior suport U
b_i = 20;       % amplada ala interior secció I
b_o = 10;       % amplada ala exterior secció I
t0  = 4;        % gruix perfil secció I (disseny inicial)
% Material [MPa]
Sut = 820;      % resistència última a tracció
Sy  = 600;      % límit de fluència
% Càrregues [N]
F0 = 3000;                          % apartats 1) i 2)
F3 = 2000;                          % apartat 3)
F4 = [1600, 2000, 2400, 2800, 3200];   % apartat 4)
% Paràmetres del mètode de la bisecció
t_min0 = 0.5;   % mm, gruix mínim de l'interval inicial
t_max0 = 10;    % mm, gruix màxim de l'interval inicial
SF_eps = 5e-3;  % tolerància pel factor de seguretat
inc    = 0.1;   % mm, resolució de fabricació del gruix

%% Apartats 1) i 2): F = F0 i t = t0
[SF_d, R1] = DiAM_f_beam_I(F0, d, r_i, r_o, t0, b_i, b_o, Sy, 1);

fprintf('Apartats 1) i 2): F = %.0f N i t = %.1f mm\n', F0, t0);
fprintf('  A = %.0f mm^2, I = %.0f mm^4, r_c = %.2f mm, r_n = %.2f mm, e = %.2f mm\n', ...
    R1.A, R1.I, R1.r_c, R1.r_n, R1.e);
fprintf('  N = %.0f N, M = %.0f N·mm\n', R1.N, R1.M);
punts = 'abcd';
for k = 1:4
    fprintf('  punt %s: sigma = %7.1f MPa, SF = %.2f\n', ...
        punts(k), R1.sigma(k), R1.SF(k));
end
fprintf('  SF de la secció = %.2f\n\n', SF_d);

%% Apartat 3): F = F3 amb el mateix SF (bisecció)
[t3, n3, SF3] = biseccio(F3, SF_d, d, r_i, r_o, b_i, b_o, Sy, ...
                         t_min0, t_max0, SF_eps, inc);
fprintf('Apartat 3): amb F = %.0f N, %d iteracions, t = %.1f mm i SF = %.2f\n\n', ...
    F3, n3, t3, SF3);

%% Apartat 4): diverses forces (bisecció)
t4  = zeros(size(F4));
n4  = zeros(size(F4));
SF4 = zeros(size(F4));
fprintf('Apartat 4): mètode de la bisecció\n');
for k = 1:length(F4)
    [t4(k), n4(k), SF4(k)] = biseccio(F4(k), SF_d, d, r_i, r_o, b_i, b_o, Sy, ...
                                      t_min0, t_max0, SF_eps, inc);
    fprintf('  F = %.0f N: %d iteracions, t = %.1f mm, SF = %.2f\n', ...
        F4(k), n4(k), t4(k), SF4(k));
end

% Comprovació amb el mètode incremental (augmentar t en passos de inc)
fprintf('\nComprovació: mètode incremental\n');
for F = F4
    t = t_min0;
    SF_trial = DiAM_f_beam_I(F, d, r_i, r_o, t, b_i, b_o, Sy, 0);
    cont = 0;
    while SF_trial < SF_d
        cont = cont + 1;
        t = round(t + inc, 6);   % evita errors d'arrodoniment
        SF_trial = DiAM_f_beam_I(F, d, r_i, r_o, t, b_i, b_o, Sy, 0);
    end
    fprintf('  F = %.0f N: %d iteracions, t = %.1f mm, SF = %.2f\n', ...
        F, cont, t, SF_trial);
end

%% Representació gràfica de la solució de l'apartat 4)
% La relació entre F i t no és lineal, comportament típic d'elements
% sotmesos principalment a flexió.
figure
hold on, grid on, set(gca, 'fontsize', 20)
plot(F4, t4, 'x', 'MarkerSize', 12, 'LineWidth', 2.8);
xlabel('$F$ [N]', 'Interpreter', 'latex', 'fontsize', 24)
ylabel('$t$ [mm]', 'Interpreter', 'latex', 'fontsize', 24)

%% Funcions emprades
% (en MatLab les funcions locals van al final de l'script; en Python, al principi)

function [r_c, r_n, A, I] = f_curved_beam_I(r_i, r_o, t, t_i, t_o, b_i, b_o)
% Propietats de la secció en I (Taula 2).
% Entrada: r_i, r_o, t, t_i, t_o, b_i, b_o (segons figura secció) [mm]
% Sortida: radi centroide r_c, radi línia neutra r_n [mm],
%          àrea A [mm^2] i moment de segon ordre I [mm^4]
h   = r_o - r_i;                                        % cantell
A   = h*t + t_i*(b_i - t) + t_o*(b_o - t);              % àrea
r_c = r_i + (1/2*h^2*t + 1/2*t_i^2*(b_i - t) ...
             + t_o*(b_o - t)*(h - t_o/2)) / A;          % radi centroide
r_n = A / (b_i*log((r_i + t)/r_i) ...
           + t*log((r_o - t_o)/(r_i + t_i)) ...
           + b_o*log(r_o/(r_o - t_o)));                 % radi línia neutra
I   = h^3*t/12           + h*t*(r_c - r_i - h/2)^2 ...
    + t_i^3*(b_i - t)/12 + t_i*(b_i - t)*(r_c - r_i - t_i/2)^2 ...
    + t_o^3*(b_o - t)/12 + t_o*(b_o - t)*(r_c - r_i - h + t_o/2)^2;
end

function [SFe, R] = DiAM_f_beam_I(F, d, r_i, r_o, t, b_i, b_o, Sy, fig)
% Tensions als punts a, b, c, d i factor de seguretat estàtic (Taules 2-6).
% Entrada: F [N], geometria [mm], gruix t [mm], Sy [MPa]
%          fig=1 dibuixa les tensions a la secció, fig=0 no ho fa
% Sortida: SFe, factor de seguretat estàtic de la secció [-]
%          R,   estructura amb resultats intermedis

% Propietats de la secció de perfil en I (t_i = t_o = t)
[r_c, r_n, A, I] = f_curved_beam_I(r_i, r_o, t, t, t, b_i, b_o);
e   = r_c - r_n;
c_i = r_n - r_i;
c_o = r_o - r_n;
y_i = r_c - r_i;
y_o = r_o - r_c;

% Sol·licitacions a la secció crítica (Taula 3)
N = F;              % N,    força axial
M = (d + y_i)*F;    % N·mm, moment flector

% Tensions als 4 punts crítics (Taula 4) [MPa]
%   'a' superior i 'b' inferior: hipòtesi biga corba
%   'c' superior i 'd' inferior: hipòtesi biga recta
sig_a = N/A + M*c_i/(A*e*r_i);
sig_b = N/A - M*c_o/(A*e*r_o);
sig_c = N/A + M*y_i/I;
sig_d = N/A - M*y_o/I;
sigma = [sig_a, sig_b, sig_c, sig_d];

% Factor de seguretat estàtic (Taula 6)
SF  = Sy./abs(sigma);
SFe = min(SF);

R = struct('h', r_o - r_i, 'A', A, 'I', I, 'r_c', r_c, 'r_n', r_n, 'e', e, ...
           'N', N, 'M', M, 'sigma', sigma, 'SF', SF);

% Distribució de tensions a la secció crítica (Taula 5)
if fig == 1
    r_ = r_i:0.1:r_o;
    sigma_rect = N/A + M*(r_c - r_)/I;
    sigma_corb = N/A + M*(r_n - r_)./(A*e*r_);
    figure;
    hold on, grid on, set(gca, 'fontsize', 20)
    plot(sigma_corb, -r_, 'LineWidth', 2.0);
    plot(sigma_rect, -r_, 'r', 'LineWidth', 2.0);
    plot([-400, 400], [-r_n, -r_n], '--b', 'LineWidth', 1.2)
    plot([-400, 400], [-r_c, -r_c], '--r', 'LineWidth', 1.2)
    plot(sigma, [-r_i, -r_o, -r_i, -r_o], 'ok', 'LineWidth', 2.8, 'MarkerSize', 8)
    legend({'$\sigma$ biga corba', '$\sigma$ biga recta', '$r_n$', '$r_c$'}, ...
        'Interpreter', 'latex', 'Location', 'Best', 'fontsize', 20)
    xlabel('$\sigma$ [MPa]', 'Interpreter', 'latex', 'fontsize', 24)
    ylabel('$r$ [mm]', 'Interpreter', 'latex', 'fontsize', 24)
end
end

function [t_inc, cont, SF_inc] = biseccio(F, SF_d, d, r_i, r_o, b_i, b_o, Sy, ...
                                          t_min, t_max, SF_eps, inc)
% Mètode de la bisecció (Algorisme 1).
% Retorna el gruix arrodonit amunt a la resolució inc, el nombre
% d'iteracions i el SF obtingut amb aquest gruix arrodonit.
SF_trial = 0;
cont = 0;
t = 0;
while abs(SF_d - SF_trial) > SF_eps         % repeteix fins a convergir
    cont = cont + 1;
    t = (t_min + t_max)/2;                  % punt mig de l'interval
    SF_trial = DiAM_f_beam_I(F, d, r_i, r_o, t, b_i, b_o, Sy, 0);
    if SF_d >= SF_trial
        t_min = t;                          % cal més gruix
    else
        t_max = t;                          % sobra gruix
    end
end
t_inc  = ceil(t/inc - 1e-9)*inc;            % arrodonit amunt a inc
SF_inc = DiAM_f_beam_I(F, d, r_i, r_o, t_inc, b_i, b_o, Sy, 0);
end
