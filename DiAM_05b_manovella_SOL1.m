%% Analisi estatica de la manovella (apartat 1)
% Factor de seguretat estatic (von Mises) a la seccio 1 (punts a i b) i a la
% seccio 2 (punts a i b), per a una distancia d donada de la forca.
% J.A. Mayugo (2026)
clear; close all; clc

%% Dades del problema
F       = 150;           % N, forca aplicada
d       = 20;            % mm, distancia entre la forca i el perfil en U
Sy      = 70;            % MPa, limit elastic de l'alumini
t       = 4;             % mm, gruix (seccions 1 i 2)

sec_name = {'Seccio 1', 'Seccio 2'};
s   = [180, 15];                     % mm, cota de cada seccio
b   = [25, 20];                      % mm, amplada de cada seccio
eta = [0.300, 0.306; 0.291, 0.306];  % coeficients de torsio de cada seccio
pts = 'ab';                          % punts analitzats

%% Calcul i resultats
for i = 1:numel(sec_name)
    h = 30 + (40-30)/200*s(i);                          % mm, cantell
    [Ix, It, c] = section_properties(b(i), h, t, eta(i,1), eta(i,2));
    for ip = 1:numel(pts)
        [V, M, T]  = loads(F, d, s(i), c, t);
        [sig, tau] = nominal_stresses(V, M, T, b(i), h, t, Ix, It, pts(ip));
        SFe = static_safety_factor(sig, tau, Sy);
        fprintf('%s, punt %c: sigma = %5.1f MPa, tau = %5.1f MPa, sigma_eq = %5.1f MPa, SFe = %4.2f\n', ...
            sec_name{i}, pts(ip), sig, tau, von_mises(sig, tau), SFe);
    end
end

%% Funcions
function [Ix, It, c] = section_properties(b, h, t, eta_a, eta_b)
% Propietats de la seccio en U
% Inputs :  amplada b, cantell h, gruix t [mm] i coeficients de torsio
%           (eta_a, eta_b)
% Outputs:  inercia a flexio Ix [mm^4], constant de torsio It [mm^4] i
%           centre de tallant c [mm]
alpha = 1.1;
Ix = b*h^3/12 - (b-t)*(h-2*t)^3/12;
It = alpha*(2*eta_a*b*t^3 + eta_b*(h-2*t)*t^3);
c  = (b-t/2)^2*(h-t)^2*t/4/Ix;
end

function [V, M, T] = loads(F, d, s, c, t)
% Sollicitacions de la seccio
% Inputs :  forca F [N], distancia d de la forca al perfil [mm],
%           cota de la seccio s [mm], centre de tallant c i gruix t [mm]
% Outputs:  tallant V [N], moment flector M [N*mm] i moment torsor T [N*mm]
V = F;
M = F*s;
T = F*(d + t/2 - c);
end

function [sig, tau] = nominal_stresses(V, M, T, b, h, t, Ix, It, point)
% Tensions nominals al punt 'a' (fibra extrema) o 'b' (eix neutre)
% Outputs:  tensio normal per flexio sig i tallant total tau [MPa]
tau_T = T*t/It;                      % tallant per torsio (igual a i b)
if point == 'a'
    sig = M*(h/2)/Ix;                % flexio
    Me  = (b-t)*t*(h-t)/2;           % moment estatic de l'area
else
    sig = 0;                         % a l'eix neutre no hi ha flexio
    Me  = b*t*(h-t)/2 + (h/2-t)*t*(h/2-t)/2;
end
tau_V = V*Me/t/Ix;                   % tallant per forca tallant
tau   = tau_V + tau_T;
end

function seqv = von_mises(sig, tau)
% Tensio equivalent de von Mises per a un estat pla (sig, tau)
seqv = sqrt(sig^2 + 3*tau^2);
end

function SFe = static_safety_factor(sig, tau, Sy)
% Factor de seguretat estatic (maxima energia de distorsio)
SFe = Sy/von_mises(sig, tau);
end
