%% Anàlisi estàtica de la manovella
% Factor de seguretat estàtic: secció 1 (punt a) i secció 2 (punts a i b)
% J.A. Mayugo 9/01/2020 (revisat 2026)
clear; close all; clc

%% Dades del problema
% Força i posició
F  = 150;     % N, força aplicada
d  = 20;      % mm, distància entre la força i el perfil en U
% Propietats del material (alumini)
Sut = 150;    % MPa, resistència última
Sy  = 70;     % MPa, límit elàstic
% Geometria de les seccions
h1 = (40-30)/200*180 + 30;   % mm, alçada del cantell (secció 1)
b1 = 25;      % mm, amplada (secció 1)
h2 = (40-30)/200*15 + 30;    % mm, alçada del cantell (secció 2)
b2 = 20;      % mm, amplada (secció 2)
t  = 4;       % mm, gruix (seccions 1 i 2)

%% Propietats de les seccions
% Secció 1
Ix1 = b1*h1^3/12 - (b1-t)*(h1-2*t)^3/12;          % mm^4, moment d'inèrcia
It1 = 1.1*(2*0.300*b1*t^3 + 0.306*(h1-2*t)*t^3);  % mm^4, constant de torsió
c1  = (b1-t/2)^2*(h1-t)^2*t/4/Ix1;                % mm, centre de tallant
% Secció 2
Ix2 = b2*h2^3/12 - (b2-t)*(h2-2*t)^3/12;
It2 = 1.1*(2*0.291*b2*t^3 + 0.306*(h2-2*t)*t^3);
c2  = (b2-t/2)^2*(h2-t)^2*t/4/Ix2;

%% Sol·licitacions (força, moments flector i torsor)
V  = F;                   % N, tallant (seccions 1 i 2)
M1 = F*180;               % N*mm, moment flector secció 1
T1 = F*(d + t/2 - c1);    % N*mm, moment torsor secció 1
M2 = F*15;                % N*mm, moment flector secció 2
T2 = F*(d + t/2 - c2);    % N*mm, moment torsor secció 2

%% Tensions nominals
% Secció 1, punt a
sigM1 = M1*h1/2/Ix1;                % MPa, tensió normal per flexió
tauT1 = T1*t/It1;                   % MPa, tallant per torsió
Me1   = (b1-t)*t*(h1-t)/2;          % mm^3, moment estàtic
tauV1 = V*Me1/t/Ix1;                % MPa, tallant per força tallant
tau1  = tauV1 + tauT1;              % MPa, tallant total

% Secció 2, punt a
sigM2 = M2*h2/2/Ix2;
tauT2 = T2*t/It2;
Me2   = (b2-t)*t*(h2-t)/2;
tauV2 = V*Me2/t/Ix2;
tau2  = tauV2 + tauT2;

% Secció 2, punt b (a l'eix neutre, sense tensió normal)
Me2b   = b2*t*(h2-t)/2 + (h2/2-t)*t*(h2/2-t)/2;
tauV2b = V*Me2b/t/Ix2;
tau2b  = tauV2b + tauT2;

%% Factor de seguretat estàtic (von Mises)
eqv1  = sqrt(sigM1^2 + 3*tau1^2);   % MPa
eqv2  = sqrt(sigM2^2 + 3*tau2^2);
eqv2b = sqrt(0^2 + 3*tau2b^2);

SFe1  = Sy/eqv1;
SFe2  = Sy/eqv2;
SFe2b = Sy/eqv2b;

%% Resultats
fprintf('Secció 1, punt a : sigma = %5.1f MPa, tau = %5.1f MPa, sigma_eq = %5.1f MPa, SFe = %4.2f\n', sigM1, tau1,  eqv1,  SFe1);
fprintf('Secció 2, punt a : sigma = %5.1f MPa, tau = %5.1f MPa, sigma_eq = %5.1f MPa, SFe = %4.2f\n', sigM2, tau2,  eqv2,  SFe2);
fprintf('Secció 2, punt b : sigma = %5.1f MPa, tau = %5.1f MPa, sigma_eq = %5.1f MPa, SFe = %4.2f\n', 0,     tau2b, eqv2b, SFe2b);
