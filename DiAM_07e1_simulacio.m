function [t, x, y, F, X_ss, F_ss] = DiAM_07e1_simulacio(m, k, zeta, y0, omega, t_end, n_per)
% Simula el moviment vertical x(t) d'una massa m (kg) sobre una molla de
% rigidesa total k (N/m) i esmorteiment relatiu zeta, quan la base es mou
% com y = y0*sin(omega*t) (y0 en m, omega en rad/s).
% Condicions inicials de repos (x = 0, dx/dt = 0).
%
% t_end : temps de simulacio en s (per defecte, 6 constants de temps)
% n_per : punts de sortida per periode d'excitacio
%
% Retorna:
%   t, x, y : temps, desplacament del motlle i de la base
%   F       : forca dinamica sobre cada suport (N), F = (k*(x-y) + c*(dx-dy))/2
%   X_ss    : amplitud de x a regim permanent (ultims 5 periodes)
%   F_ss    : amplitud de F a regim permanent (ultims 5 periodes)

wn = sqrt(k/m);
c = 2*zeta*sqrt(k*m);
if nargin < 6 || isempty(t_end), t_end = 6/(zeta*wn); end
if nargin < 7 || isempty(n_per), n_per = 60; end
T = 2*pi/omega;
t = (0:T/n_per:t_end).';

f = @(t, z) [z(2); ...
    (-c*(z(2) - y0*omega*cos(omega*t)) - k*(z(1) - y0*sin(omega*t)))/m];
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-11);
[t, z] = ode45(f, t, [0; 0], opts);

x = z(:, 1);
v = z(:, 2);
y = y0*sin(omega*t);
yd = y0*omega*cos(omega*t);
F = (k*(x - y) + c*(v - yd))/2;
ss = t >= t(end) - 5*T;
X_ss = max(abs(x(ss)));
F_ss = max(abs(F(ss)));
end
