%%% Suport amb eix fix
%%% J.A. Mayugo 2015
clear, close all

%% Dades problema
l = 60; % mm , distancia canvi seccio eix-politja
R = 100; % mm , radi politja
D = 25; % mm , diametre D eix politja
d = 20; % mm , diametre d eix politja
r = 4; % mm , R4, radi acord entre D i d

F1 = 100; % N , força corretja 1
F2 = 400; % N , força corretja 2

%% Seccio A, punt a
% Càlcul de sol·licitacions 
M=(l-r)*(F1+F2);
T=R*(F2-F1);
V=F1+F2;

% Càlcul de factors de concentració teòrics K_t
K_t = K_t_rodo_flexio(D/d,r/d,1);

% Càlcul de la tensió nomal 
sigma=(32*M)/(pi*d^3)*K_t  % MPa