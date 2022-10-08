%%% Disseny anell de mesura de forces 'proving ring'
%%% J.A. Mayugo 10/01/2021
clear, close all

%% Dades problema test
P = 400;   % N

E = 190e3; % MPa
nu= 0.29;  % -
Sy= 600;   % MPa

Dm= 20  % mm 
h = 2    % mm
b = 6    % mm

eps_sens  = 1e-6; % mm/mm

%% Analisi del disseny
ro = Dm/2+h/2;       % radi exterior
ri = Dm/2-h/2;       % radi interior
r_n = h/log(ro/ri);    % linia neutra
r_c=Dm/2;           % linia centroide

Mo = P*r_n*(1/2-1/pi); % moment flector

sig_mo = Mo *(r_n-ro)/(b*h*ro*(r_c-r_n)) ;
sig_mi = Mo *(r_n-ri)/(b*h*ri*(r_c-r_n)) ;
sig_a = (P/2)/(b*h);

sig_o = sig_mo - sig_a ;    % outer point
sig_i = sig_mi + sig_a ;    % inner point

SFe = [abs(Sy/sig_i) abs(Sy/sig_o)];

F_sens = 0; %% pendent de programar

%% Resultats
disp(['Disseny']);
disp(['  Dm = ',num2str(Dm),' mm']);
disp(['  h  = ',num2str(h),' mm']);
disp(['  b  = ',num2str(b),' mm']);
disp([' ']);
disp(['Factor de seguretat']);
disp(['  SF = ',num2str(min(SFe)),' ']);
disp([' ']);
disp(['Valor de la sensibilitat']);
disp(['  Sensibilitat_F = ',num2str(F_sens),' N']);