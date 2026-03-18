%%% Disseny de una cel·la càrrega a tracció
%%% J.A. Mayugo 8/01/2022
clear, close all

%% Dades problema
F = 10e3;           % N, força nominal / rang

E = 190e3;          % MPa, elasticitat acer inoxidable
nu= 0.29;           % -  , coeficient de Poisson
Sy= 240;            % MPa, límit elàstic

Vo_Vi_sens = 2e-6;  % V/V, sensibilitat de l'equip (+-)
K = 2            ;  % -  , factor de galga

SF_d = 3;           % -    , SF estàtic mínim
sens_F_d = 50;      % N, sensibilitat dessitjada

d = 100;            % mm, diàmetre test inicial

disp([newline,'Determinació de la geometria ', ...
    'amb una sensibilitat dessitjada de ',num2str(round(sens_F_d,1)),newline ]);

%% MÈTODE de la biseccio per obtenir mides òptimes
sens_F_e = sens_F_d*2; cont = 0; sens_F_eps = sens_F_d/1000;
d_min = d/100; d_max = d*10; % mm,  diàmetre inicial mínim i màxim
% Generació del bucle 'while' per a calcular les dimensions màximes per
% aconseguir la sensibilitat dessitjada
while (abs(sens_F_d-sens_F_e) > sens_F_eps)        % repeteix mentre no troba solucio
    cont = cont +1;
    d  = (d_min+d_max)/2; % punt mig
    [sens_F_e,SF_e,factible] = sens_F_trial(d,F,E,nu,Sy,SF_d,Vo_Vi_sens,K);
    if sens_F_d >= sens_F_e
        d_min=d;
    else
        d_max=d;
    end
    disp (['# iter.: ',num2str(cont),' el diàmetre és de ',num2str(d)...
        'mm. La sensibilitat obtinguda és ',num2str(round(sens_F_e,1))...
        '.El SF obtingut és ',num2str(round(SF_e,2)),' (',factible,')']);
end
%% Resultat arrodonit a geometria de 0.1 mm
inc=0.1;d=floor(d/inc)*inc;
[sens_F_e,SF_e,factible] = sens_F_trial(d,F,E,nu,Sy,SF_d,Vo_Vi_sens,K);

disp([newline,'Resultat a la iteració: ',num2str(cont),'.',...
    newline,'Amb força ',num2str(F),...
    'N, el nominal nominal és de d=',num2str(ceil(d/inc)*inc),'mm.',...
    newline,'El SF obtingut és ',num2str(round(SF_e,2)),...
    ' (',factible,')',newline]);

%% Resultats
disp(['Mida del disseny final']);
disp(['  d = ',num2str(d),' mm',newline]);
disp(['Factor de seguretat']);
disp(['  SF = ',num2str(min(SF_e)),' (',factible,')',newline]);
disp(['Valor de la sensibilitat i resolució']);
disp(['  Sensibilitat_F = ',num2str(round(sens_F_e,1)),' N']);
disp(['  Resolució      = ',num2str(floor(F/sens_F_e)),' mesures']);

%% Funcions emprades

function [sens_F_e,SF_e,factible] = sens_F_trial(d,F,E,nu,Sy,SF_d,Vo_Vi_sens,K)
% Càlcul de la sensibilitat del disseny
% Input:    d (mm, diàmetre), F (N, força), E i nu (prop. elàstiques), 
%           Sy i SF_d (prop. resistència), Vo_Vi_sens i K (extensiometria)
% Output:   sens_F_e (N, sensibilitat del disseny)
%           SF_e (factor seguretat)
%           factible (etiqueta si el disseny és o no és factible 

A=pi*d^2/4;                     % àrea
sig =  F/A;                     % tensió longitudinal
eps = [1/E -nu/E -nu/E;-nu/E 1/E -nu/E;-nu/E -nu/E 1/E]*[sig 0 0]';
eps_wheastone = 2 * (eps(1));   % 2 galgues actives en direcció longitudinal
eps_sens = 4/K*Vo_Vi_sens;      % mm/mm, sensibilitat de l'equip (+-)
sens_F_e = eps_sens/eps_wheastone*F;

SF_e = Sy/sig;                  % factor de seguretat
if SF_e >= SF_d
    factible = 'disseny factible';
else
    factible = 'disseny NO factible';
end
end
