%%% Disseny anell de mesura de forces 'proving ring'
%%% J.A. Mayugo 10/01/2021
clear, close all

%% Dades problema
P = 400;   % N, força nominal / rang

E = 190e3; % MPa, elasticitat acer inoxidable
nu= 0.29;  % -  , coeficient de Poisson
Sy= 600;   % MPa, límit elàstic

%%
% Geometria 'test'
Dm= 20;  % mm
h = 2;    % mm
b = 6;    % mm

eps_sens = 1e-6; % mm/mm, sensibilitat de l'equip (+-)
SF_d = 3;         % -    , SF estàtic desitjat

disp([newline,'Apartat 3) determinació de la geometria ', ...
    'amb un SF dessitjat de ',num2str(round(SF_d,2)),newline ]);

%% METODE de la biseccio per obtenir mides òptimes
SF_e = 0; cont = 0; SF_eps = 0.05;
h_min = h/10; h_max = h*10; % mm,  gruix perfil mínim i màxim
% Generació del bucle 'while' per a calcular les dimensions mínimes per aconseguir el SF desitjat
while (abs(SF_d-SF_e)>SF_eps)        % repeteix mentre no troba solucio
    cont = cont +1;
    h  = (h_min+h_max)/2; % punt mig
    Dm = h*10;
    b  = h*3;
    [SF_e,sig_o,sig_i] = SF_trial(Dm,h,b,P,Sy);
    if SF_d>=SF_e
        h_min=h;
    else
        h_max=h;
    end
    disp (['# iter.: ',num2str(cont),' el gruix es de ',num2str(h)...
        'mm. El SF obtingut és ',num2str(round(SF_e,2))]);
end
%% Resultat arrodonit a geometria de 0.01 mm
inc=0.01;h=ceil(h/inc)*inc;Dm=h*10;b=h*3;
[SF_e,sig_o,sig_i] = SF_trial(Dm,h,b,P,Sy);

disp([newline,'Resultat a la iteració: ',num2str(cont),'.',...
      newline,'Amb força ',num2str(P),...
      'N, el gruix nominal és de h=',num2str(ceil(h/inc)*inc),'mm.',...
      newline,'El SF obtingut és ',num2str(round(SF_e,2)),newline]);

%% Valoració de la sensibilitat del disseny òptim
eps_i = [1/E -nu/E -nu/E;-nu/E 1/E -nu/E;-nu/E -nu/E 1/E]*[sig_i 0 0]';
eps_o = [1/E -nu/E -nu/E;-nu/E 1/E -nu/E;-nu/E -nu/E 1/E]*[0 0 sig_o]';
eps_wheastone = 2 * (eps_i(1) - eps_o(3));
F_sens = eps_sens/eps_wheastone*P;

%% Resultats
disp(['Mides del disseny final']);
disp(['  Dm = ',num2str(Dm),' mm']);
disp(['  h  = ',num2str(h),' mm']);
disp(['  b  = ',num2str(b),' mm',newline]);
disp(['Factor de seguretat']);
disp(['  SF = ',num2str(min(SF_e)),newline]);
disp(['Valor de la sensibilitat']);
disp(['  Sensibilitat_F = ',num2str(F_sens),' N']);
disp(['  Resolució      = ',num2str(round(P/F_sens)),' mesures']);

%% Estat de tensió i deformació en els punts de col·locació galgues

disp([newline,'Apartat 5) Cercles de Mohr',newline ]);

cercle_Mohr(sig_i,0,0,'Estat de tensions galga interior','sigma',[0, 0.4470, 0.7410])
cercle_Mohr(eps_i(1),eps_i(2),eps_i(3),'Estat de deformacions galga interior','epsilon',[0, 0.4470, 0.7410])
cercle_Mohr(0,0,sig_o,'Estat de tensions galga exterior','sigma',[0, 0.4470, 0.7410])
cercle_Mohr(eps_o(1),eps_o(2),eps_o(3),'Estat de deformacions galga exterior','epsilon',[0, 0.4470, 0.7410])

%% Funcions emprades

function [SFe,sig_o,sig_i]=SF_trial(Dm,h,b,P,Sy)
%% Analisi del disseny
ro = Dm/2+h/2;      % radi exterior
ri = Dm/2-h/2;      % radi interior
r_n= h/log(ro/ri);  % línia neutra
r_c=Dm/2;           % línia centroide

Mo = P*r_n*(1/2-1/pi); % moment flector

sig_mo = Mo *(r_n-ro)/(b*h*ro*(r_c-r_n)) ;
sig_mi = Mo *(r_n-ri)/(b*h*ri*(r_c-r_n)) ;
sig_a = (P/2)/(b*h);

sig_o = sig_mo - sig_a ;    % outer point
sig_i = sig_mi + sig_a ;    % inner point

SFe = max([abs(Sy/sig_i) abs(Sy/sig_o)]);
end

function sigma= cercle_Mohr(sigma_I,sigma_II,sigma_III,title_,labels,color_,xylim)
%% Representa el cercle de Mohr d'un estat de tensió plana
% inputs:   sigma_I,    tensió principal màxima
%           sigma_II,   tensió principal mitja
%           sigma_III,  tensió principal mínima
%           title_,     títol del gràfic
%           labels,     de tensions 'sigma' o de deformacions 'epsilon'
%           color,      color
%           xylim,      límits del gràfic predefinits
% output:   sigma,  tensions principals

if isequal(labels, 'epsilon')
    sigma = sort([sigma_I,sigma_II,sigma_III])*1e6;
    pabs = '$\varepsilon \times 10^{-6}$';
    pord = '$\gamma/2 \times 10^{-6}$';
    pI =    '$\varepsilon_{I}$';
    pII=    '$\varepsilon_{II}$';
    pIII=   '$\varepsilon_{III}$';
else
    sigma = sort([sigma_I,sigma_II,sigma_III]);
    pabs = '$\sigma$ [MPa]';
    pord = '$\tau$ [MPa]';
    pI =    '$\sigma_{I}$';
    pII=    '$\sigma_{II}$';
    pIII=   '$\sigma_{III}$';  
end
sigma_I = sigma(3); sigma_II = sigma(2); sigma_III = sigma(1);

%% Cercle de Mohr 
f0=figure;hold on;grid on,axis equal;set(gca,'FontSize',18)
theta = 0 : 0.05 : 2*pi;

%cercle sigmaI-sigmaII-sigmaIII
plot((sigma_I-sigma_II)/2 * cos(theta) + (sigma_I+sigma_II)/2,...
    (sigma_I-sigma_II)/2 * sin(theta),'Color',color_,'LineWidth',1.6);
plot((sigma_I-sigma_III)/2 * cos(theta) + (sigma_I+sigma_III)/2,...
    (sigma_I-sigma_III)/2 * sin(theta),'Color',color_,'LineWidth',1.6);
plot((sigma_II-sigma_III)/2 * cos(theta) + (sigma_II+sigma_III)/2,...
    (sigma_II-sigma_III)/2 * sin(theta),'Color',color_,'LineWidth',1.6);

if ~exist('xylim','var')
     % third parameter does not exist, so default it to something
      x_lim = xlim;y_lim = ylim;
      xylim = [x_lim y_lim];
else
      xlim(xylim(1:2));ylim(xylim(3:4));
end
title(title_);
% dibuixa eixos
line([0 0], [xylim(3:4)],'Color','k','LineWidth',1.2);  %y-axis
line([xylim(1:2)], [0 0],'Color','k','LineWidth',1.2);  %x-axis
xlabel(pabs,'interpreter','latex');ylabel(pord,'interpreter','latex');  
% dibuixa punts
plot([sigma_I sigma_II sigma_III],[0 0 0],'o','MarkerEdgeColor',color_,...
     'MarkerFaceColor',color_,'MarkerSize',12)
% Escriu anotacions
%y_lim = ylim;
gap = xylim(4)/12;
if sigma_I == sigma_II 
    signe = 1;
else
    signe = 0;
end
text(sigma_I+gap*0.3,   +gap,pI,    'FontSize',24,'interpreter','latex');
text(sigma_II+gap*0.3 - signe*gap*3.3,  +gap,pII,   'FontSize',24,'interpreter','latex');
text(sigma_III-gap*3.2, +gap,pIII,  'FontSize',24,'interpreter','latex');
saveas(gcf,['../figures/',title_,'.png'])
end