%% Calcul de tensions: balancí
% DiAM, 2024

clear, close all
%% Dades del problema
%
l1 = 105;   % [mm] distància de la força F_1 a l'eix de rotació
l2 = 65;    % [mm] distància de la força F_2 a l'eix de rotació
lA = 80;    % [mm] distància de la força F_1  a la secció A-A
lB = 45;    % [mm] distància de la força F_1 al pla mig del balancí en la direcció $x$

h = 22;     % [mm] cantell de la secció a A-A
b = 10;     % [mm] gruix de la secció a A-A

F1 = 1000;  % [N]  força aplicada F1

%% *Càlcul de les sol.licitacions sobre la secció A-A*
% *Situació de la secció A-A* 
% 
% Les sol·licitacions a la secció d'estudi són:

N =   0;    % N, forca axial a la seccio
V_y = F1;   % N, forca tallant en l'eix y a la seccio
V_z = 0;    % N, forca tallant en l'eix z a la seccio
M_x = lA*F1;% N·mm, moment creat en l'eix z
M_y = 0;    % N·mm, moment creat en l'eix y
T = lB*F1;  % N·mm, torsor creat en l'eix z

%%    
% Càlcul de les propietats de la secció A-A
%
I = (1/12)*(h^3*b);  % Moment de segon ordre de la seccio circular
[eta1,eta2,eta3] = f_coef_torsio_rec(h/b);

%% Estudi dels punts crítics
% S'identifiquen tres punts crítics a la secció: 'a', 'b' i 'c'. El punt 'a' 
% és el superior de la secció, el punt 'c' està situat a la línia neutra i 
% el punt 'c' és el central a la part inferior.
% 
% 
%% Tensions en el punt 'a'
% Rep un esforç normal causat pel moment.

sig_M = M_x*(h/2)/I;    % Tensio normal del moment flector en el punt a
tau_T = 0;              % Tensio tallant del torsor en el punt a

[sig_A(1),sig_B(1),tau_max(1),sig_VM(1)]=f_invariants2D(sig_M,0,tau_T);
%% Tensions en el punt 'b'
% Rep un esforç normal causat pel moment, i un esforç tallant resultant del 
% moment torsor.

sig_M = -M_x*(h/2)/I;       % Tensio normal del moment flector en el punt b
tau_T = eta1/eta2*T/(b^2*h);% Tensio tallant del torsor en el punt b

[sig_A(2),sig_B(2),tau_max(2),sig_VM(2)]=f_invariants2D(sig_M,0,tau_T);
%% Tensions en el punt 'c'
% Està situat en el punt de màxim esforç causat pel moment torsor. 

sig_M = 0;                   % Tensio normal del moment flector en el punt b
tau_T = 1/eta2*T/(b^2*h);    % Tensio tallant del torsor en el punt b

[sig_A(3),sig_B(3),tau_max(3),sig_VM(3)]=f_invariants2D(sig_M,0,tau_T);
%% Cercles de Mohr en cada punt (del primer i últim gruix calculat)

desc_punt =['a','b','c'];
for i_punt = [1:3]
    f_cercle_Mohr2D(sig_A(i_punt),sig_B(i_punt),desc_punt(i_punt));
end

%% Definició de funcions

function [sig_A,sig_B,tau_max,sig_VM] =f_invariants2D (sigX, sigY, tauXY)
%f_invariants2D Aquesta funció calcula 4 invariants d'un estat de tensió
%               plana (bidimensional):
%               sig_A  : tensió principal màxima en el pla
%               sig_B  : tensió principal mínima en el pla
%               tau_max: tallant màxim
%               sig_VM : tensió equivalent de von Mises
% Els resultats surten ordenats a la variable 'invariants'

sig = sort(eig([sigX     tauXY;
                tauXY    sigY ]),'descend'); 
sig_A = sig(1); % Tensio principal A
sig_B = sig(2); % Tensio principal B
tau_max= 0.5*(sig(1)-sig(2)); % Tallant maxim
sig_VM = sqrt(((sig(1))^2)+((sig(2))^2)-(sig(1)*sig(2))); % Tensio Von Mises
end

function f_cercle_Mohr2D(sigma_A,sigma_B,desc_punt)
% Representa el cercle de Mohr d'un estat de tensió plana (bidimensional)
% inputs:   sigma_A  : tensió principal màxima en el pla
%           sigma_B  : tensió principal mínima en el pla
%           desc_punt, etiqueta del punt representat

H=figure;hold on;grid on;axis equal;set(gca,'FontSize',18);
theta = 0 : 0.05 : 2*pi;
%cercle sigmaI-sigmaII
plot((sigma_A-sigma_B)/2 * cos(theta) + (sigma_A+sigma_B)/2,...
    (sigma_A-sigma_B)/2 * sin(theta),'Color',[0, 0.4470, 0.7410],'LineWidth',1.6);
% dibuixa eixos
line([0 0], ylim,'Color','k','LineWidth',1.2);  %y-axis
line(xlim, [0 0],'Color','k','LineWidth',1.2);  %x-axis
% dibuixa punts
plot([sigma_A sigma_B],[0 0],'o','MarkerEdgeColor',[0, 0.4470, 0.7410],...
        'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerSize',12);
% Escriu anotacions
y_lim = ylim;gap = y_lim(2)/10;
text(sigma_A,+gap,'\sigma_{A}','FontSize',24);
text(sigma_B,+gap,'\sigma_{B}','FontSize',24);
xlabel('\sigma [MPa]');ylabel('-\tau [MPa]');
title(['Cercle de Mohr del punt ',desc_punt]); 
%saveas(H,['.\figures\BAL_punt_',desc_punt,'.svg'])
end

function [eta1,eta2,eta3] = f_coef_torsio_rec(h_b)
%% h_b data
h_b_ =[1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 100];

%% eta values
eta_1 = [1    ,0.858,0.796,0.775,0.753,0.743,0.743,0.743,0.743,0.743,0.743,0.743,0.743];
eta_2 = [0.208,0.238,0.256,0.269,0.278,0.290,0.298,0.303,0.307,0.310,0.312,0.314,0.333];
eta_3 = [0.141,0.196,0.229,0.249,0.263,0.281,0.291,0.299,0.303,0.307,0.310,0.313,0.333];

eta1 = interp1(h_b_, eta_1, min(h_b,100));
eta2 = interp1(h_b_, eta_2, min(h_b,100));
eta3 = interp1(h_b_, eta_3, min(h_b,100));
end