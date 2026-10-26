%% Suport amb eix fix
% Autor: J.A. Mayugo, 2019-2025
clear, close all

%% Dades problema
%
% Geometria
l = 60;  % mm , distancia canvi seccio eix-politja
R = 100; % mm , radi politja
D = 25;  % mm , diametre D eix politja
d = 20;  % mm , diametre d eix politja
%r_ = 4;   % mm , R4, radi acord entre D i d
r_ =[1:0.1:5]; % mm , radi acord entre D i d

%%
%
%  Forces F1 i F2 valen
F1 = 100; % N , tensio corretja 1
F2 = 400; % N , tensio corretja 2

%%
%
%  Les propietats del material son
Sut = 750; % MPa , resistencia ultima
Sy  = 550; % MPa , limit elastic

%% Obtencio de les sol.licitacions i de les tensions
%
i=0;
for r = r_
    i=i+1;
    %% Seccio A, punt a
    %
    % Calcul de sol.licitacions
    M=(l-r)*(F1+F2);
    T=R*(F2-F1);
    V=F1+F2;
    %%
    %
    % Calcul de factors de concentracio teorics, i amb sensibilitat entalla
    K_t  = K_t_rodo_flexio(D/d,r/d,1);
    K_ts = K_t_rodo_torsio(D/d,r/d);
    
    k_t_(i)  = K_t;
    k_ts_(i) = K_ts;
    
    K_t_g  = K_t_rodo_flexio_grafica(D/d,r/d);
    K_ts_g = K_t_rodo_torsio_grafica(D/d,r/d);
    
    k_t_g_(i)  = K_t_g;
    k_ts_g_(i) = K_ts_g;
    
    q  = q_sensibilitat_entalla(r, Sut, 0);
    qs = q_sensibilitat_entalla(r, Sut, 1);
    
    q_(i) =q;
    qs_(i)=qs;
    
    K_f  = 1 + q*(K_t -1);
    K_fs = 1 + qs*(K_ts -1);
    
    K_f_g  = 1 + q*(K_t_g -1);
    K_fs_g = 1 + qs*(K_ts_g -1);
    %%
    %
    % Calcul de la tensio nomal i tallant
    sigma_0=(32*M)/(pi*d^3);  % MPa
    tau_0  =(16*T)/(pi*d^3);  % MPa
    %%
    %
    % Calcul tensio von Mises en els casos 5a, 5b i 5c
    sig_VM5a(i) = ((sigma_0*K_t)^2 + 3* (tau_0*K_ts)^2)^0.5;
    sig_VM5ag(i)= ((sigma_0*K_t_g)^2 + 3* (tau_0*K_ts_g)^2)^0.5;
    
    sig_VM5b(i) = ((sigma_0*K_f)^2 + 3* (tau_0*K_fs)^2)^0.5;
    sig_VM5bg(i)= ((sigma_0*K_f_g)^2 + 3* (tau_0*K_fs_g)^2)^0.5;
    
    sig_VM5c(i) = ((sigma_0)^2 + 3* (tau_0)^2)^0.5;
end



%% Resultats
%
% En la següent taula es mostren els valors obtinguts
%
fprintf('R4     k_t   k_ts    q    q_s      Tensió amb k_t  Tensió amb k_f  Tensió amb q=0\n');
fprintf('[mm]                               [MPa]           [MPa]           [MPa]\n');
fprintf('----------------------------------------------------------------------------------\n');
for i = 1:length(r_)
    fprintf('%3.1f    %3.2f  %3.2f    %3.2f  %3.2f    %3.1f            %3.1f            %3.1f\n',...
        r_(i), k_t_(i), k_ts_(i), q_(i), qs_(i), sig_VM5a(i), sig_VM5b(i), sig_VM5c(i) );
end
%
% Diagrames de les tensions en funcio del radi 'R4'
%
figure() % primer figura amb dos gràfics

subplot(2,1,1) % --- Subplot 1: k_t i k_ts ---
hold on; % Serveix per poder dibuixar en el mateix grafic varies corbes
grid on; % Comanda per activar la malla dels eixos
plot(r_, k_t_, '-b', 'LineWidth', 1.5);
plot(r_, k_t_g_, '--b', 'LineWidth', 1.5);
plot(r_, k_ts_, '-r', 'LineWidth', 1.5);
plot(r_, k_ts_g_, '--r', 'LineWidth', 1.5);
xlabel('r [mm]', 'FontSize', 16); ylabel('k', 'FontSize', 16);
legend({'k_t per fórmula','k_t per gràfic','k_{ts} per fórmula','k_{ts} per gràfic'}, 'FontSize', 14, 'Location', 'best');
title('k_t i k_{ts}', 'FontSize', 16);
set(gca, 'FontSize', 14); % augmenta mida dels nombres dels eixos

subplot(2,1,2) % --- Subplot 2: q i q_s ---
hold on; % Serveix per poder dibuixar en el mateix grafic varies corbes
grid on; % Comanda per activar la malla dels eixos
plot(r_, q_, '-b', 'LineWidth', 1.5);
plot(r_, qs_, '-r', 'LineWidth', 1.5);
xlabel('r [mm]', 'FontSize', 16); ylabel('q', 'FontSize', 16);
legend({'q','q_{s}'}, 'FontSize', 14, 'Location', 'best');
title('q i q_s', 'FontSize', 16);
set(gca, 'FontSize', 14); % augmenta mida dels nombres dels eixos

set(gcf, 'Position', [100 100 700 800]); % augmenta la mida de la figura

figure()
hold on; % Serveix per poder dibuixar en el mateix grafic varies corbes
grid on; % Comanda per activar la malla dels eixos
plot(r_,sig_VM5a,'-b',  'LineWidth', 1.5);  
plot(r_,sig_VM5ag,'--b','LineWidth', 1.5);  
plot(r_,sig_VM5b,'-r',  'LineWidth', 1.5);  
plot(r_,sig_VM5bg,'--r','LineWidth', 1.5); 
plot(r_,sig_VM5c,'-k',  'LineWidth', 1.8);  
xlabel('r [mm]', 'FontSize', 16);ylabel('\sigma [MPa]', 'FontSize', 16);
legend('amb K_t i K_{ts} per formula','amb K_t i K_{ts} per grafic',...
       'amb K_f i K_{fs} per formula','amb K_f i K_{fs} per grafic',...
       'amb q=0','FontSize', 14,'Location','northeast');
title('Tensio von Mises del punt a','FontSize', 16);
set(gca, 'FontSize', 14); % augmenta mida dels nombres dels eixos

%% Funcions 

function [ K_t] = K_t_rodo_flexio( D_d , r_d, i_flag )
%Calcula la K_t d'un rodo sotmes a flexio amb un canvi de seccio de D a d:
%   D = diametre gran; d= diametre petit; r=radi de l'entalla
%   D_d = D/d ; r_d = r/D
%   i_flag = 1, si es vol K_t per tensio normal
%   i_flag <> 1, si es vol K_t per von Mises

if i_flag==1 
    K_t = 0.632+0.377*D_d^(-4.4)+r_d^(-0.5)*...
      sqrt((-0.14-0.363*D_d^2+0.503*D_d^4)/(1-2.39*D_d^2+3.368*D_d^4));
else 
    K_t = 0.622+0.38*D_d^(-4.3)+r_d^(-0.5)*...
      sqrt((-0.322-0.277*D_d^2+0.599*D_d^4)/(1-2.55*D_d^2+5.27*D_d^4));
end
end

function [ K_t] = K_t_rodo_torsio( D_d , r_d)
%Calcula la K_t d'un rodó sotmés a torsió amb un canvi de secció de D a d:
%   D = diàmetre gran; d= diàmetre petit; r=radi de l'entalla
%   D_d = D/d ; r_d = r/D
 
    K_t = 0.78+0.2*D_d^(-10)+r_d^(-0.46)*...
      sqrt((0.002-0.125*D_d^2+0.123*D_d^4)/(1-2.75*D_d^2+2.55*D_d^4));

end

function [K_t] = K_t_rodo_flexio_grafica(D_d,r_d)
    r_d = max(min(r_d,0.3),0.025);
    D_d = max(min(D_d,3),1.02);

    %% r_d data
    r_d_vec =[	0.025	,	0.050	,	0.075	,	0.100	,	0.125	,	0.150	,	0.175	,	0.200	,	0.225	,	0.250	,	0.275	,	0.300	];

    %% r_d data flexi?
    D_d_vec_f =[	1.02	,	1.05	,	1.1     ,	1.5 	,	3	];

    %% K_t values flexi?
    data_f = [1.913	,	1.667	,	1.537	,	1.450	,	1.388	,	1.341	,	1.304	,	1.268	,	1.239	,	1.228	,	1.209	,	1.198	;
            2.140	,	1.804	,	1.645	,	1.544	,	1.471	,	1.410	,	1.359	,	1.326	,	1.297	,	1.275	,	1.257	,	1.238	;
            2.256	,	1.884	,	1.699	,	1.594	,	1.515	,	1.464	,	1.413	,	1.366	,	1.366	,	1.366	,	1.366	,	1.366	;
            2.591	,	2.075	,	1.826	,	1.677	,	1.572	,	1.507	,	1.460	,	1.416	,	1.400	,	1.400	,	1.400	,	1.400	;
            2.852	,	2.317	,	2.002	,	1.803	,	1.677	,	1.587	,	1.510	,	1.452	,	1.402	,	1.362	,	1.325	,	1.307	];
    data_f_t = data_f';

    K_t = interp2(D_d_vec_f,r_d_vec,data_f_t,D_d,r_d);
end

function [K_ts] = K_t_rodo_torsio_grafica(D_d,r_d)
    r_d = max(min(r_d,0.3),0.025);
    D_d = max(min(D_d,2),1.09);
    
    %% r_d data
    r_d_vec =[	0.025	,	0.050	,	0.075	,	0.100	,	0.125	,	0.150	,	0.175	,	0.200	,	0.225	,	0.250	,	0.275	,	0.300	];

    %% r_d data torsi?
    D_d_vec_t =[	1.09	,	1.20	,	1.33     ,	2 	];

    %% K_ts values torsi?
    data_t = [1.489	,	1.290	,	1.209	,	1.165	,	1.140	,	1.125	,	1.114	,	1.106	,	1.100	,	1.100	,	1.100	,   1.100;
            1.915	,	1.592	,	1.434	,	1.345	,	1.283	,	1.239	,	1.195	,	1.162	,	1.140   ,	1.125	,	1.118	,   1.105;
            2.033	,	1.695	,	1.518	,	1.426	,	1.353	,	1.305	,	1.254	,	1.224	,	1.195	,	1.180	,	1.173	,   1.170;
            2.154	,	1.768	,	1.573	,	1.467	,	1.397	,	1.342	,	1.301	,	1.261	,	1.228	,	1.202	,	1.191	,   1.180];
    data_t_t = data_t';

    K_ts = interp2(D_d_vec_t,r_d_vec,data_t_t,D_d,r_d);
end

function [q] = q_sensibilitat_entalla(r, Sut, iflag)
% Calcula la sensibilitat a l'entalla (q) per acers d'elemets sota 
% flexio alterna o forca axial alterna, o tallant altern
%   r: radi de l'entalla, r input en [mm], formula en [in]
%                   1 in = 25.4 mm
%   Sut: resitencia ultima, Sut input en [MPa), formula en [kpsi]
%                   1 kpsi = 6.9 MPa
%   iflag = 0 flexio/axial, iflag <> 0 tallant

Sut_kpsi=Sut/6.9; % [kpsi]
r_in = r/25.4;    % [in]

if iflag == 0
    sqr_a=0.246-3.08e-3*Sut_kpsi+1.51e-5*Sut_kpsi^2-2.67e-8*Sut_kpsi^3;
else
    sqr_a=0.190-2.51e-3*Sut_kpsi+1.35e-5*Sut_kpsi^2-2.67e-8*Sut_kpsi^3;
end

q=1/(1+sqr_a/sqrt(r_in));
end

