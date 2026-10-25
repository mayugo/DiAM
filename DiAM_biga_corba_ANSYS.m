%% Solució cas d'estudi biga corba: analític vs numèric
% J.A. Mayugo, 2020

clear all; close all;
%% Dades problema
% Dades d'ANSYS

path = [0 0.83333 1.6667 2.5 3.3333 4.1667 5 5.8333 6.6667 7.5 8.3333 9.1667 ...
   10 10.833 11.667 12.5 13.333 14.167 15 15.833 16.667 17.5 18.333 19.167 ...
   20 20.833 21.667 22.5 23.333 24.167 25 25.833 26.667 27.5 28.333 29.167 ...
   30 30.833 31.667 32.5 33.333 34.167 35 35.833 36.667 37.5 38.333 39.167 ...
   40];

sigmaX_A = [
   -34.6 -33.065 -31.529 -29.993 -28.444 -26.865 -25.285 -23.706 -22.184 ...
   -20.698 -19.212 -17.726 -16.311 -14.901 -13.491 -12.091 -10.729 -9.3662 ...
   -8.0036 -6.6562 -5.3216 -3.9871 -2.6526 -1.3257 -0.00020224 1.3253 2.6522 ...
   3.9869 5.3216 6.6562 8.0037 9.3661 10.729 12.091 13.49 14.901 16.311 ...
   17.727 19.201 20.676 22.15 23.661 25.23 26.799 28.368 30.06 31.806 33.553 ...
   35.299];

sigmaX_B = [ 
   -28.012 -27.18 -26.349 -25.518 -24.671 -23.786 -22.901 -22.016 -21.09 ...
   -20.138 -19.186 -18.234 -17.211 -16.182 -15.153 -14.104 -12.987 -11.871 ...
   -10.754 -9.5853 -8.3708 -7.1562 -5.9417 -4.6336 -3.3084 -1.9832 -0.63851 ...
   0.8135 2.2655 3.7175 5.2381 6.8387 8.4393 10.04 11.779 13.559 15.339 ...
   17.137 19.139 21.141 23.143 25.252 27.53 29.808 32.086 34.592 37.199 ...
   39.806 42.414];

sigmaX_C = [ 
   -30.925 -29.83 -28.736 -27.641 -26.525 -25.358 -24.192 -23.026 -21.84 ...
   -20.643 -19.446 -18.248 -17.022 -15.793 -14.565 -13.327 -12.059 -10.791 ...
   -9.5228 -8.2305 -6.9174 -5.6042 -4.2911 -2.9333 -1.5673 -0.20141 1.174 ...
   2.6019 4.0297 5.4575 6.919 8.4198 9.9206 11.421 12.989 14.576 16.164 ...
   17.758 19.445 21.131 22.818 24.551 26.36 28.169 29.978 31.928 33.941 ...
   35.954 37.967];

sigmaX_D = [
   -33.977 -32.532 -31.088 -29.644 -28.203 -26.771 -25.338 -23.906 -22.481 ...
   -21.061 -19.641 -18.221 -16.811 -15.402 -13.993 -12.585 -11.183 -9.7802 ...
   -8.3777 -6.9765 -5.5764 -4.1762 -2.7761 -1.3745 0.027375 1.4292 2.8317 ...
   4.2376 5.6435 7.0494 8.4569 9.8662 11.276 12.685 14.094 15.502 16.911 ...
   18.319 19.72 21.121 22.522 23.917 25.303 26.688 28.074 29.446 30.811 ...
   32.177 33.543];
%% 
% Estat de càrregues, i geometria

F = 4500;   % N, forces verticals
d = 40;     % mm, distància entre forces

h = 40;     % mm, cantell biga
b = 20;     % mm, gruix biga

ro = 70;    % mm, radi exterior biga corba
%% Propietats de la seccio
% biga recta

A = h*b;            % àrea secció
I = 1/12 * h^3*b;   % moment segon ordre secció
%% 
% biga corba

ri = ro-h;          % mm, radi interior biga corba
rc = (ro+ri)/2;     % mm, radi centroide

rn = h /log(ro/ri); % mm, radi línia neutre
e = rc-rn;
%% Càlculs

M = F*d;            % N·mm, moment flector

sigmaX_recta = [];
sigmaX_corba = [];
for path_  = path
    y_recta = path_ - h/2;
    y_corba = path_ - (h/2 + e);

    sigmaX_recta = [sigmaX_recta   M*y_recta/I];
    sigmaX_corba = [sigmaX_corba   M*y_corba/(A*e*(rn-y_corba))];
end
%% Representació gràfica de la solució

figure
x0=10;
y0=10;
width=1000;
height=800;
set(gcf,'position',[x0,y0,width,height])

subplot(2,2,1)
hold on, grid on, set(gca,'fontsize', 20)
p1=plot(sigmaX_A,path,'x','MarkerSize',6,'linewidth',1.8);
p2=plot(sigmaX_recta,path,'linewidth',1.8);
p3=plot(sigmaX_corba,path,'linewidth',1.8);
plot([0 0],[path(1) path(end)],'k')
xlabel('$\sigma_L$ [MPa]','Interpreter','latex','fontsize', 24)
ylabel('$y$ [mm]','Interpreter','latex','fontsize', 24)
legend([p1,p2,p3],{'Secció A numèric','biga recta','biga corba'},'Location', 'best')

subplot(2,2,2)
hold on, grid on, set(gca,'fontsize', 20)
p1=plot(sigmaX_B,path,'x','MarkerSize',6,'linewidth',1.8);
p2=plot(sigmaX_recta,path,'linewidth',1.8);
p3=plot(sigmaX_corba,path,'linewidth',1.8);
plot([0 0],[path(1) path(end)],'k')
xlabel('$\sigma_L$ [MPa]','Interpreter','latex','fontsize', 24)
ylabel('$y$ [mm]','Interpreter','latex','fontsize', 24)
legend([p1,p2,p3],{'Secció B numèric','biga recta','biga corba'},'Location', 'best')

subplot(2,2,3)
hold on, grid on, set(gca,'fontsize', 20)
p1=plot(sigmaX_C,path,'x','MarkerSize',6,'linewidth',1.8);
p2=plot(sigmaX_recta,path,'linewidth',1.8);
p3=plot(sigmaX_corba,path,'linewidth',1.8);
plot([0 0],[path(1) path(end)],'k')
xlabel('$\sigma_L$ [MPa]','Interpreter','latex','fontsize', 24)
ylabel('$y$ [mm]','Interpreter','latex','fontsize', 24)
legend([p1,p2,p3],{'Secció C numèric','biga recta','biga corba'},'Location', 'best')

subplot(2,2,4)
hold on, grid on, set(gca,'fontsize', 20)
p1=plot(sigmaX_D,path,'x','MarkerSize',6,'linewidth',1.8);
p2=plot(sigmaX_recta,path,'linewidth',1.8);
p3=plot(sigmaX_corba,path,'linewidth',1.8);
plot([0 0],[path(1) path(end)],'k')
xlabel('$\sigma_L$ [MPa]','Interpreter','latex','fontsize', 24)
ylabel('$y$ [mm]','Interpreter','latex','fontsize', 24)
legend([p1,p2,p3],{'Secció D numèric','biga recta','biga corba'},'Location', 'best')