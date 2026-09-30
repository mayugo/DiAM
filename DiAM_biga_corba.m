%% Biga corba vs biga recta a flexió (solució analítica)
clear; close all;

%% Dades
F = 4500;           % N, forces verticals
d = 40;             % mm, distància entre forces
h = 40;             % mm, cantell biga
b = 20;             % mm, gruix biga
ro = 70;            % mm, radi exterior

M = F*d;            % N·mm, moment flector

%% Propietats de la secció
A  = h*b;           % mm^2
I  = b*h^3/12;      % mm^4
ri = ro - h;        % mm, radi interior
rc = (ro+ri)/2;     % mm, radi del centroide
rn = h/log(ro/ri);  % mm, radi de la línia neutra
e  = rc - rn;       % mm, excentricitat

%% Tensions (path: 0 = fibra exterior, h = fibra interior)
path = linspace(0, h, 200);

y_recta = path - h/2;           % des del centroide
y_corba = path - (h/2 + e);     % des de la línia neutra

sigma_recta = M*y_recta/I;
sigma_corba = M*y_corba./(A*e*(rn - y_corba));

%% Posició de les fibres neutres i tensions als extrems
path_n_recta = h/2;             % mm des de la fibra exterior
path_n_corba = h/2 + e;

ext = [1 numel(path)];          % índexs: fibra exterior, fibra interior
s_rec = sigma_recta(ext);
s_cor = sigma_corba(ext);
dif_pct = (abs(s_cor) - abs(s_rec)) ./ abs(s_rec) * 100;   % > 0: corba més tensionada

fprintf('Desplaçament fibra neutra (e) = %.2f mm\n', e);
fprintf('Fibra exterior: recta %.1f MPa | corba %.1f MPa | %+.1f %%\n', s_rec(1),   s_cor(1),   dif_pct(1));
fprintf('Fibra interior: recta %.1f MPa | corba %.1f MPa | %+.1f %%\n', s_rec(end), s_cor(end), dif_pct(end));

%% Gràfic
figure('Position',[10 10 800 600])
hold on, grid on, set(gca,'fontsize',14)

p1 = plot(sigma_recta, path, 'linewidth',1.8);
p2 = plot(sigma_corba, path, 'linewidth',1.8);
plot([0 0], [0 h], 'k')

xm = 1.9*max(abs([sigma_recta sigma_corba]));
xlim([-xm xm]), ylim([-3 h+3])

% fibres neutres
plot([-xm xm], path_n_recta*[1 1], '--', 'color', p1.Color)
plot([-xm xm], path_n_corba*[1 1], '--', 'color', p2.Color)
text(-0.98*xm, path_n_recta, sprintf('F. neutra recta (y = %.1f mm)', path_n_recta), ...
    'color', p1.Color, 'VerticalAlignment','top', 'fontsize',14)
text(-0.98*xm, path_n_corba, sprintf('F. neutra corba (y = %.1f mm, e = %.2f mm)', path_n_corba, e), ...
    'color', p2.Color, 'VerticalAlignment','bottom', 'fontsize',14)

% extrems
plot(s_rec, path(ext), 'o', 'color', p1.Color, 'MarkerFaceColor', p1.Color)
plot(s_cor, path(ext), 'o', 'color', p2.Color, 'MarkerFaceColor', p2.Color)

txt = @(i) sprintf('recta: %.1f MPa\ncorba: %.1f MPa\n\\Delta = %+.1f %%', s_rec(i), s_cor(i), dif_pct(i));
text(min(s_rec(1),s_cor(1))   - 0.03*xm, path(1),   txt(1), 'HorizontalAlignment','right', 'fontsize',14)
text(max(s_rec(2),s_cor(2))   + 0.03*xm, path(end), txt(2), 'HorizontalAlignment','left',  'fontsize',14)

xlabel('$\sigma_L$ [MPa]','Interpreter','latex','fontsize',18)
ylabel('$y$ [mm]  (0 = fibra exterior)','Interpreter','latex','fontsize',16)
legend([p1 p2], {'Biga recta','Biga corba'}, 'Location','northwest','Interpreter','latex','fontsize',16)