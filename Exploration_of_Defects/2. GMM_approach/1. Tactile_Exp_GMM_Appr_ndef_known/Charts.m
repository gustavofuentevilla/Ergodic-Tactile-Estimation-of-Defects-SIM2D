%% Charts

close all
clearvars -except M
clc

%22,28,38,42 close def
%35,57,58,70,73,78,85,92,94,99,100 overlaped defects
%86 overlaped defects. Failure?
load("Results100/output_86.mat") 

%%

fig1h = figure(1);
layout1h = tiledlayout(fig1h, 3, 1);
nexttile
plot(t_total, X_e_total, 'LineWidth', 2)
title("Position States")
xlabel('Time [s]')
ylabel('Position [m]')
legend('$x_1$', '$x_2$')
grid on
nexttile
plot(t_total, X_e_dot_total, 'LineWidth', 2)
title("Velocity States")
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('$\dot{x}_1$', '$\dot{x}_2$')
grid on
nexttile
plot(t_total, u_total, 'LineWidth', 2) %stairs
title("Control actions")
xlabel('Time [s]')
ylabel('Force [N]')
legend('$u_1$', '$u_2$')
grid on

set(findall(fig1h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig1h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig1h, "-property", "FontSize"), "FontSize", 18)

fig2h = figure(2);
plot(t_total, Varepsilon_total, "LineWidth",2)
title("Ergodic Metric")
xlabel('Time [s]')
ylabel('$\varepsilon \left( \mathbf{X_e}(t), \Phi(\mathbf{x}) \right) $')
grid on

set(findall(fig2h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig2h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig2h, "-property", "FontSize"), "FontSize", 18)

%% 2-D Charts

if n_iter <= 2
    aux = 1;
elseif n_iter > 2 && n_iter <= 5
    aux = 2;
elseif n_iter > 5 && n_iter <= 8
    aux = 3;
elseif n_iter > 8 && n_iter <= 11
    aux = 4;
elseif n_iter > 11 && n_iter <= 14
    aux = 5;
end

nbDrawingSeg = 100;
tmp_vec = linspace(-pi, pi, nbDrawingSeg)';
Elipse_Phi = zeros(height(tmp_vec), 2, n_def); %Elipse
for j = 1:n_def
    Elipse_Phi(:,:,j) = [cos(tmp_vec), sin(tmp_vec)] * real(Sigma_ast_Phi(:,:,j)) + repmat(Mu(j,:),nbDrawingSeg,1);
end

contornos = 25;

fig3h = figure(3);
subplot(aux,3,1)
pcolor(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "FaceColor","interp","EdgeColor","none")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Real PDF")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
axis equal tight
grid on
colormap("default")
hold on
for j = 1:n_def
    plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "-.k", "LineWidth", 1.3)
end
plot(Mu(:,1), Mu(:,2), ".", 'MarkerSize', 8)
hold off
legend('$\Phi(\mathbf{x})$','Location', 'best') %'northeastoutside')
for i = 1:n_iter
    subplot(aux,3,i+1)
    contour(x_1_grid, x_2_grid, reshape(Phi_hat_x_reg(:,:,i), length(x_2), length(x_1)) )%, contornos)
    %pcolor(x_1_grid, x_2_grid, reshape(Phi_hat_x_reg(:,:,i), length(x_2), length(x_1)),...
        % "FaceColor","interp","EdgeColor","none")
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    title("Estimated PDF, iteration " + i)
    xlabel('$x_1$ [m]')
    ylabel('$x_2$ [m]')
    axis equal
    grid on
    hold on
    plot(X_e_reg(:,1,i), X_e_reg(:,2,i),'LineWidth',1.5)
    plot(X_e_reg(1,1,i), X_e_reg(1,2,i),'ksq','MarkerSize',7,'LineWidth',1.5)
    for j = 1:n_def
        plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "-.k", "LineWidth", 1)
    end
    plot(Mu(:,1),Mu(:,2),'.k','MarkerSize',8)
    legend('$\hat{\Phi}(\mathbf{x})$','$\mathbf{X_e}(t)$', '$\mathbf{X_e}(0)$',...
        "Real" + newline + "defects",'Location','best')%'northeastoutside')
    hold off
end

set(findall(fig3h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig3h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig3h, "-property", "FontSize"), "FontSize", 16)

%% Gráficas con adición de puntos Spline

figure(4)
subplot(3,1,1)
plot(t_spline_total, X_e_spline_total, 'LineWidth', 1.5)
title("Position States",'Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Position [m]','Interpreter','latex')
legend('$x_1$', '$x_2$','Interpreter','latex')
grid on
subplot(3,1,2)
plot(t_spline_total, X_e_dot_spline_total, 'LineWidth', 1.5)
title("Velocity States",'Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Velocity [m/s]','Interpreter','latex')
legend('$\dot{x}_1$', '$\dot{x}_2$','Interpreter','latex')
grid on
subplot(3,1,3)
plot(t_spline_total, u_spline_total, 'LineWidth', 1.5) %stairs
title("Control actions",'Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Force [N]','Interpreter','latex')
legend('$u_1$', '$u_2$','Interpreter','latex')
grid on

figure(5)
subplot(3,1,1)
plot(t_spline_total, V_Xe_total, 'LineWidth', 1.5)
title("Sensor measurement on time",'Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Force [N]','Interpreter','latex')
legend('$V(t)$','Interpreter','latex')
grid on

for i = 1:n_iter

    subplot(3,1,2)
    hold on
    scatter(X_e_spline_reg(:,1,i), V_Xe_reg(:,:,i), 14, 'filled',...
        'DisplayName', "$V(x_1)$, iteration " + i, "MarkerFaceAlpha", 0.4)
    title("Sensor measurement on spatial domain",'Interpreter','latex')
    xlabel('$x_1$ [m]','Interpreter','latex')
    ylabel('Force [N]','Interpreter','latex')
    grid on
    hold off
    legend('Interpreter','latex')
    subplot(3,1,3)
    hold on
    scatter(X_e_spline_reg(:,2,i), V_Xe_reg(:,:,i), 14, 'filled',...
        'DisplayName', "$V(x_2)$, iteration " + i, "MarkerFaceAlpha", 0.4)
    title("Sensor measurement on spatial domain",'Interpreter','latex')
    xlabel('$x_2$ [m]','Interpreter','latex')
    ylabel('Force [N]','Interpreter','latex')
    grid on
    hold off
    legend('Interpreter','latex')

end


%% Charts about the GMM training

if Estim_sol(n_iter).flag_done
    limit_id = n_iter + 1;
else
    limit_id = n_iter;
end

for i = 1:limit_id
    figure(21)
    subplot(aux,3,i)
    surf(x_1_grid, x_2_grid, reshape(Phi_hat_x_reg(:,:,i), length(x_2), length(x_1)),...
        'EdgeColor','interp','FaceColor','interp')
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    title("Estimated PDF: GMM, iteration " + i,'Interpreter','latex')
    xlabel('$x_1$ [m]','Interpreter','latex')
    ylabel('$x_2$ [m]','Interpreter','latex')
    grid on
    legend('$\hat{\Phi}(\mathbf{x})$','Interpreter','latex','Location','best')
end

if n_iter <= 3
    aux_2 = 1;
elseif n_iter > 3 && n_iter <= 6
    aux_2 = 2;
elseif n_iter > 6 && n_iter <= 9
    aux_2 = 3;
elseif n_iter > 9 && n_iter <= 12
    aux_2 = 4;
elseif n_iter > 12 && n_iter <= 15
    aux_2 = 5;
end

for i = 1:n_iter
    figure(22)
    subplot(aux_2,3,i)
    hist3(Estim_sol(i).Data_Xe_hist_V,'CDataMode','auto','FaceColor','interp',...
        "EdgeColor","none",'Nbins',[length(x_1)-1, length(x_2)-1])
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    title("GMM fit, iteration " + i,'Interpreter','latex')
    legend("Training data",'Interpreter','latex', 'Location','best')
    xlabel('$x_1$ [m]','Interpreter','latex')
    ylabel('$x_2$ [m]','Interpreter','latex')
    zlabel('Measurement int','Interpreter','latex')
    grid on

    figure(23)
    subplot(aux_2,3,i)
    scatter3(X_e_spline_reg(:,1,i), X_e_spline_reg(:,2,i), V_Xe_reg(:,:,i), 10, "black")
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    title("Measurements on spatial domain, iteration " + i,'Interpreter','latex')
    xlabel('$x_1$ [m]','Interpreter','latex')
    ylabel('$x_2$ [m]','Interpreter','latex')
    zlabel('$V_k$ [N]','Interpreter','latex')
    grid on
end

%% Real PDF  vs  Estimated PDF: Resultados

GMM_hat = gmdistribution(Mu_found, Sigma_found);
Phi_hat_final = pdf(GMM_hat, Omega);
n_def_found = size(Sigma_found, 3);

nbDrawingSeg = 100;
tmp_vec = linspace(-pi, pi, nbDrawingSeg)';
stdev_Phi_hat = zeros(size(Sigma_found));
Sigma_ast_Phi_hat = zeros(size(Sigma_found));
Elipse_Phi_hat = zeros(height(tmp_vec), 2, n_def_found); %Elipse
for j = 1:n_def_found
    stdev_Phi_hat(:,:,j) = sqrtm(Sigma_found(:,:,j));
    Sigma_ast_Phi_hat(:,:,j) = 3*stdev_Phi_hat(:,:,j);
    Elipse_Phi_hat(:,:,j) = [cos(tmp_vec), sin(tmp_vec)] * real(Sigma_ast_Phi_hat(:,:,j)) +...
        repmat(Mu_found(j,:), nbDrawingSeg, 1);
end

% Si el algoritmo no pudo encontrar todos los defectos, se calculan las
% elipses de las estimaciones de los defectos no encontrados
if ~Estim_sol(end).flag_done 
    Sigma_not_found = Estim_sol(end).GMModel.Sigma;
    Mu_not_found = Estim_sol(end).GMModel.mu;
    n_def_not_found = size(Sigma_not_found, 3);

    stdev_notfound = zeros(size(Sigma_not_found));
    Sigma_ast_not_found = zeros(size(Sigma_not_found));
    Elipse_not_found = zeros(height(tmp_vec), 2, n_def_not_found); %Elipse
    for j = 1:n_def_not_found
        stdev_notfound(:,:,j) = sqrtm(Sigma_not_found(:,:,j));
        Sigma_ast_not_found(:,:,j) = 3*stdev_notfound(:,:,j);
        Elipse_not_found(:,:,j) = [cos(tmp_vec), sin(tmp_vec)] * real(Sigma_ast_not_found(:,:,j)) +...
            repmat(Mu_not_found(j,:), nbDrawingSeg, 1);
    end
end

FoundDef_color = [0.8500 0.3250 0.0980];
NotFoundDef_color = [0.4660 0.6740 0.1880];

fig24h = figure(24);
tiledlayout(2,3)

nexttile([2, 2])
contour(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    contornos, "DisplayName","Real Defects");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Real PDF vs Estimated PDF",'Interpreter','latex')
xlabel('$x_1$ [m]','Interpreter','latex')
ylabel('$x_2$ [m]','Interpreter','latex')
axis equal
grid on
hold on
foundplot = 0;
notfoundplot = 0;
%Grafica el primer defecto encontrado (si lo hay)
if n_def_found >= 1
    patch(Elipse_Phi_hat(:,1,1), Elipse_Phi_hat(:,2,1), FoundDef_color,...
            'LineWidth', 1.5, 'EdgeColor', FoundDef_color, "FaceAlpha",0.2,...
            "DisplayName", "Estimated found defects")
    foundplot = 1;
end
%Grafica los defectos no encontrados (si los hay)
if ~Estim_sol(end).flag_done 
        for i = 1:n_def_not_found
            patch(Elipse_not_found(:,1,i), Elipse_not_found(:,2,i), NotFoundDef_color,...
                'LineWidth', 1.5, 'EdgeColor', NotFoundDef_color, "FaceAlpha",0.2,...
                "DisplayName", "Estimation of Not found defects")
        end
        notfoundplot = 1;
end
%Termina de graficar los defectos encontrados (si hay más de uno)
if n_def_found > 1
    for j = 2:n_def_found
        patch(Elipse_Phi_hat(:,1,j), Elipse_Phi_hat(:,2,j), FoundDef_color,...
                'LineWidth', 1.5, 'EdgeColor', FoundDef_color, "FaceAlpha",0.2)
    end
end
%Grafica las elipses de defectos reales
for j = 1:n_def
    plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "-.k", "LineWidth", 1)
end
%Grafica los centroides
plot(Mu(:,1),Mu(:,2),'.k','MarkerSize',8)
plot(Mu_found(:,1), Mu_found(:,2), '+', 'LineWidth', 1.5, 'color', FoundDef_color);
if ~Estim_sol(end).flag_done 
    plot(Mu_not_found(:,1), Mu_not_found(:,2), '+', 'LineWidth', 1.5, 'color', NotFoundDef_color);
end
hold off
lgd = legend;
lgd.String(1 + foundplot + notfoundplot + 1:end) = [];


nexttile(3)
surf(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "EdgeColor", "none", "FaceColor", "interp")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Real PDF")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
grid on

nexttile(6)
surf(x_1_grid, x_2_grid, reshape(Phi_hat_final, length(x_2), length(x_1)),...
    "EdgeColor", "none", "FaceColor", "interp")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Estimated PDF")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
grid on

set(findall(fig24h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig24h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig24h, "-property", "FontSize"), "FontSize", 20)

%% Minimum variation constraint

iter_vec = (1:n_iter)';
variaciones = zeros(size(iter_vec));
for r = 1:n_iter
    variaciones(r) = Estim_sol(r).MinVariation;
end

fig25h = figure(25);
plot(iter_vec, variaciones, "-o", "LineWidth", 2, "MarkerSize", 10)
yline(Par_PDF.Thres_Variation, "-", "Variation threshold for finding a defect",...
    "LineWidth", 2)
title("Minimum variation constraint")
xlabel('Iteration')
ylabel('Distance [m]')
xlim([1, n_iter])
grid on

set(findall(fig25h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig25h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig25h, "-property", "FontSize"), "FontSize", 20)