%% Charts

close all
clearvars -except Erg_traj_ipopt
clc

load("Test3/Case2/7Def/X0_15/output_3.mat")

%% Reconstrucción de distribución empírica y métrica ergódica

Varepsilon_reg = zeros(length(t), 1, n_iter);
C_x_reg = zeros(height(Omega), length(t), n_iter);

for r = 1:n_iter
    c_k = zeros(K^n, 1);
    C_x = zeros(height(Omega), 1);
    for i = 1:length(t)

        % Compute Fourier Functions and coefficients on the new position
        f_k_traj = prod(cos( K_cal'.*pi.*(X_e_reg(i,:,r) - L_i_l)./(L_i_u - L_i_l) ), 2) ./ h_k_reg ;
        c_k = c_k + (f_k_traj*T_s)/(t_f) ; %i*T_s %t_f

        %Ergodic metric
        Varepsilon = sum( Lambda_k .* (c_k - phi_k_REG(:,:,r)).^2 );

        %Empirical distribution reconstruction
        C_x_i = zeros(height(Omega), 1);
        for j = 1:K^n
            C_x_i = C_x_i + c_k(j)*f_k_reg(:,j);
        end

        %Se suman todas las distribuciones generadas en cada muestra
        C_x = C_x + C_x_i;

        %Se registra
        C_x_reg(:,i,r) = C_x;
        Varepsilon_reg(i,:,r) = Varepsilon;

    end

end


%% Pre-Procesing for charts

t_total = zeros(n_iter*(length(t)-1) + 1, 1);
X_e_total = zeros(n_iter*(length(t)-1) + 1, 2);
X_e_dot_total = zeros(n_iter*(length(t)-1) + 1, 2);
Varepsilon_total = zeros(n_iter*(length(t)-1) + 1, 1);
u_total = zeros(n_iter*(length(t)-1) + 1, 2);

t_spline_total = zeros(n_iter*(length(t_spline)-1) + 1, 1);
X_e_spline_total = zeros(n_iter*(length(t_spline)-1) + 1, 2);
X_e_dot_spline_total = zeros(n_iter*(length(t_spline)-1) + 1, 2);
u_spline_total = zeros(n_iter*(length(t_spline)-1) + 1, 2);
V_Xe_total = zeros(n_iter*(length(t_spline)-1) + 1, 1);

for i = 1:n_iter

    id_init = ((i - 1)*(length(t)-1) + 1); %1,101,..
    id_last = (i*(length(t)-1) + 1);        %101, 201,...

    t_total( id_init:id_last ) = (i - 1)*t(end) + t;
    X_e_total( id_init:id_last, : ) = X_e_reg(:,:,i);
    X_e_dot_total( id_init:id_last, : ) = X_e_dot_reg(:,:,i);
    Varepsilon_total( id_init:id_last, : ) = Varepsilon_reg(:,:,i);
    u_total( id_init:id_last, : ) = [u_reg(:,:,i); [NaN, NaN]];
    
    id_init_spline = ((i - 1)*(length(t_spline)-1) + 1);
    id_last_spline = (i*(length(t_spline)-1) + 1);
    t_spline_total( id_init_spline:id_last_spline ) = (i - 1)*t_spline(end) + t_spline;
    X_e_spline_total( id_init_spline:id_last_spline, : ) = X_e_spline_reg(:,:,i);
    X_e_dot_spline_total( id_init_spline:id_last_spline, : ) = X_e_dot_spline_reg(:,:,i);
    u_spline_total( id_init_spline:id_last_spline, : ) = [u_spline_reg(:,:,i); [NaN, NaN]];
    V_Xe_total( id_init_spline:id_last_spline, : ) = V_Xe_reg(:,:,i);

end

%% Definitions
% For Fig3
nbDrawingSeg = 100;
tmp_vec = linspace(-pi, pi, nbDrawingSeg)';
Elipse_Phi = zeros(height(tmp_vec), 2, n_def); %Elipse
for j = 1:n_def
    Elipse_Phi(:,:,j) = [cos(tmp_vec), sin(tmp_vec)] * ...
                        real(3*sqrtm(Sigma(:,:,j))) + ...
                        repmat(Mu(j,:),nbDrawingSeg,1);
end

% For fig24
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
    Elipse_Phi_hat(:,:,j) = [cos(tmp_vec), sin(tmp_vec)]* ...
                            real(Sigma_ast_Phi_hat(:,:,j)) +...
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
        Elipse_not_found(:,:,j) = [cos(tmp_vec), sin(tmp_vec)]* ...
                                real(Sigma_ast_not_found(:,:,j)) +...
                                repmat(Mu_not_found(j,:), nbDrawingSeg, 1);
    end
end

FoundDef_color = [0.8500 0.3250 0.0980];
NotFoundDef_color = [0.4660 0.6740 0.1880];
contornos = 25;

% For fig25
iter_vec = (1:n_iter)';
variaciones = zeros(size(iter_vec));
for r = 1:n_iter
    variaciones(r) = Estim_sol(r).MinVariation;
end

D_KL_reg = zeros(1, n_iter);
for i = 1:n_iter
    D_KL_reg(i) = Estim_sol(i).D_KL; 
end

%For fig26
D_KL_bar_u = Par_PDF.D_KL_bar_u;
D_KL_tmp = (0:0.01:D_KL_bar_u+2)';
Variation_thres_eps = Par_PDF.Thres_Variation - 0.001;

MaxVarCons = Par_PDF.MaxVarCons;
% nu_p = Par_PDF.nu_p;

func = D_KL_tmp*nu_p*MaxVarCons/D_KL_bar_u;

func(func <= Variation_thres_eps) = Variation_thres_eps;
func(func >= nu_p*MaxVarCons) = nu_p*MaxVarCons;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHARTS xd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
set(findall(fig1h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig1h, "-property", "FontSize"), "FontSize", 18)

fig2h = figure(2);
plot(t_total, Varepsilon_total, "LineWidth",2)
title("Ergodic Metric")
xlabel('Time [s]')
ylabel('$\varepsilon \left( \mathbf{X_e}(t), \Phi(\mathbf{x}) \right) $')
grid on

set(findall(fig2h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig2h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig2h, "-property", "FontSize"), "FontSize", 18)

%% 2-D Charts

columnas = 3;
filas = ceil((n_iter + 1)/columnas);

fig3h = figure(3);
subplot(filas, columnas, 1)
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
    plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "-.r", "LineWidth", 1.3)
end
plot(Mu(:,1), Mu(:,2), ".", 'MarkerSize', 8)
hold off
legend('$\Phi(\mathbf{x})$','Location', 'northeastoutside')
for i = 1:n_iter
    subplot(filas, columnas, i+1)
    pcolor(x_1_grid, x_2_grid, ...
            reshape(Phi_hat_x_reg(:,:,i), length(x_2), length(x_1)),...
            "EdgeColor","none", "FaceColor","interp")
    %pcolor(x_1_grid, x_2_grid, ...
    %       reshape(Phi_hat_x_reg(:,:,i), length(x_2), length(x_1)),...
        %  "FaceColor","interp","EdgeColor","none")
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    title("Estimated PDF, iteration " + i)
    xlabel('$x_1$ [m]')
    ylabel('$x_2$ [m]')
    axis equal tight
    grid on
    hold on
    plot(X_e_reg(:,1,i), X_e_reg(:,2,i),"black",'LineWidth',1.5)
    plot(X_e_reg(1,1,i), X_e_reg(1,2,i),'ksq', ...
        'MarkerSize',7,'LineWidth',1.5)
    for j = 1:n_def
        plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "-.r", "LineWidth", 1)
    end
    plot(Mu(:,1),Mu(:,2),'.r','MarkerSize',8)
    legend('$\hat{\Phi}(\mathbf{x})$','$\mathbf{X_e}(t)$', ...
        '$\mathbf{X_e}(0)$',...
        "Real" + newline + "defects",'Location','northeastoutside')
    hold off
end

set(findall(fig3h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig3h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
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

columnas = 3;
if Estim_sol(n_iter).flag_done
    limit_id = n_iter + 1;
    filas = ceil((n_iter + 1)/columnas);
else
    limit_id = n_iter;
    filas = ceil((n_iter)/columnas);
end

for i = 1:limit_id
    fig21h = figure(21);
    subplot(filas, columnas, i)
    surf(x_1_grid, x_2_grid, ...
        reshape(Phi_hat_x_reg(:,:,i), length(x_2), length(x_1)),...
        'EdgeColor','interp','FaceColor','interp')
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    title("Estimated PDF: GMM, iteration " + i,'Interpreter','latex')
    xlabel('$x_1$ [m]','Interpreter','latex')
    ylabel('$x_2$ [m]','Interpreter','latex')
    grid on
    legend('$\hat{\Phi}(\mathbf{x})$',...
            'Interpreter','latex','Location','best')
end

columnas = 3;
filas = ceil((n_iter)/columnas);

for i = 1:n_iter
    fig22h = figure(22);
    subplot(filas, columnas, i)
    hist3(Estim_sol(i).Data_Xe_hist_V,'CDataMode','auto', ...
        'FaceColor','interp',...
        "EdgeColor","none",'Nbins',[length(x_1)-1, length(x_2)-1])
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    title("Histogram for GMM fitting, iteration " + i,'Interpreter','latex')
    legend("Training data",'Interpreter','latex', 'Location','best')
    xlabel('$x_1$ [m]','Interpreter','latex')
    ylabel('$x_2$ [m]','Interpreter','latex')
    % zlabel('Measurement int','Interpreter','latex')
    grid on

    fig23h = figure(23);
    subplot(filas, columnas, i)
    scatter3(X_e_spline_reg(:,1,i), X_e_spline_reg(:,2,i), ...
            V_Xe_reg(:,:,i), 10, "black")
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    title("Measurements on spatial domain, iteration " + i, ...
          'Interpreter','latex')
    xlabel('$x_1$ [m]','Interpreter','latex')
    ylabel('$x_2$ [m]','Interpreter','latex')
    zlabel('$V_k$ [N]','Interpreter','latex')
    grid on
end

set(findall(fig21h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig21h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig21h, "-property", "FontSize"), "FontSize", 20)

set(findall(fig22h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig22h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig22h, "-property", "FontSize"), "FontSize", 20)

set(findall(fig23h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig23h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig23h, "-property", "FontSize"), "FontSize", 20)

%% Real PDF  vs  Estimated PDF: Resultados

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
            'LineWidth', 1.5, 'EdgeColor', FoundDef_color, ...
            "FaceAlpha",0.2,"DisplayName", "Estimated found defects")
    foundplot = 1;
end
%Grafica los defectos no encontrados (si los hay)
if ~Estim_sol(end).flag_done 
        for i = 1:n_def_not_found
            patch(Elipse_not_found(:,1,i), Elipse_not_found(:,2,i), ...
                NotFoundDef_color,'LineWidth', 1.5, 'EdgeColor', ...
                NotFoundDef_color, "FaceAlpha",0.2,...
                "DisplayName", "Estimation of Not found defects")
        end
        notfoundplot = 1;
end
%Termina de graficar los defectos encontrados (si hay más de uno)
if n_def_found > 1
    for j = 2:n_def_found
        patch(Elipse_Phi_hat(:,1,j), Elipse_Phi_hat(:,2,j), ...
            FoundDef_color,...
            'LineWidth', 1.5, 'EdgeColor', FoundDef_color, "FaceAlpha",0.2)
    end
end
%Grafica las elipses de defectos reales
for j = 1:n_def
    plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "-.k", "LineWidth", 1)
end
%Grafica los centroides
plot(Mu(:,1),Mu(:,2),'.k','MarkerSize',8)
plot(Mu_found(:,1), Mu_found(:,2), '+', ...
    'LineWidth', 1.5, 'color', FoundDef_color);
if ~Estim_sol(end).flag_done 
    plot(Mu_not_found(:,1), Mu_not_found(:,2), '+', ...
        'LineWidth', 1.5, 'color', NotFoundDef_color);
end
hold off
lgd = legend;
lgd.String(1 + foundplot + notfoundplot + 1:end) = [];
lgd.Location = "best";


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
surf(x_1_grid, x_2_grid, ...
    reshape(Phi_hat_final, length(x_2), length(x_1)),...
    "EdgeColor", "none", "FaceColor", "interp")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Estimated PDF")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
grid on

set(findall(fig24h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig24h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig24h, "-property", "FontSize"), "FontSize", 20)

%% Minimum variation constraint

fig25h = figure(25);

yyaxis left
VarCons_plot = plot(iter_vec, variaciones, "b-o", ...
                    "LineWidth", 2, "MarkerSize", 10);
yline(Par_PDF.Thres_Variation, "b--", ...
    "Threshold for finding a defect = " + Par_PDF.Thres_Variation,...
    "LineWidth", 2, "LabelHorizontalAlignment","left")
xlabel('Iteration')
ylabel('Distance [m]')
xlim([1, n_iter])
grid on
ax = gca;
ax.YColor = 'b';

yyaxis right
D_KL_plot = plot(iter_vec, D_KL_reg, "-o", ...
                "LineWidth", 2, "MarkerSize", 10);
ylabel("Information [nat]")
grid on

legend([VarCons_plot, D_KL_plot],...
       {"Minimum Variation Constraint", ...
       "$D_{KL}$ measure from Prior to Posterior PDF"})

set(findall(fig25h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig25h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig25h, "-property", "FontSize"), "FontSize", 20)

%% Function that maps D_KL to Variation Constraint

fig26h = figure(26);
plot(D_KL_tmp, func, "LineWidth", 2)
yline(Par_PDF.Thres_Variation, "-", ...
    "Variation threshold for finding a defect",...
    "LineWidth", 2)
xline(D_KL_bar_u, "-r", "$\bar{D}_{KL}$ = " + D_KL_bar_u, ...
     "LineWidth", 2, "LabelVerticalAlignment","middle")
title("Minimum variation constraint as a function of $KL$ divergence")
xlabel("Information [nat]")
ylabel("Distance [m]")
xtickformat('%.1f')
ytickformat('%.2f')
grid on
% hold on
% plot(D_KL_bar_u,nu_p*MaxVarCons,"sq","MarkerSize",14,"LineWidth",2.5)
% hold off

set(findall(fig26h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig26h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig26h, "-property", "FontSize"), "FontSize", 20)

%% Publication Figs

% Calcula el número de filas dependiendo del número de iteraciones
% suponiendo un número de columnas totales
FoundDef_color = hex2rgb("#238b45"); %hex2rgb("#d94801");
NotFoundDef_color = "yellow";
RealDef_color = "black";
Trayectory_color = hex2rgb("#d94801");  %hex2rgb("#045a8d"); %"black";

columnas = 6;
filas = ceil((n_iter + 1)/columnas);

fig30h = figure(30);
layout30h = tiledlayout(fig30h, filas, columnas);

for i = 1:n_iter
    SigF_tmp(i).Sigma_found = Estim_sol(i).Sigma_found;
    MuF_tmp(i).Mu_found = Estim_sol(i).Mu_found;
end

for i = 1:n_iter

    nexttile(layout30h)

    probmap_ax = pcolor(x_1_grid, x_2_grid, ...
            reshape(Phi_hat_x_reg(:,:,i), length(x_2), length(x_1)),...
            "EdgeColor","none", "FaceColor","interp");
    title("Iteration " + i)
    % xlabel('$x_1$ [m]')
    % ylabel('$x_2$ [m]')
    xtickformat('%.1f')
    ytickformat('%.1f')
    axis equal tight
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    grid on

    if i == n_iter
        c = colorbar;
        c.Ticks = [c.Limits(1), c.Limits(2)];
        c.TickLabels = {'0', '1'};
        % c.Label.String = "Probability";
    end

    hold on
    
    % For equal level of severities on the defects
    %
    for j = 1:n_def
        realdef_ax(j) = plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j),...
                            "-.", "LineWidth", 3,...
                            "Color", RealDef_color);
    end
    plot(Mu(:,1),Mu(:,2),'.','MarkerSize',15, "Color", RealDef_color)

    % For different level of severities on the defects
    %
    % for j = 1:n_def
    %     realdef_ax(j) = plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j),...
    %                         "-.", "LineWidth", 25*proportions(j),...
    %                         "Color", RealDef_color);
    %     plot(Mu(j,1),Mu(j,2),'.','MarkerSize',90*proportions(j),...
    %         "Color", RealDef_color)
    % end    

    % Plotting \pi_d text
    %
    % if i == 1
    %     for j = 1:height(Mu)
    %         text(Mu(j,1) + 0.04, Mu(j,2), "$\pi_" + j + "$",...
    %              'FontSize', 8)
    %     end
    % end

    % ----Código para graficar los defectos ya encontrados en la i-esima it
    if i > 1
        nbDrawingSeg = 100;
        tmp_vec = linspace(-pi, pi, nbDrawingSeg)';
        Sigma_tmp = cat(3, SigF_tmp(1:i-1).Sigma_found);
        Mu_tmp = cat(1, MuF_tmp(1:i-1).Mu_found);
        if ~isempty(Sigma_tmp)
            stdev_tmp = zeros(size(Sigma_tmp));
            Sigma_ast_tmp = zeros(size(Sigma_tmp));
            Elipse_tmp = zeros(height(tmp_vec), 2, size(Sigma_tmp, 3));
            for j = 1:size(Sigma_tmp, 3)
                stdev_tmp(:,:,j) = sqrtm(Sigma_tmp(:,:,j));
                Sigma_ast_tmp(:,:,j) = 3*stdev_tmp(:,:,j);
                Elipse_tmp(:,:,j) = [cos(tmp_vec), sin(tmp_vec)]* ...
                                        real(Sigma_ast_tmp(:,:,j)) +...
                                        repmat(Mu_tmp(j,:), nbDrawingSeg, 1);
            end
            for j = 1:size(Sigma_tmp, 3)
                realdef_ax(j) = plot(Elipse_tmp(:,1,j), Elipse_tmp(:,2,j),...
                                    "-.", "LineWidth", 3,...
                                    "Color", FoundDef_color);
            end
            plot(Mu_tmp(:,1),Mu_tmp(:,2),'.',...
                'MarkerSize',15, "Color", FoundDef_color)
            for j = 1:size(Sigma_tmp, 3)
                F_def_ax(j) = patch(Elipse_tmp(:,1,j), Elipse_tmp(:,2,j),...
                                    FoundDef_color, 'LineWidth', 3,...
                                    'EdgeColor', FoundDef_color,...
                                    "FaceAlpha",0.2);
            end
        end
    end

    traj_ax = plot(X_e_reg(:,1,i), X_e_reg(:,2,i),...
                    "Color", Trayectory_color,'LineWidth',3);
    traj0_ax = plot(X_e_reg(1,1,i), X_e_reg(1,2,i),'sq', "Color",...
                    Trayectory_color, 'MarkerSize',7,'LineWidth',10);
    
    % lgd = legend([probmap_ax, realdef_ax(1), traj_ax, traj0_ax],...
    %         {"$\hat{\Phi}(\mathbf{x})$",...
    %          "Real" + newline + "defects",...
    %          "$\mathbf{X_e}(t)$",...
    %          "$\mathbf{X_e}(0)$"});
    % lgd.Location = "northeastoutside";

    hold off

    % Ticks
    xticks([L_1_l, L_1_u])
    yticks([L_2_l, L_2_u])
end

nexttile(layout30h)

% pcolor(x_1_grid, x_2_grid, reshape(Phi_hat_x_reg(:,:,n_iter+1),...
%        length(x_2), length(x_1)),...
%        "EdgeColor","none", "FaceColor","interp")
title("Result",'Interpreter','latex')
xticks([L_1_l, L_1_u])
yticks([L_2_l, L_2_u])
xtickformat('%.1f')
ytickformat('%.1f')
axis square
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
hold on

%Grafica las elipses de defectos reales

% For equal level of severities
%
for j = 1:n_def
    R_def_ax(j) = plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "-.",...
                        "LineWidth", 3, "Color", RealDef_color);
end
plot(Mu(:,1),Mu(:,2),'.','MarkerSize',15,"Color",RealDef_color)

% For different level of severities
%
% for j = 1:n_def
%     R_def_ax(j) = plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j),...
%                         "-.", "LineWidth", 25*proportions(j),...
%                         "Color", RealDef_color);
%     plot(Mu(j,1),Mu(j,2),'.','MarkerSize',90*proportions(j),...
%         "Color", RealDef_color)
% end 

% Plot estimated centers
%
plot(Mu_found(:,1), Mu_found(:,2), '+', ...
    'LineWidth', 3, 'color', FoundDef_color);
if ~Estim_sol(end).flag_done 
    plot(Mu_not_found(:,1), Mu_not_found(:,2), '+', ...
        'LineWidth', 3, 'color', NotFoundDef_color);
end
hold off

%Grafica los defectos encontrados (si lo hay)
if n_def_found >= 1
    for j = 1:n_def_found
        F_def_ax(j) = patch(Elipse_Phi_hat(:,1,j), Elipse_Phi_hat(:,2,j), ...
            FoundDef_color,...
            'LineWidth', 3, 'EdgeColor', FoundDef_color, "FaceAlpha",0.2);
    end
end

%Grafica los defectos no encontrados (si los hay)
if ~Estim_sol(end).flag_done 
        for i = 1:n_def_not_found
            NF_def_ax(i) = patch(Elipse_not_found(:,1,i), Elipse_not_found(:,2,i), ...
                NotFoundDef_color,'LineWidth', 3, 'EdgeColor', ...
                NotFoundDef_color, "FaceAlpha",0.2);
        end
        notfoundplot = 1;
end

% Asignar leyendas
% if ~Estim_sol(end).flag_done
%     lgd = legend([R_def_ax(1)  F_def_ax(1) NF_def_ax(1)],...
%             {'Real Defects','Found Defects', 'Not Found Defects'});
% else
%     lgd = legend([R_def_ax(1)  F_def_ax(1)],...
%             {"Real" + newline + "defects","Found" + newline + "Defects"});
% end
% 
% lgd.Location = 'northeastoutside';

layout30h.TileSpacing = 'tight';
layout30h.Padding = 'tight';

% Remueve los números del eje Y en las últimas 5 gráficas

layout30h.Children(1).YTick = [];
% layout30h.Children(2).YTick = [];
layout30h.Children(3).YTick = [];
layout30h.Children(4).YTick = [];
layout30h.Children(5).YTick = [];
layout30h.Children(6).YTick = [];

xlabel(layout30h, '$x_1$ [m]','Interpreter','latex', "FontSize", 22)
ylabel(layout30h, '$x_2$ [m]','Interpreter','latex', "FontSize", 22)

set(findall(fig30h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig30h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig30h, "-property", "FontSize"), "FontSize", 22)
%%
colormap(brewermap(15,"-Blues"))