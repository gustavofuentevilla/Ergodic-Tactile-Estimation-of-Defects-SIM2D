close all
clear
clc
%%
load("BO_unconsTest/Test2_new/Case1/7Def/X0_15/output_3.mat")

%% Plot results

FoundDef_color = hex2rgb("#238b45");
NotFoundDef_color = "yellow";
RealDef_color = "black";
Trayectory_color = hex2rgb("#d94801");

GMM_hat = gmdistribution(Mu_found, Sigma_found);
if GMM_hat.NumComponents > 0
    Phi_hat_final = pdf(GMM_hat, Omega);
else
    Phi_hat_final = zeros(height(Omega), 1);
end
n_def_found = size(Sigma_found, 3);

nbDrawingSeg = 100;
tmp_vec = linspace(-pi, pi, nbDrawingSeg)';
Elipse_Phi = zeros(height(tmp_vec), 2, n_def); %Elipse
for j = 1:n_def
    Elipse_Phi(:,:,j) = [cos(tmp_vec), sin(tmp_vec)] * ...
                        real(3*sqrtm(Sigma(:,:,j))) + ...
                        repmat(Mu(j,:),nbDrawingSeg,1);
end

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

fig10h = figure(10);
layout10h = tiledlayout(fig10h, 1, 1);

nexttile
title("Result",'Interpreter','latex')
xticks([L_1_l, L_1_u])
yticks([L_2_l, L_2_u])
xtickformat('%.1f')
ytickformat('%.1f')
axis square
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
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
hold off

%Grafica los defectos encontrados (si los hay)
if n_def_found >= 1
    for j = 1:n_def_found
        F_def_ax(j) = patch(Elipse_Phi_hat(:,1,j), Elipse_Phi_hat(:,2,j), ...
            FoundDef_color,...
            'LineWidth', 3, 'EdgeColor', FoundDef_color, "FaceAlpha",0.2);
    end
end

xlabel(layout10h, '$x_1$ [m]','Interpreter','latex', "FontSize", 22)
ylabel(layout10h, '$x_2$ [m]','Interpreter','latex', "FontSize", 22)

set(findall(fig10h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig10h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig10h, "-property", "FontSize"), "FontSize", 22)

%% Performance plots

t_plot = linspace(0, t_exec(end), height(X_e_sp))';

fig11h = figure(11);
tiledlayout(3,1)

nexttile(1)
plot(t_plot, X_e_sp, "LineWidth", 3)
xlabel('Time [s]')
ylabel('Position [m]')
title("Position")
legend("$x_{e_1}$", "$x_{e_2}$")

v_normplot = sqrt(sum((diff(X_e_sp)/T_s_interp).^2, 2));
v_normplot(end + 1) = v_normplot(end);

nexttile(2)
plot(t_plot, v_normplot, "LineWidth", 3)
xlabel('Time [s]')
ylabel('Velocity [m/s]')
title('Velocity norm')
grid on;

samples = BO_samples;
nexttile(3)
tmp = (1:samples)';
loglikelihood_plot = plot(tmp, loglikelihood(end-samples+1:end), "LineWidth", 2);
title("Log-likelihood")
xlabel('Iteration')
grid on

set(findall(fig11h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig11h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig11h, "-property", "FontSize"), "FontSize", 20)

%% GPR scope

fig12h = figure(12);
tiledlayout(2,2)

nexttile(1)
TrueF_plot = surf(x_1_grid, x_2_grid, reshape(V_real, length(x_2), length(x_1)),...
                 "EdgeColor", "none", "FaceColor", "interp");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Ground Truth Function to be Optimized")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
grid on

nexttile(2)
Pred_plot = surf(x_1_grid, x_2_grid, reshape(V_pred, length(x_2), length(x_1)),...
                "EdgeColor", "none", "FaceColor", "interp");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("GPR Prediction")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
grid on

nexttile(3)
sd_plot = surf(x_1_grid, x_2_grid, reshape(sd, length(x_2), length(x_1)),...
        "EdgeColor", "none", "FaceColor", "interp");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title('Uncertainty')

nexttile(4)
EI_plot = surf(x_1_grid, x_2_grid, reshape(EI, length(x_2), length(x_1)),...
                "EdgeColor", "none", "FaceColor", "interp");
hold on
xPosEI_plot = plot3(xEI([1, 1]),xrange(2,:),eimax*[1 1],'--k','LineWidth',1.5);
yPosEI_plot = plot3(xrange(1,:),xEI([2, 2]),eimax*[1 1],'--k','LineWidth',1.5); 
hold off; 
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title('Expected Improvement'); 
grid on;

set(findall(fig12h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig12h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig12h, "-property", "FontSize"), "FontSize", 22)

%% Plane trajectories and results

fig13h = figure(13);
tiledlayout(1,4)

nexttile(1)
pcolor(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "EdgeColor", "none", "FaceColor", "interp")
c = colorbar;
c.Ticks = [c.Limits(1), c.Limits(2)];
c.TickLabels = {'0', '1'};
hold on
Traj_plot = plot(X_e(:,1), X_e(:,2),'Color', Trayectory_color,...
                'LineWidth', 2, 'Marker','*', 'MarkerSize',8);
plot(X_e(1,1), X_e(1,2), 'ksq', 'LineWidth',8, 'MarkerSize',10)
plot(X_e(end,1), X_e(end,2), 'ko', 'LineWidth',8, 'MarkerSize',10)
hold off
axis equal
xticks([L_1_l, L_1_u])
yticks([L_2_l, L_2_u])
xtickformat('%.1f')
ytickformat('%.1f')
axis square
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
xlabel("$x_1$ [m]")
ylabel("$x_2$ [m]")

% nexttile(2)
% Est_plot = pcolor(x_1_grid, x_2_grid, reshape(V_pred, length(x_2), length(x_1)),...
%                  "EdgeColor", "none", "FaceColor", "interp");
% axis equal
% xticks([L_1_l, L_1_u])
% yticks([L_2_l, L_2_u])
% xtickformat('%.1f')
% ytickformat('%.1f')
% axis square
% xlim([L_1_l, L_1_u])
% ylim([L_2_l, L_2_u])
% 
% nexttile(3)
% xticks([L_1_l, L_1_u])
% yticks([L_2_l, L_2_u])
% xtickformat('%.1f')
% ytickformat('%.1f')
% axis square
% xlim([L_1_l, L_1_u])
% ylim([L_2_l, L_2_u])
% hold on

% For equal level of severities on the defects
%
% for j = 1:n_def
%     realdef_ax(j) = plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j),...
%                         "-.", "LineWidth", 3,...
%                         "Color", RealDef_color);
% end
% plot(Mu(:,1),Mu(:,2),'.','MarkerSize',15, "Color", RealDef_color)

% For different level of severities on the defects
%
% for j = 1:n_def
%     realdef_ax(j) = plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j),...
%                         "-.", "LineWidth", 25*proportions(j),...
%                         "Color", RealDef_color);
%     plot(Mu(j,1),Mu(j,2),'.','MarkerSize',90*proportions(j),...
%         "Color", RealDef_color)
% end    


%Grafica las elipses de defectos reales

% For equal level of severities
%
% for j = 1:n_def
%     R_def_ax(j) = plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "-.",...
%                         "LineWidth", 3, "Color", RealDef_color);
% end
% plot(Mu(:,1),Mu(:,2),'.','MarkerSize',15,"Color",RealDef_color)

% For different level of severities
%
% for j = 1:n_def
%     R_def_ax(j) = plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j),...
%                         "-.", "LineWidth", 25*proportions(j),...
%                         "Color", RealDef_color);
%     plot(Mu(j,1),Mu(j,2),'.','MarkerSize',90*proportions(j),...
%         "Color", RealDef_color)
% end 

% plot(X_e(:,1), X_e(:,2),'Color', Trayectory_color,...
%                 'LineWidth', 2, 'Marker','*', 'MarkerSize',8);
% 
% plot(X_e(1,1), X_e(1,2), 'ksq', 'LineWidth',4, 'MarkerSize',12)
% plot(X_e(end,1), X_e(end,2), 'ko', 'LineWidth',4, 'MarkerSize',12)
% hold off



nexttile(2)
title("Result",'Interpreter','latex')
xticks([L_1_l, L_1_u])
yticks([L_2_l, L_2_u])
xtickformat('%.1f')
ytickformat('%.1f')
axis square
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
xlabel("$x_1$ [m]")
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
hold off

%Grafica los defectos encontrados (si los hay)
if n_def_found >= 1
    for j = 1:n_def_found
        F_def_ax(j) = patch(Elipse_Phi_hat(:,1,j), Elipse_Phi_hat(:,2,j), ...
            FoundDef_color,...
            'LineWidth', 3, 'EdgeColor', FoundDef_color, "FaceAlpha",0.2);
    end
end

xlabel(layout10h, '$x_1$ [m]','Interpreter','latex', "FontSize", 22)
ylabel(layout10h, '$x_2$ [m]','Interpreter','latex', "FontSize", 22)

set(findall(fig13h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig13h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig13h, "-property", "FontSize"), "FontSize", 22)

colormap(brewermap(15,"-Blues"))

