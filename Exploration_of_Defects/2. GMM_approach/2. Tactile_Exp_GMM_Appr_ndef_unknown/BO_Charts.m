%% Plot results

RealDef_color = "black";
FoundDef_color = [0.8500 0.3250 0.0980];

GMM_hat = gmdistribution(Mu_found, Sigma_found);
Phi_hat_final = pdf(GMM_hat, Omega);
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

fig3h = figure(3);
layout3h = tiledlayout(fig3h, 1, 1);

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
% for j = 1:n_def
%     realdef_ax(j) = plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j),...
%                         "-.", "LineWidth", 3,...
%                         "Color", RealDef_color);
% end
% plot(Mu(:,1),Mu(:,2),'.','MarkerSize',15, "Color", RealDef_color)

% For different level of severities on the defects
%
for j = 1:n_def
    realdef_ax(j) = plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j),...
                        "-.", "LineWidth", 25*proportions(j),...
                        "Color", RealDef_color);
    plot(Mu(j,1),Mu(j,2),'.','MarkerSize',90*proportions(j),...
        "Color", RealDef_color)
end    


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
for j = 1:n_def
    R_def_ax(j) = plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j),...
                        "-.", "LineWidth", 25*proportions(j),...
                        "Color", RealDef_color);
    plot(Mu(j,1),Mu(j,2),'.','MarkerSize',90*proportions(j),...
        "Color", RealDef_color)
end 

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

xlabel(layout3h, '$x_1$ [m]','Interpreter','latex', "FontSize", 22)
ylabel(layout3h, '$x_2$ [m]','Interpreter','latex', "FontSize", 22)

set(findall(fig3h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig3h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig3h, "-property", "FontSize"), "FontSize", 22)

