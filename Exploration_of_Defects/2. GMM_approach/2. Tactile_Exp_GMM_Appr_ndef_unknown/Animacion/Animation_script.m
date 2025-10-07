close all
clear
clc

%% Load data
% 10 -- good results
% 86 -- Not so good
load("Results/output_30.mat");

%% Preliminary computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Real defect ellipses
nbDrawingSeg = 100;
tmp_vec = linspace(-pi, pi, nbDrawingSeg)';
Elipse_Phi = zeros(height(tmp_vec), 2, n_def);
for j = 1:n_def
    Elipse_Phi(:,:,j) = [cos(tmp_vec), sin(tmp_vec)] * real(Sigma_ast_Phi(:,:,j)) + repmat(Mu(j,:),nbDrawingSeg,1);
end

% Si el algoritmo no pudo encontrar todos los defectos, se calculan las
% elipses de las estimaciones de los defectos no encontrados
if ~Estim_sol(end).flag_done 
    Sigma_not_found = Estim_sol(end).GMModel.Sigma;
    Mu_not_found = Estim_sol(end).GMModel.mu;
    n_def_not_found = size(Sigma_not_found, 3);

    SD_notfound = zeros(size(Sigma_not_found));
    Sigma_ast_DefNotFound = zeros(size(Sigma_not_found));
    Ellipse_DefNotFound = zeros(height(tmp_vec), 2, n_def_not_found); %Elipse
    for r = 1:n_def_not_found
        SD_notfound(:,:,r) = sqrtm(Sigma_not_found(:,:,r));
        Sigma_ast_DefNotFound(:,:,r) = 3*SD_notfound(:,:,r);
        Ellipse_DefNotFound(:,:,r) = [cos(tmp_vec), sin(tmp_vec)] * real(Sigma_ast_DefNotFound(:,:,r)) +...
            repmat(Mu_not_found(r,:), nbDrawingSeg, 1);
    end
end

%% Phi_hat result (for last frame)
GMM_hat = gmdistribution(Mu_found, Sigma_found);
Phi_hat_final = pdf(GMM_hat, Omega);

FoundDef_color = [0.8500 0.3250 0.0980];
NotFoundDef_color = [0.4660 0.6740 0.1880];

%% First Frame

figh = figure;
figh.Units = "centimeters";
figh.Position = [5, 5, 50, 30];
figh.Units = "pixels";

layouth = tiledlayout(figh, 3, 4);
colormap default

% %%%%%%%%%%%%%%% Gráfica en el plano %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plane_tile = nexttile(layouth, [3, 3]);
Est_PDF_plot = pcolor(x_1_grid, x_2_grid,...
                reshape(Phi_hat_x_reg(:,:,1), length(x_2), length(x_1)),...
                "FaceColor","interp","EdgeColor","none");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
%title("Real PDF vs Estimated PDF")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
axis equal tight
grid on
hold on
%Grafica las elipses de defectos reales
for j = 1:n_def
    RealDef_plot = plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "-.k", "LineWidth", 1);
end
%Grafica los centroides
plot(Mu(:,1),Mu(:,2),'.k','MarkerSize',8)
trayectoria_plot = plot(X_e_reg(1,1,1), X_e_reg(1,2,1),...
                        "Color", "black", 'LineWidth',2);
X_e_0_plot = plot(X_e_reg(1,1,1), X_e_reg(1,2,1),'ksq',...
                'MarkerSize',15,'LineWidth',2);
X_e_act_plot = plot(X_e_reg(1,1,1), X_e_reg(1,2,1),"ok",...
                    'MarkerSize',15,'LineWidth',2);
hold off
legend([Est_PDF_plot, RealDef_plot, trayectoria_plot, X_e_0_plot],...
       {"PDF Estimation, $\hat{\Phi}(\mathbf{x})$",...
       "Real defects",...
       "$\mathbf{X_e}(t)$",...
       "$\mathbf{X_e}(0)$"},...
       'Location','northeast')

% %%%%%%%%%%%%%%%%%%%% Gráficas de superficies %%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile(layouth, 4)
surf(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "EdgeColor","none", "FaceColor","interp")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$\Phi(\mathbf{x})$')
title("Real PDF")
grid on
axis tight

nexttile(layouth, 8)
Est_PDF_surf_plot = surf(x_1_grid, x_2_grid,...
    reshape(Phi_hat_x_reg(:,:,1), length(x_2), length(x_1)),...
    "FaceColor","interp","EdgeColor","none");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$\hat{\Phi}(\mathbf{x})$')
title("PDF Estimation")
grid on
axis tight

nexttile(layouth, 12)
Emp_dist_surf_plot = surf(x_1_grid, x_2_grid,...
    reshape(C_x_reg(:,1,1), length(x_2), length(x_1)),...
    "EdgeColor","none", "FaceColor","interp");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$C(\mathbf{x})$')
title("Trajectory's empirical distribution")
grid on
axis tight

% %%%%%%%%%%%%%% Setting general plot properties %%%%%%%%%%%%%%%%%%%%%%%%%%
set(findall(figh,'-property','Interpreter'),'Interpreter','latex') 
set(findall(figh,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(figh, "-property", "FontSize"), "FontSize", 20)

% %%%%%%%%%%%%%% Frames Pre-allocation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loops = length(t)*n_iter;
MovieVector(loops) = struct("cdata", [], "colormap", []);

%% Loop

for j = 1:n_iter

    title(layouth, "Defects Localization Using GMM, iteration " + j,...
    "interpreter", "latex", "FontSize", 30)
    % Actualizar Phi_hat
    Est_PDF_plot.CData = reshape(Phi_hat_x_reg(:,:,j), length(x_2), length(x_1));
    Est_PDF_surf_plot.ZData = reshape(Phi_hat_x_reg(:,:,j), length(x_2), length(x_1));
    % Actualizar punto inicial de trayectoria
    X_e_0_plot.XData = X_e_reg(1,1,j);
    X_e_0_plot.YData = X_e_reg(1,2,j);

    for i = 1:length(t)
        % Delay (if needed)
        % pause(0.0001)

        % Update title to have the evolution of Time
        tiempo_act = t(i) + t(end)*(j - 1);
        title(plane_tile,"Time: " + num2str(tiempo_act, "%.2f") + " sec",...
            "interpreter", "latex", "FontSize", 26)
        
        % Actualizar distribucion empirica en el plano
        % Emp_dist_surf_plot.ZData = reshape(C_x_reg(:,i,j), length(x_2), length(x_1));
        
        % Aumentar los datos de trayectoria
        if i == 1 % Si es el inicio, limpia la trayectoria y grafica solo el punto inicial
            trayectoria_plot.XData = X_e_reg(i,1,j);
            trayectoria_plot.YData = X_e_reg(i,2,j);
        else
            trayectoria_plot.XData(end + 1) = X_e_reg(i,1,j);
            trayectoria_plot.YData(end + 1) = X_e_reg(i,2,j);
        end
        
        % Actualizar la posición actual
        X_e_act_plot.XData = X_e_reg(i,1,j);
        X_e_act_plot.YData = X_e_reg(i,2,j);
        % Actualizar distribucion empirica
        Emp_dist_surf_plot.ZData = reshape(C_x_reg(:,i,j), length(x_2), length(x_1));

        % Graficar estimación de defectos encontrados en la iteración actual
        if ~isempty(Estim_sol(j).Mu_found) && (i == length(t))
            n_def_found = size(Estim_sol(j).Sigma_found, 3);
            nbDrawingSeg = 100;
            tmp_vec = linspace(-pi, pi, nbDrawingSeg)';
            SD_DefFound = zeros(size(Estim_sol(j).Sigma_found));
            Sigma_ast_DefFound = zeros(size(Estim_sol(j).Sigma_found));
            Ellipse_DefFound = zeros(height(tmp_vec), 2, n_def_found);
            for r = 1:n_def_found
                SD_DefFound(:,:,r) = sqrtm(Estim_sol(j).Sigma_found(:,:,r));
                Sigma_ast_DefFound(:,:,r) = 3*SD_DefFound(:,:,r);
                Ellipse_DefFound(:,:,r) = [cos(tmp_vec), sin(tmp_vec)] * real(Sigma_ast_DefFound(:,:,r)) +...
                        repmat(Estim_sol(j).Mu_found(r,:), nbDrawingSeg, 1);
            end
            figh.CurrentAxes = plane_tile;
            hold on
            for r = 1:n_def_found
                p1 = patch(Ellipse_DefFound(:,1,r), Ellipse_DefFound(:,2,r), FoundDef_color,...
                        'LineWidth', 1.5, 'EdgeColor', FoundDef_color, "FaceAlpha",0.2);
                plot(Estim_sol(j).Mu_found(r,1), Estim_sol(j).Mu_found(r,2),...
                    '+','LineWidth', 1.5, 'color', FoundDef_color);
            end
            hold off
            legend([Est_PDF_plot, RealDef_plot, trayectoria_plot, X_e_0_plot, p1],...
                   {"PDF Estimation, $\hat{\Phi}(\mathbf{x})$",...
                   "Real defects",...
                   "$\mathbf{X_e}(t)$",...
                   "$\mathbf{X_e}(0)$",...
                   "Found Defects"},...
                   'Location','northeast')
        end
        
        % Grafica los defectos No encontrados pero estimados, si los hay
        if (j == n_iter) && (i == length(t)) && ~Estim_sol(end).flag_done
            figh.CurrentAxes = plane_tile;
            hold on
            for r = 1:n_def_not_found
                p2 = patch(Ellipse_DefNotFound(:,1,r), Ellipse_DefNotFound(:,2,r),...
                    NotFoundDef_color,'LineWidth', 1.5, 'EdgeColor',...
                    NotFoundDef_color, "FaceAlpha",0.2);
            end
            plot(Mu_not_found(:,1), Mu_not_found(:,2), '+',...
                'LineWidth', 1.5, 'color', NotFoundDef_color);
            hold off
            legend([Est_PDF_plot, RealDef_plot, trayectoria_plot,...
                   X_e_0_plot, p1, p2],...
                   {"PDF Estimation, $\hat{\Phi}(\mathbf{x})$",...
                   "Real defects",...
                   "$\mathbf{X_e}(t)$",...
                   "$\mathbf{X_e}(0)$",...
                   "Found Defects",...
                   "Not Found Defects"},...
                   'Location','northeast')
        end
        
        % Grafica el último frame para Phi_hat como zeros en caso de  tener 
        % todos los defectos encontrados
        if (j == n_iter) && (i == length(t)) && Estim_sol(end).flag_done
            Est_PDF_surf_plot.ZData = reshape(zeros(size(Phi_x)), length(x_2), length(x_1));
        end
        
        % Guardar frames
        MovieVector(i + length(t)*(j-1)) = getframe(figh);
    
    end

end 

% Para reproducir "n" número de veces vez (n = 1)
% movie(MovieVector, 1, 60)

%% Crear y guardar animación

% *Solution for "All 'cdata' fields in FRAMES must be the same size" error
% *Use the following custom function to resize
% MovieVector = MakeMovieVectorFramesSameSize(MovieVector);

% *Save animation
% myWriter = VideoWriter("Animacion", "Motion JPEG AVI");
% myWriter.FrameRate = 60;
% myWriter.Quality = 95;
% open(myWriter);
% writeVideo(myWriter, MovieVector);
% close(myWriter);