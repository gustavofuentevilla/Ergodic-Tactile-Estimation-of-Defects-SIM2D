close all
clear
clc

%% Cargar datos
load("data_for_animation.mat");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figh = figure;
figh.Units = "centimeters";
figh.Position = [5, 5, 50, 30];
figh.Units = "pixels";

layouth = tiledlayout(figh,2,3);
colormap default

%Primera iteración
nexttile
pcolor(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "FaceColor","interp","EdgeColor","none")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Reference PDF")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
axis equal tight
grid on
legend("$\Phi(\mathbf{x})$")

nexttile
Est_PDF_plot = pcolor(x_1_grid, x_2_grid, reshape(Phi_hat_x_reg(:,:,1), length(x_2), length(x_1)),...
    "FaceColor","interp","EdgeColor","none");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Estimated PDF")
xlabel('$x_1$')
ylabel('$x_2$')
axis equal tight
grid on
legend("$\hat{\Phi}(\mathbf{x})$")

nexttile
Emp_dist_plot = pcolor(x_1_grid, x_2_grid, reshape(C_x_reg(:,1,1), length(x_2), length(x_1)),...
    "FaceColor","interp","EdgeColor","none");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Empirical distribution reconstruction")
xlabel('$x_1$')
ylabel('$x_2$')
axis equal tight
grid on
hold on
trayectoria_plot = plot(X_e_reg(1,1,1), X_e_reg(1,2,1),"Color", "black", 'LineWidth',2);
X_e_0_plot = plot(X_e_reg(1,1,1), X_e_reg(1,2,1),'ksq','MarkerSize',15,'LineWidth',2);
X_e_act_plot = plot(X_e_reg(1,1,1), X_e_reg(1,2,1),"ok", 'MarkerSize',15,'LineWidth',2);
hold off
legend([Emp_dist_plot, trayectoria_plot, X_e_0_plot],...
    {"$C(\mathbf{x})$", "$\mathbf{X_e}(t)$", "$\mathbf{X_e}(0)$"},...
    'Location','northeastoutside')

nexttile
surf(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "EdgeColor","none", "FaceColor","interp")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$\Phi(\mathbf{x})$')
grid on
axis tight

nexttile
Est_PDF_surf_plot = surf(x_1_grid, x_2_grid, reshape(Phi_hat_x_reg(:,:,1), length(x_2), length(x_1)),...
    "FaceColor","interp","EdgeColor","none");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$\Phi(\mathbf{x})$')
grid on
axis tight

nexttile
Emp_dist_surf_plot = surf(x_1_grid, x_2_grid, reshape(C_x_reg(:,1,1), length(x_2), length(x_1)),...
    "EdgeColor","none", "FaceColor","interp");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$C(\mathbf{x})$')
grid on
axis tight

set(findall(figh,'-property','Interpreter'),'Interpreter','latex') 
set(findall(figh,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(figh, "-property", "FontSize"), "FontSize", 20)

%Frames Pre-allocation
loops = length(t)*n_iter;
MovieVector(loops) = struct("cdata", [], "colormap", []);

for j = 1:n_iter

    title(layouth, "Single Defect Localization, iteration " + j,...
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
        
        % Actualizar distribucion empirica en el plano
        Emp_dist_plot.CData = reshape(C_x_reg(:,i,j), length(x_2), length(x_1));
        
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
        Emp_dist_surf_plot.ZData = reshape(C_x_reg(:,i,j), length(x_2), length(x_1));
        
        % Guardar frames
        MovieVector(i + length(t)*(j-1)) = getframe(figh);
    
    end

end 

% Para reproducir "n" número de veces vez (n = 1)
% movie(MovieVector, 1, 60)

%% Crear y guardar animación

% Solution for "All 'cdata' fields in FRAMES must be the same size" error
% Use the following custom function to resize
% MovieVector = MakeMovieVectorFramesSameSize(MovieVector);

% Save animation
myWriter = VideoWriter("Animacion.avi");
myWriter.FrameRate = 30;
myWriter.Quality = 95;
open(myWriter);
writeVideo(myWriter, MovieVector);
close(myWriter);