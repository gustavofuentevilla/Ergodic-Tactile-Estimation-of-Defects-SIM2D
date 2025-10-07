%% Gráficas animadas (descomente el comando pause() para visualizar)

load("data_for_animation");

figh = figure;
figh.Units = "centimeters";
figh.Position = [5, 5, 40, 22];
figh.Units = "pixels";
% fig4h.Visible = "off"; %Para no mostrar la animación durante el For

layout4h = tiledlayout(figh,1,2);
title(layout4h, "Reconstrucci\'on de la distribuci\'on emp\'irica",...
    "interpreter", "latex", "FontSize", 20) %El comando set de abajo no afecta este título

% Primera iteración
superficie_tile = nexttile(layout4h);
surf_plot = surf(x_1_grid, x_2_grid, reshape(C_x_reg(:,1), length(x_2), length(x_1)),...
    "EdgeColor","none", "FaceColor","interp");
view([-38 63])
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$C(\mathbf{x})$')
axis tight
shading flat
grid on

plano_tile = nexttile(layout4h);
% using pcolor()
% plane_plot = pcolor(x_1_grid, x_2_grid, reshape(C_x_reg(:,1), length(x_2), length(x_1)),...
%     "EdgeColor","none", "FaceColor","interp"); 
% using contour()
[~, plane_plot] = contour(x_1_grid, x_2_grid, reshape(C_x_reg(:,1), length(x_2), length(x_1)),...
    "LineWidth", 2);
% using imagesc()
% plane_plot = imagesc([L_1_l, L_1_u], [L_2_l, L_2_u], reshape(C_x_reg(:,1), length(x_2), length(x_1)));
axis([L_1_l, L_1_u, L_2_l, L_2_u])
xlabel('$x_1$')
ylabel('$x_2$')
axis equal tight
axis xy
grid on
hold on
X_init_plot4 = plot(x_t(1), y_t(1), "sk", "LineWidth", 3, "MarkerSize", 18);
X_t_plot4 = plot(x_t(1), y_t(1), "LineWidth", 3, "Color", "black");
X_circle_plot4 = plot(x_t(1), y_t(1), "ok", "LineWidth", 3, "MarkerSize", 18);
hold off

legend([X_init_plot4, X_t_plot4], {"$\mathbf{X}_e(0)$", "$\mathbf{X}_e(t)$"})
    %"Location", "best")

set(findall(figh,'-property','Interpreter'),'Interpreter','latex') 
set(findall(figh,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(figh, "-property", "FontSize"), "FontSize", 20)

%Frames Pre-allocation
loops = length(t);
MovieVector(loops) = struct("cdata", [], "colormap", []);

MovieVector(1) = getframe(figh);

for i = 2:length(t)
    % Delay (if needed)
    % pause(0.0001)
    
    % Actualizar surf
    surf_plot.ZData = reshape(C_x_reg(:,i), length(x_2), length(x_1));

    % Actualizar distribución en el plano con pcolor()
    % note que la actualización es en CData y no en ZData 
    % (ZData es constante, pcolor es una vista planar de surf)
    % plane_plot.CData = reshape(C_x_reg(:,i), length(x_2), length(x_1));

    % con contour()
    plane_plot.ZData = reshape(C_x_reg(:,i), length(x_2), length(x_1));

    % con imagesc()
    % plane_plot.CData = reshape(C_x_reg(:,i), length(x_2), length(x_1));
    
    % Aumentar los datos de trayectoria y actualizar la posición actual
    X_t_plot4.XData(end + 1) = x_t(i);
    X_t_plot4.YData(end + 1) = y_t(i);
    X_circle_plot4.XData = x_t(i);
    X_circle_plot4.YData = y_t(i);

    MovieVector(i) = getframe(figh);

end

% fig4h.Visible = "on";

% Para reproducir "n" número de veces vez (n = 1)
% movie(MovieVector, 1, 60)

%% Crear y guardar animación

% Solution for "All 'cdata' fields in FRAMES must be the same size" error
% Use the following custom function to resize
% MovieVector = MakeMovieVectorFramesSameSize(MovieVector);

%% Save animation
% myWriter = VideoWriter("Animacion");
% myWriter.FrameRate = 60;
% myWriter.Quality = 95;
% open(myWriter);
% writeVideo(myWriter, MovieVector);
% close(myWriter);