close all
clear
clc

%% Espacio de búsqueda y definición de la distribución de referencia

n = 2; % Número de dimensiones

%Parámetros del espacio de búsqueda Omega = [L_1_l, L_1_u] \times [L_2_l, L_2_u]
% con \mathbf{x} = [x_1 x_2]^T
L_1_l = 0.5;
dx_1 = 0.01;
L_1_u = 1.5;

L_2_l = 0.5;
dx_2 = 0.01;
L_2_u = 1.5;

x_1 = (L_1_l:dx_1:L_1_u)';
x_2 = (L_2_l:dx_2:L_2_u)';

%x_vect = [x_1, x_2];

%vector de límites inferior y superiores de las dimensiones
L_i_l = [L_1_l, L_2_l];
L_i_u = [L_1_u, L_2_u];

% Gaussian Mixture distribution (PDF de REFERENCIA)
mu_1 = [0.9, 0.9]; %Media del primer Gaussiano
mu_2 = [0.8, 1.2];
mu_3 = [1.3, 1];

Mu = [mu_1; mu_2; mu_3];

Cov_1 = [0.01, 0.004;
         0.004, 0.01];
Cov_2 = [0.005, -0.003;
         -0.003, 0.005];
Cov_3 = [0.008, 0.0;
         0.0, 0.004];

Sigma = cat(3,Cov_1,Cov_2,Cov_3);

proporciones = [0.5, 0.2, 0.3]; % Pesos sobre cada Gaussiano (la suma de los pesos debe dar 1)

gm_dist = gmdistribution(Mu, Sigma, proporciones); %Genera el objeto de la distribución

[x_1_grid, x_2_grid] = meshgrid(x_1, x_2);

Omega = [reshape(x_1_grid,[],1), reshape(x_2_grid,[],1)]; %Describe todos las coordenadas en el espacio de búsqueda

Phi_x = pdf(gm_dist, Omega);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%Cálculo de los coeficientes de Fourier para la PDF de referencia%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 25; %Número de coeficientes por dimensión

%Al final, se tendrán K^n coeficientes en total para toda la función de
%referencia \Phi(\mathbf{x})

%Cálculo de los valores que puede tomar k_i
k_1 = (0:K-1)';
k_2 = (0:K-1)';

[k_1_grid, k_2_grid] = meshgrid(k_1, k_2);

%K_cal es el conjunto de todas las ^j\mathbf{k}, con (K^n) elementos
K_cal = [reshape(k_1_grid,1,[]); reshape(k_2_grid,1,[])];

%% Cálculo de coeficientes usando un For
%
% %Registro para guardar los coeficientes
% phi_k_reg = zeros(size(K_cal,2), 1);
% 
% %Registro para guardar las funciones evaluadas en \mathbf{x}, cada columna
% %corresponde a la j-ésima combinacion de \mathbf{k}. Las filas corresponden
% %a cada punto discreto en el espacio, y las columnas a la combinación de k.
% f_k_reg = zeros(height(U), size(K_cal,2));
% 
% %Registro para el término normalizador (h_{\mathbf{k}})
% h_k_reg = zeros(size(K_cal,2),1);
% 
% % Este for recorre la j-ésima combinación para k
% for j = 1:size(K_cal,2)
% 
%     k_vect_j = K_cal(:,j)';
% 
%     % función ortogonal \tilde{f}_{\mathbf{k}}
%     f_tilde_k_j = prod(cos( k_vect_j.*pi.*(U - L_i_l)./(L_i_u - L_i_l) ), 2);
% 
%     % Término normalizador (h_{\mathbf{k}}) = un escalar por cada j-ésima
%     % combinación de \mathbf{k}
%     h_k_j = sqrt(sum(f_tilde_k_j .^2) * dx_1 * dx_2);
% 
%     % función ortonormal f_{\mathbf{k}} correspondiente a la j-ésima
%     % combinación del vector k (^j \mathbf{k}), evaluada en todo el espacio U.
%     f_k_j = f_tilde_k_j ./ h_k_j;
% 
%     % Coeficientes \phi_{^j \mathbf{k}}, la integral se aproxima con Riemann
%     phi_k_j = sum(Phi_x .* f_k_j)*dx_1*dx_2;
% 
%     % Se guardan los coeficientes, las funciones ortonormales y los
%     % términos normalizadores en registros
%     h_k_reg(j) = h_k_j;
%     f_k_reg(:,j) = f_k_j;
%     phi_k_reg(j) = phi_k_j;
% 
% end

%%%%%%%% Calculo de coeficientes SIN un For loop %%%%%%%%%
f_tilde_k = prod(cos(repmat(reshape(K_cal,1,n,K^n), [height(Omega), 1, 1]).*...
    repmat(pi.*(Omega - L_i_l)./(L_i_u - L_i_l), [1, 1, K^n])), 2);
h_k = sqrt(sum(f_tilde_k.^2) * dx_1 * dx_2);
f_k = f_tilde_k ./ repmat(h_k, [height(Omega),1,1]);
phi_k = sum( repmat(Phi_x,[1,1,K^n]) .* f_k)*dx_1*dx_2;

%Registros
h_Kcal = reshape(h_k, K^n, 1);
f_Kcal = reshape(f_k, height(Omega), K^n);
phi_Kcal = reshape(phi_k, K^n, 1);
h_k_reg = h_Kcal;
f_k_reg = f_Kcal;
phi_k_reg = phi_Kcal;

%% Reconstrucción de la PDF de referencia para verificación de coeficientes

Phi_x_reconstructed = zeros(height(Omega),1);

for i = 1:height(phi_k_reg)
    Phi_x_reconstructed = Phi_x_reconstructed + phi_k_reg(i)*f_k_reg(:,i);
end

%% Trayectoria aleatoria (hypotrochoid)
      
T = 0.01;                  % Paso
t_f = 10;     % Tiempo de simulación
t = (0:T:t_f)';

D = 0.3;
R = D;
r = 3*D/5;
x_c = (L_1_u - L_1_l)/2 + L_1_l;
y_c = (L_2_u - L_2_l)/2 + L_2_l;

x_t = x_c + (R - r)*cos(pi*t) + D*cos((R - r)/r * pi*t);
y_t = y_c + (R - r)*sin(pi*t) - D*sin((R - r)/r * pi*t);

%Trayectoria Circulo

% r_d = 0.4;
% x_c = (L_1_u - L_1_l)/2 + L_1_l;
% y_c = (L_2_u - L_2_l)/2 + L_2_l;
%
% x_t = x_c + r_d*cos(delta_t);
% y_t = y_c + r_d*sin(delta_t);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%Cálculo de los coeficientes de Fourier de la distribución empírica que
%genera la trayectoria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Registro para guardar los valores de las funciones de Fourier sobre la
%trayectoria, las filas corresponden a las combinaciones de \mathbf{k}, y
%las columnas corresponden a los pasos temporales
f_k_traj_reg = zeros(size(K_cal, 2), length(t));

%Registro para guardar los vectores de coeficientes en todo t
c_k_reg = zeros(size(K_cal, 2), length(t));

%Inicialización de los coeficientes para la trayectoria, se generan los K^n
%coeficientes en cada paso de muestreo t
c_k_i = zeros(size(K_cal, 2), 1);

for i = 1:length(t)

    X_e = [x_t(i), y_t(i)]; %Vector fila, para que los cálculos coincidan en dimensión

    %función ortonormal para el par de puntos actuales (x(t),y(t)) sobre
    %todas las combinaciones de \mathbf{k}
    f_k_traj_i = prod(cos( K_cal'.*pi.*(X_e - L_i_l)./(L_i_u - L_i_l) ), 2) ./ h_k_reg ;

    %Generación y actualización de los coeficientes: la normalización se
    %hace con respecto al tiempo final de la integral actual (t_act + T)
    c_k_i = c_k_i + (f_k_traj_i*T)/t_f ;

    %Registros
    f_k_traj_reg(:,i) = f_k_traj_i;
    c_k_reg(:,i) = c_k_i;

end

%% Cálculo de los coeficientes de la trayectoria calculando h_k en línea
% 
% h_k_traj_reg = zeros(length(t), 1);
% 
% for i = 1:length(t)
%     X_e = [x_t(i), y_t(i)]; %Vector fila, para que los cálculos coincidan en dimensión
% 
%     % función ortogonal \tilde{f}_{\mathbf{k}} para el par de puntos actuales (x(t),y(t))
%     f_tilde_k_traj_i = prod(cos( K_cal'.*pi.*(X_e - L_i_l)./(L_i_u - L_i_l) ), 2);
% 
%     % Término normalizador (h_{\mathbf{k}}) = un escalar por cada j-ésima
%     % combinación de \mathbf{k}
%     h_k_traj_i = sqrt(sum(f_tilde_k_traj_i .^2) * dx_1 * dx_2);
% 
%     % función ortonormal f_{\mathbf{k}} correspondiente a la j-ésima
%     % combinación del vector k (^j \mathbf{k}), evaluada en todo el espacio U.
%     f_k_traj_i = f_tilde_k_traj_i ./ h_k_traj_i;
% 
%     %Generación y actualización de los coeficientes: la normalización se
%     %hace con respecto al tiempo final de la integral actual (t_act + T)
%     c_k_i = c_k_i + (f_k_traj_i*T)/t_f ;
% 
%     %Registros
%     f_k_traj_reg(:,i) = f_k_traj_i;
%     c_k_reg(:,i) = c_k_i;
%     h_k_traj_reg(i) = h_k_traj_i;
% end

%% Reconstrucción de la distribución empírica que genera la trayectoria C(\mathbf{x})

C_x = zeros(height(Omega), 1);
C_x_reg = zeros(height(Omega), length(t));

for i = 1:length(t)
    
    %Genera la distribución en t(i), en total se tienen tantas
    %distribuciones como muestras de tiempo
    C_x_i = zeros(height(Omega), 1);
    for j = 1:size(K_cal, 2)
        C_x_i = C_x_i + c_k_reg(j,i)*f_k_reg(:,j);
    end
    
    %Se suman todas las distribuciones generadas en cada muestra
    C_x = C_x + C_x_i;

    %Se registra
    C_x_reg(:,i) = C_x;

end

%% Ergodic Metric

%También se puede generar justo después de ir generando los coeficientes
%para la trayectoria empírica, en el For de arriba

Epsilon_t_reg = zeros(length(t), 1);

%Cálculo de todas las Lambdas_{\mathbf{k}}
p = 2; %norma 2
Lambda_k = (1 + vecnorm(K_cal, p, 1)').^(-(n + 1)/2);

for i = 1:length(t)

    %Métrica calculada en cada t
    Epsilon_t = sum( Lambda_k .* (c_k_reg(:,i) - phi_k_reg).^2 ); 

    %Registro de la métrica en todo el historial de tiempo
    Epsilon_t_reg(i) =  Epsilon_t;

end

%% Save workspace variables for animation

save("Animacion/data_for_animation.mat")


%% Gráficas

% Distribución de referencia

% Distribución de Referencia

fig1h = figure(1);
layout1h = tiledlayout(fig1h,2,2);
title(layout1h, "Reconstrucci\'on de PDF con funciones base de Fourier",...
    "interpreter", "latex", "FontSize", 18)
colormap default

nexttile(layout1h)
surf(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "EdgeColor","none", "FaceColor","interp")
axis([L_1_l, L_1_u, L_2_l, L_2_u])
title("PDF de Referencia",'Interpreter','latex')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
zlabel('$\Phi(\mathbf{x})$','Interpreter','latex')
grid on
axis tight

nexttile(layout1h)
surf(x_1_grid, x_2_grid, reshape(Phi_x_reconstructed, length(x_2), length(x_1)),...
    "EdgeColor","none", "FaceColor","interp")
axis([L_1_l, L_1_u, L_2_l, L_2_u])
title("Reconstrucci\'on")
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$\Phi(\mathbf{x})$')
axis tight
grid on

nexttile(layout1h)
pcolor(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "EdgeColor","none", "FaceColor","interp")
axis([L_1_l, L_1_u, L_2_l, L_2_u])
xlabel('$x_1$')
ylabel('$x_2$')
axis equal tight
grid on

nexttile(layout1h)
pcolor(x_1_grid, x_2_grid, reshape(Phi_x_reconstructed, length(x_2), length(x_1)),...
    "EdgeColor","none", "FaceColor","interp")
axis([L_1_l, L_1_u, L_2_l, L_2_u, 0, 10])
xlabel('$x_1$')
ylabel('$x_2$')
axis equal tight
grid on

set(findall(fig1h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig1h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig1h, "-property", "FontSize"), "FontSize", 18)

%% Distribución empírica y trayectoria

fig2h = figure(2);
layout2h = tiledlayout(fig2h,1,2);
title(layout2h, "Reconstrucci\'on de la distribuci\'on emp\'irica",...
    "interpreter", "latex", "FontSize", 20)

ax_surf = nexttile(layout2h);
surf(x_1_grid, x_2_grid, reshape(C_x, length(x_2), length(x_1)),...
    "EdgeColor","none", "FaceColor","interp")
view([-44 56])
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$C(\mathbf{x})$')
grid on
axis tight
colormap(ax_surf, brewermap([],'PuBuGn'))

ax_plane = nexttile(layout2h);
pcolor(x_1_grid, x_2_grid, reshape(C_x, length(x_2), length(x_1)),...
    "EdgeColor","none", "FaceColor","interp")
axis([L_1_l, L_1_u, L_2_l, L_2_u])
xlabel('$x_1$')
ylabel('$x_2$')
axis equal tight
grid on
colormap(ax_plane, brewermap([],'PuBuGn'))
hold on
X_t_plot = plot(x_t, y_t, "LineWidth", 3, "Color", "black");
X_init_plot = plot(x_t(1), y_t(1), "sk", "LineWidth", 3, "MarkerSize", 20);
hold off
legend([X_init_plot, X_t_plot], {"$\mathbf{X}_e(0)$", "$\mathbf{X}_e(t)$"},...
    "Location", "best")

set(findall(fig2h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig2h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig2h, "-property", "FontSize"), "FontSize", 20)

%% Ergodic Metric

fig3h = figure(3);
plot(t, Epsilon_t_reg, "LineWidth",3)
title("Ergodic Metric",'Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('$\varepsilon \left( \mathbf{X_e}(t), \Phi(\mathbf{x}) \right) $','Interpreter','latex')
grid on

set(findall(fig3h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig3h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig3h, "-property", "FontSize"), "FontSize", 20)
