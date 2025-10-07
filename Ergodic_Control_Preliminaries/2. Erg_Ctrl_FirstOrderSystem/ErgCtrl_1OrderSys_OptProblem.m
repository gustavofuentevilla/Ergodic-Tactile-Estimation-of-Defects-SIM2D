% Ejecute primero estas dos lineas en consola
% import casadi.*
% M = Function.load('M.casadi');

close all
clearvars -except M
clc


%% Parámetros del espacio de búsqueda U = [L_1_l, L_1_u] \times [L_2_l, L_2_u]
L_1_l = 0;
dx_1 = 0.01;
L_1_u = 1;

L_2_l = 0;
dx_2 = 0.01;
L_2_u = 1;

% Dimensiones \mathbf{x} = [x_1 x_2]^T
x_1 = (L_1_l:dx_1:L_1_u)';
x_2 = (L_2_l:dx_2:L_2_u)';

%vector de límites inferior y superiores de las dimensiones
L_i_l = [L_1_l, L_2_l];
L_i_u = [L_1_u, L_2_u];

%% Gaussian Mixture distribution (PDF de REFERENCIA)

%Medias de los Gaussianos
mu_1 = [0.5, 0.7]; 
mu_2 = [0.6, 0.3];

Mu = [mu_1; mu_2];

% Matrices de Covarianza
Cov_1 = [0.0500, 0.0150;
         0.0150, 0.0100];
Cov_2 = [0.0130, 0.0060;
         0.0060, 0.0220];

Sigma = cat(3,Cov_1,Cov_2);

%Pesos sobre cada Gaussiano
proporciones = [0.5, 0.5];

gm_dist = gmdistribution(Mu, Sigma, proporciones);

[x_1_grid, x_2_grid] = meshgrid(x_1, x_2);

%Espacio de búsqueda discretizado
Omega = [reshape(x_1_grid,[],1), reshape(x_2_grid,[],1)]; 

%PDF de referencia
Phi_x = pdf(gm_dist, Omega);

%% Cálculo de los coeficientes de Fourier para la PDF de referencia

% Coeficientes por dimensión
K = 8;

% Conjunto de valores para k_i
k_1 = (0:K-1)';
k_2 = (0:K-1)';

[k_1_grid, k_2_grid] = meshgrid(k_1, k_2);

% Conjunto de vectores índice
K_cal = [reshape(k_1_grid,1,[]); reshape(k_2_grid,1,[])];

% Registro para guardar los coeficientes
phi_k_reg = zeros(size(K_cal,2), 1);

% Registro para guardar las funciones de Fourier ortonormales
f_k_reg = zeros(height(Omega), size(K_cal,2));

% Registro para el término normalizador
h_k_reg = zeros(size(K_cal,2),1);

% Cálculo de los coeficientes
for j = 1:size(K_cal,2)
    
    k_vect_j = K_cal(:,j)';

    % función ortogonal de Fourier
    f_tilde_k_j = prod(cos( k_vect_j.*pi.*(Omega - L_i_l)./(L_i_u - L_i_l) ), 2);

    % Término normalizador = escalar por cada k
    h_k_j = sqrt(sum(f_tilde_k_j .^2) * dx_1 * dx_2);

    % función ortonormal de Fourier
    f_k_j = f_tilde_k_j ./ h_k_j;

    % Coeficientes de Fourier, la integral se aproxima con Riemann
    phi_k_j = sum(Phi_x .* f_k_j)*dx_1*dx_2;

    % Registros
    h_k_reg(j) = h_k_j;
    f_k_reg(:,j) = f_k_j;
    phi_k_reg(j) = phi_k_j;

end

%% Reconstrucción de la PDF de referencia para verificación de coeficientes

% Phi_x_reconstructed = zeros(height(Omega),1);
% 
% for i = 1:height(phi_k_reg)
%     Phi_x_reconstructed = Phi_x_reconstructed + phi_k_reg(i)*f_k_reg(:,i);
% end
% 
% % Graficación
% 
% figure(1)
% subplot(2,2,1);
% surf(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)))
% xlim([L_1_l, L_1_u])
% ylim([L_2_l, L_2_u])
% title("PDF de Referencia",'Interpreter','latex')
% xlabel('$x_1$','Interpreter','latex')
% ylabel('$x_2$','Interpreter','latex')
% zlabel('$\Phi(\mathbf{x})$','Interpreter','latex')
% %axis equal
% grid on
% 
% subplot(2,2,2);
% contour(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)))
% xlim([L_1_l, L_1_u])
% ylim([L_2_l, L_2_u])
% title("PDF de Referencia",'Interpreter','latex')
% xlabel('$x_1$','Interpreter','latex')
% ylabel('$x_2$','Interpreter','latex')
% axis equal
% grid on
% 
% subplot(2,2,3);
% surf(x_1_grid, x_2_grid, reshape(Phi_x_reconstructed, length(x_2), length(x_1)))
% xlim([L_1_l, L_1_u])
% ylim([L_2_l, L_2_u])
% title("Reconstrucci\'on",'Interpreter','latex')
% xlabel('$x_1$','Interpreter','latex')
% ylabel('$x_2$','Interpreter','latex')
% zlabel('$\Phi(\mathbf{x})$','Interpreter','latex')
% %axis equal
% grid on
% 
% subplot(2,2,4);
% contour(x_1_grid, x_2_grid, reshape(Phi_x_reconstructed, length(x_2), length(x_1)))
% xlim([L_1_l, L_1_u])
% ylim([L_2_l, L_2_u])
% title("Reconstrucci\'on",'Interpreter','latex')
% xlabel('$x_1$','Interpreter','latex')
% ylabel('$x_2$','Interpreter','latex')
% axis equal
% grid on


%% Condiciones Iniciales y parámetros
T_s = 0.01;                  % Tiempo de muestreo
t_f = 10;     % Tiempo de simulación
t = (0:T_s:t_f)';

n = 2; % Número de dimensiones espaciales

N = 250; % Horizonte de predicción %250

% Peso sobre controles
R = [1e-8, 0;
     0, 1e-8];

% Peso sobre métrica ergódica
gamma = 1;

% Estado inicial z = [z_1; z_2] = [x_1; x_2]
z_0 = [0.1; 0.8]; 

% Velocidad deseada
u_d = [0; 0];

% Coeficientes iniciales
X_e_0 = [z_0(1), z_0(2)];
f_k_traj_0 = prod(cos( K_cal'.*pi.*(X_e_0 - L_i_l)./(L_i_u - L_i_l) ), 2) ./ h_k_reg ;
c_k_0 = (f_k_traj_0*T_s)/t_f ;

% Ergodic metric inicial
p = 2; %norma 2
Lambda_k = (1 + vecnorm(K_cal, p, 1)').^(-(n + 1)/2);
Varepsilon_0 = sum( Lambda_k .* (c_k_0 - phi_k_reg).^2 );

%% Problem Setup

% import casadi.*
% 
% %Ecuaciones x_1_dot = u_1;    x_2_dot = u_2;
% z = MX.sym('z', 2); %states z = [z(1); z(2)] = [x_1; x_2]
% u = MX.sym('u', 2); %controls u = [u(1); u(2)]
% 
% %ODE construction: z_dot = [u(1); u(2)]
% z_dot = [u(1); u(2)];
% f = Function('f', {z,u}, {z_dot}, {'z', 'u'}, {'z_dot'});
% 
% %DAE problem structure
% intg_options = struct;
% intg_options.tf = T_s;   %Integration time (one step ahead)
% intg_options.simplify = true;
% intg_options.number_of_finite_elements = 4; %intermediate steps on the integration horizon
% 
% dae = struct;
% 
% dae.x = z;  % states (formalized)
% dae.p = u;  % parameter, fixed during integration horizon (just one step ahead)
% dae.ode = f(z,u); % symbolic dynamics
% 
% intg = integrator('intg', 'rk', dae, intg_options); %RK4 integration method
% 
% %One step integration (symbolic)
% res = intg('x0', z, 'p', u);
% z_next = res.xf;
% 
% F = Function('F', {z,u}, {z_next}, {'z','u'}, {'z_next'});
% 
% %% Multiple Shooting for one prediction horizon with N+1 samples
% opti = casadi.Opti();
% 
% z = opti.variable(2, N+1);
% u = opti.variable(2, N);
% z_0_sym = opti.parameter(2, 1);  % parameter (not optimized over): initial condition
% phi_k_sym = opti.parameter(size(K_cal,2), 1);
% 
% % Symbolic Fourier functions, coefficients and ergodic metric with casadi
% X_e_sym = [z(1,:)', z(2,:)'];     %Position [x_1, x_2] for all N samples
% 
% c_k_sym = 0;
% J = 0;
% f_k_traj_sym = opti.variable(size(K_cal,2),1);
% for i = 1:N
%     for j = 1:size(K_cal,2)
%         %problems using(.*) with casadi when K_cal is a matrix
%         temp = cos(K_cal(:,j)'.*pi.*(X_e_sym(i,:) - L_i_l)./(L_i_u - L_i_l));   
%         %problems using prod() function
%         %f_k_traj_sym_j = temp(1)*temp(2)/h_k_reg(j);
%         % f_tilde_k_sym_j = temp(1)*temp(2);
%         % h_k_traj_i = sqrt(sum(f_tilde_k_sym_j .^2) * dx_1 * dx_2);
%         % f_k_traj_sym_j = f_tilde_k_sym_j ./ h_k_traj_i;
%         f_k_traj_sym(j,1) = temp(1)*temp(2)/h_k_reg(j); 
%     end
%     c_k_sym = c_k_sym + (f_k_traj_sym*T_s)/(N*T_s);%/i*T_s /t_f
%     Varepsilon_sym = sum( Lambda_k.*(c_k_sym - phi_k_reg).^2 );
% 
%     % Objetive function
%     J = J + gamma*Varepsilon_sym + (u(:,i) - u_d)'*R*(u(:,i) - u_d);
% end 
% 
% opti.minimize( J );
% 
% % Equality Constraints
% for k = 1:N
%     opti.subject_to( z(:,k+1) == F( z(:,k),u(:,k) ) );
% end
% opti.subject_to( z(:,1) == z_0_sym );
% 
% % Inequality Constraints
% opti.subject_to( -100 <= u <= 100 );
% opti.subject_to( L_1_l <= z(1,:) <= L_1_u );    % x_1 boundaries
% opti.subject_to( L_1_l <= z(2,:) <= L_1_u );    % x_2 boundaries
% 
% % Solver definition
% solver_opts = struct;
% solver_opts.expand = true;
% solver_opts.print_time = false;
% solver_opts.ipopt.max_iter = 2000;
% %solver_opts.ipopt.print_level = 0;
% solver_opts.ipopt.acceptable_tol = 1e-8;
% solver_opts.ipopt.acceptable_obj_change_tol = 1e-6;
% 
% opti.solver('ipopt', solver_opts);
% 
% opti.set_value(z_0_sym, z_0); 
% opti.set_value(phi_k_sym, phi_k_reg);
% 
% solution = opti.solve();
% M = opti.to_function('M', {z_0_sym, phi_k_sym}, {z, u}, {'z_0','phi_k'}, {'z','u'});

%% Saving and loading casadi function object

% M.save('M.casadi');

% M = Function.load('M.casadi');

%% Charts: One solution for horizon with N steps

t_1 = (0:T_s:N*T_s)';

[Z, U] = M(z_0, phi_k_reg); %Soluciones
estados = full(Z)';
u_tmp = full(U)';
controles = [u_tmp; [NaN, NaN]];

metrica = zeros(length(t_1), 1);
Varepsilon_temp = Varepsilon_0;
c_k_temp = c_k_0;

C_x = zeros(height(Omega), 1);
C_x_reg = zeros(height(Omega), length(t_1));

for i = 1:length(t_1)

    metrica(i) =  Varepsilon_temp;

    % Compute Fourier Functions and coefficients on the new position
    z_temp = estados(i,:);
    X_e_temp = [z_temp(1), z_temp(2)];   % Position [x_1, x_2]
    f_k_traj_temp = prod(cos( K_cal'.*pi.*(X_e_temp - L_i_l)./(L_i_u - L_i_l) ), 2) ./ h_k_reg ;
    % f_tilde_k_temp_j = prod(cos( K_cal'.*pi.*(X_e_temp - L_i_l)./(L_i_u - L_i_l) ), 2);
    % h_k_traj_i = sqrt(sum(f_tilde_k_temp_j .^2) * dx_1 * dx_2);
    % f_k_traj_temp = f_tilde_k_temp_j ./ h_k_traj_i;
    c_k_temp = c_k_temp + (f_k_traj_temp*T_s)/(N*T_s) ; %/t_f t_1(i) + T_s

    %Ergodic metric
    Varepsilon_temp = sum( Lambda_k .* (c_k_temp - phi_k_reg).^2 ); 
    
    %Empirical distribution reconstruction
    C_x_i = zeros(height(Omega), 1);
    for j = 1:size(K_cal, 2)
        C_x_i = C_x_i + c_k_temp(j)*f_k_reg(:,j);
    end
    
    %Se suman todas las distribuciones generadas en cada muestra
    C_x = C_x + C_x_i;

    %Se registra
    C_x_reg(:,i) = C_x;

end

%% Save workspace variables

save("Animacion/data_for_animation.mat", "-regexp", "^(?!(M)$).");

%% Gráficas
fig1h = figure(1);
subplot(2,1,1)
plot(t_1, estados, 'LineWidth', 2)
title("Position states")
xlabel('Time [s]')
ylabel('Position [m]')
legend('$x_1$', '$x_2$')
grid on
subplot(2,1,2)
stairs(t_1, controles, 'LineWidth', 2)
title("Control actions")
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('$u_1$', '$u_2$')
grid on

set(findall(fig1h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig1h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig1h, "-property", "FontSize"), "FontSize", 20)

fig2h = figure(2);
subplot(2,2,1);
pcolor(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "FaceColor","interp","EdgeColor","none")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Reference PDF")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
axis equal tight
grid on
legend('$\Phi(\mathbf{x})$','Location','northeastoutside')
subplot(2,2,2);
pcolor(x_1_grid, x_2_grid, reshape(C_x, length(x_2), length(x_1)),...
    "FaceColor","interp","EdgeColor","none")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Empirical distribution reconstruction")
xlabel('$x_1$')
ylabel('$x_2$')
axis equal tight
grid on
hold on
plot(estados(:,1),estados(:,2),"Color", "black", 'LineWidth',2)
plot(estados(1,1),estados(1,2),'ksq','MarkerSize',15,'LineWidth',2)
legend('$C(\mathbf{x})$','$\mathbf{X_e}(t)$', '$\mathbf{X_e}(0)$',...
    'Location','northeastoutside')
hold off
subplot(2,2,3);
surf(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)))
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$\Phi(\mathbf{x})$')
grid on
subplot(2,2,4);
surf(x_1_grid, x_2_grid, reshape(C_x, length(x_2), length(x_1)))
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$C(\mathbf{x})$')
grid on

set(findall(fig2h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig2h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig2h, "-property", "FontSize"), "FontSize", 20)

fig3h = figure(3);
plot(t_1, metrica, "LineWidth",2)
title("Ergodic Metric")
xlabel('Time [s]')
ylabel('$\varepsilon \left( \mathbf{X_e}(t), \Phi(\mathbf{x}) \right) $')
grid on

set(findall(fig3h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig3h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig3h, "-property", "FontSize"), "FontSize", 20)



