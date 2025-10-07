% Ejecute primero estas dos lineas en consola
% import casadi.*
% M = Function.load('/home/gustavo-fuentevilla/MATLAB/Tactile_Defects_Localization/Casadi_Formulation_ExplTask/M_N200.casadi');

close all
clearvars -except M
clc

%% Parámetros del espacio de búsqueda U = [L_1_l, L_1_u] \times [L_2_l, L_2_u]
L_1_l = 0.5;
dx_1 = 0.02;
L_1_u = 1.5;

L_2_l = 0.5;
dx_2 = 0.02;
L_2_u = 1.5;

% Dimensiones \mathbf{x} = [x_1 x_2]^T
x_1 = (L_1_l:dx_1:L_1_u)';
x_2 = (L_2_l:dx_2:L_2_u)';

%vector de límites inferior y superiores de las dimensiones
L_i_l = [L_1_l, L_2_l];
L_i_u = [L_1_u, L_2_u];

%% Gaussian Mixture distribution (PDF de REFERENCIA)

%Medias de los Gaussianos
mu_1 = [0.9, 0.9]; 
mu_2 = [0.8, 1.2];
mu_3 = [1.3, 1];

Mu = [mu_1; mu_2; mu_3];

% Matrices de Covarianza
Cov_1 = [0.01, 0.004;
         0.004, 0.01];
Cov_2 = [0.005, -0.003;
         -0.003, 0.005];
Cov_3 = [0.008, 0.0;
         0.0, 0.004];

Sigma = cat(3,Cov_1,Cov_2,Cov_3);

%Pesos sobre cada Gaussiano
proporciones = [0.5, 0.2, 0.3];

gm_dist = gmdistribution(Mu, Sigma, proporciones);

[x_1_grid, x_2_grid] = meshgrid(x_1, x_2);

%Espacio de búsqueda discretizado
Omega = [reshape(x_1_grid,[],1), reshape(x_2_grid,[],1)]; 

%PDF de referencia
Phi_x = pdf(gm_dist, Omega);

%% Cálculo de los coeficientes de Fourier para la PDF de referencia

% Coeficientes por dimensión
K = 12;

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

%% Condiciones Iniciales y parámetros

n = 2; % Número de dimensiones espaciales

N = 200; % Horizonte de predicción

t_f = 10;           %Tiempo final por iteración
T_s = t_f/N;                  % Tiempo de muestreo
t = (0:T_s:t_f)';   %Vector de tiempo por iteración

% Peso sobre controles
R = [7e-5, 0;
     0, 7e-5]*(1/T_s);

% Peso sobre métrica ergódica
gamma = 1;

% Estado inicial z = [z_1; z_2; z_3; z_4] = [x_1; x_1_dot; x_2; x_2_dot]
z_0 = [0.5; 0; 0.5; 0]; 

% Pre-cálculo de Lambda
p = 2; %norma 2
Lambda_k = (1 + vecnorm(K_cal, p, 1)').^(-(n + 1)/2);

%% Problem Setup
% import casadi.*
%
% % Ecuaciones x_1_ddot = u_1;    x_2_ddot = u_2;
% z = MX.sym('z', 4); %states z = [z(1); z(2); z(3); z(4)] = [x_1; x_1_dot; x_2; x_2_dot]
% u = MX.sym('u', 2); %controls u = [u(1); u(2)]
% 
% %ODE construction: z_dot = [z(2); u(1); z(4); u(2)]
% z_dot = [z(2); u(1); z(4); u(2)];
% f = Function('f', {z,u}, {z_dot}, {'z', 'u'}, {'z_dot'});
% 
% %DAE problem structure
% 
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
% %One step integration (numerical)
% % res = intg('x0', z_0, 'p', [0.5; 0]);  %z_0 initial condition, p = u controls
% % z_next = res.xf;
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
% z = opti.variable(4, N+1);
% u = opti.variable(2, N);
% z_0_sym = opti.parameter(4, 1);  % parameter (not optimized over): initial condition
% %u_0_sym = opti.parameter(2, 1);
% phi_k_sym = opti.parameter(K^n, 1);
% %u_d_sym = opti.parameter(2, 1);
% 
% % Symbolic Fourier functions, coefficients and ergodic metric with casadi
% X_e_sym = [z(1,:)', z(3,:)'];     %Position [x_1, x_2] for all N samples
% 
% c_k_sym = 0;
% f_k_traj_sym = opti.variable(K^n,1);
% J = 0;
% for i = 1:N
%     for j = 1:K^n
%         %problems using(.*) with casadi when K_cal is a matrix
%         temp = cos(K_cal(:,j)'.*pi.*(X_e_sym(i,:) - L_i_l)./(L_i_u - L_i_l));   
%         %problems using prod() function
%         f_k_traj_sym(j,1) = temp(1)*temp(2)/h_k_reg(j);
%     end
%     c_k_sym = c_k_sym + (f_k_traj_sym*T_s)/(t_f);%/i*T_s /t_f
%     Varepsilon_sym = sum( Lambda_k.*(c_k_sym - phi_k_sym).^2 );
% 
%     % Objetive function
%     J = J + gamma*Varepsilon_sym + u(:,i)'*R*u(:,i)*T_s; %(u(:,i) - u_d)   u(:,i)'*R*u(:,i)
% end
% 
% % f_k_traj_sym = prod(cos( K_cal'.*pi.*(X_e_sym(i,:) - L_i_l)./(L_i_u - L_i_l) ), 2) ./ h_k_reg ;
% % c_k_sym = c_k_sym + (f_k_traj_act*T_s)/t_f ;
% % Varepsilon_act = sum( Lambda_k .* (c_k_act - phi_k_reg).^2 ); 
% 
% opti.minimize( J );
% 
% % Equality Constraints
% for k = 1:N
%     opti.subject_to( z(:,k+1) == F( z(:,k),u(:,k) ) );
% end
% opti.subject_to( z(1,1) == z_0_sym(1) ); % posiciones iniciales
% opti.subject_to( z(3,1) == z_0_sym(3) );
% opti.subject_to( z(2,1:2) == 0 ); % velocidades iniciales cero
% opti.subject_to( z(4,1:2) == 0 );
% opti.subject_to( z(2,end) == 0 ); % velocidades finales cero
% opti.subject_to( z(4,end) == 0 );
% %opti.subject_to( u(:,1) == u_0_sym ); % controles (aceleraciones) iniciales
% 
% % Inequality Constraints
% opti.subject_to( -50 <= u <= 50 );
% opti.subject_to( L_1_l <= z(1,:) <= L_1_u );    % x_1 boundaries
% opti.subject_to( L_1_l <= z(3,:) <= L_1_u );    % x_2 boundaries
% 
% % Solver definition
% % p_opts = struct('expand',true);
% % s_opts = struct('max_iter',100);
% solver_opts = struct;
% solver_opts.expand = true;
% solver_opts.ipopt.max_iter = 2000;
% opts.ipopt.print_level = 0;
% %opts.ipopt.acceptable_tol = 1e-8;
% %opts.ipopt.acceptable_obj_change_tol = 1e-6;
% 
% opti.solver('ipopt', solver_opts);
% 
% %opti.solver('sqpmethod', struct('qpsol','qrqp'));
% 
% %setting the initial conditions
% opti.set_value(z_0_sym, z_0);
% %opti.set_value(u_0_sym, u_0);
% opti.set_value(phi_k_sym, phi_k_reg);
% %opti.set_value(u_d_sym, u_d);
% 
% %Solving the optimization problem over the horizon N
% solution = opti.solve();
% 
% %Function mapping from initial condition z0 to the optimal control action
% %(first u of the N control actions). M contains IPOPT method embedded and
% %integration method RK4
% M = opti.to_function('M', {z_0_sym, phi_k_sym}, {z, u}, {'z_0','phi_k'}, {'z','u'});
% 
% %M = opti.to_function('M', {z_0_sym, u_0_sym, phi_k_sym}, {z, u}, {'z_0','u_0','phi_k'}, {'z','u'});

%% Saving and loading casadi function object

% M.save('M_N200.casadi');

% M = Function.load('M_N200.casadi');


%% Optimal Problem Solution

% Initializations
z_act = z_0;
phi_k_act = phi_k_reg;

[Z, U] = M(z_act, phi_k_act); %Soluciones
Z = full(Z)';
U = full(U)';
%%
u_plot = [U; [NaN, NaN]];

X_e = [Z(:,1), Z(:,3)]; %X_e = [x_1, x_2] System Position
X_e_dot = [Z(:,2), Z(:,4)];

%% Reconstrucción de distribución empírica y métrica ergódica

Varepsilon_reg = zeros(length(t), 1);

C_x_reg = zeros(height(Omega), length(t));
c_k = zeros(K^n, 1);
C_x = zeros(height(Omega), 1);

for i = 1:length(t)

    % Compute Fourier Functions and coefficients on the new position
    f_k_traj = prod(cos( K_cal'.*pi.*(X_e(i,:) - L_i_l)./(L_i_u - L_i_l) ), 2) ./ h_k_reg ;
    c_k = c_k + (f_k_traj*T_s)/t_f ;

    %Ergodic metric
    Varepsilon = sum( Lambda_k .* (c_k - phi_k_reg).^2 );
    
    %Empirical distribution reconstruction
    C_x_i = zeros(height(Omega), 1);
    for j = 1:K^n
        C_x_i = C_x_i + c_k(j)*f_k_reg(:,j);
    end
    
    %Se suman todas las distribuciones generadas en cada muestra
    C_x = C_x + C_x_i;

    %Se registra
    C_x_reg(:,i) = C_x;
    Varepsilon_reg(i) =  Varepsilon;

end

%% Save workspace variables

% save("data_for_animation.mat", "-regexp", "^(?!(M)$).");

%% Gráficas

figure(1)
subplot(3,1,1)
plot(t, X_e, 'LineWidth', 1.5)
title("Position States",'Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Position [m]','Interpreter','latex')
legend('$x_1$', '$x_2$','Interpreter','latex')
grid on
subplot(3,1,2)
plot(t, X_e_dot, 'LineWidth', 1.5)
title("Velocity States",'Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Velocity [m/s]','Interpreter','latex')
legend('$\dot{x}_1$', '$\dot{x}_2$','Interpreter','latex')
grid on
subplot(3,1,3)
stairs(t, u_plot, 'LineWidth', 1.5)
title("Control actions",'Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Force [N]','Interpreter','latex')
legend('$u_1$', '$u_2$','Interpreter','latex')
grid on

figure(2)
plot(t, Varepsilon_reg, "LineWidth",1.5)
title("Ergodic Metric",'Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('$\varepsilon \left( \mathbf{X_e}(t), \Phi(\mathbf{x}) \right) $','Interpreter','latex')
grid on

fig3h = figure(3);
layout3h = tiledlayout(2,2);
title(layout3h, "Ergodic Control", "interpreter", "latex", "FontSize", 20)

layout3h.TileSpacing = "compact";
layout3h.Padding = "compact";

nexttile
pcolor(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "EdgeColor", "none", "FaceColor", "interp")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Search Space")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
axis equal tight
grid on
legend('$\hat{\Phi}(\mathbf{x})$','Location','best')

nexttile
pcolor(x_1_grid, x_2_grid, reshape(C_x, length(x_2), length(x_1)),...
    "EdgeColor", "none", "FaceColor", "interp")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Empirical distribution reconstruction")
xlabel('$x_1$')
ylabel('$x_2$')
axis equal tight
grid on
hold on
plot(X_e(:,1), X_e(:,2), "k",'LineWidth',2)
plot(X_e(1,1), X_e(1,2),'ksq','MarkerSize',15,'LineWidth',2)
legend('$C(\mathbf{x})$','$\mathbf{X_e}(t)$', '$\mathbf{X_e}(0)$',...
    'Location','best')
hold off

nexttile
surf(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "EdgeColor", "none", "FaceColor", "interp")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Reference PDF",'Interpreter','latex')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
zlabel('$\Phi(\mathbf{x})$','Interpreter','latex')
grid on

nexttile
surf(x_1_grid, x_2_grid, reshape(C_x, length(x_2), length(x_1)),...
    "EdgeColor", "none", "FaceColor", "interp")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Empirical distribution reconstruction",'Interpreter','latex')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
zlabel('$C(\mathbf{x})$','Interpreter','latex')
grid on

set(findall(fig3h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig3h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig3h, "-property", "FontSize"), "FontSize", 20)