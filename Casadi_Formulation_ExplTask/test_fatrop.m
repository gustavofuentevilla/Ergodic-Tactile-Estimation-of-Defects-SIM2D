close all
clear
% clearvars -except Erg_traj_ipopt
clc

import casadi.*

%% Parámetros del espacio de búsqueda U = [L_1_l, L_1_u] \times [L_2_l, L_2_u]

n = 2; % Número de dimensiones espaciales

L_1_l = 1.0;
L_1_u = 1.5;
dx_1 = (L_1_u - L_1_l)/100;

L_2_l = 1.0;
L_2_u = 1.5;
dx_2 = (L_2_u - L_2_l)/100;

% Dimensiones \mathbf{x} = [x_1 x_2]^T
x_1 = (L_1_l:dx_1:L_1_u)';
x_2 = (L_2_l:dx_2:L_2_u)';

%vector de límites inferior y superiores de las dimensiones
L_i_l = [L_1_l, L_2_l];
L_i_u = [L_1_u, L_2_u];

[x_1_grid, x_2_grid] = meshgrid(x_1, x_2);

%Espacio de búsqueda discretizado
Omega = [reshape(x_1_grid,[],1), reshape(x_2_grid,[],1)]; 

%% Uniform PDF as First distribution

Phi_hat_x_1 = unifpdf(x_1, L_1_l, L_1_u);
Phi_hat_x_2 = unifpdf(x_2, L_2_l, L_2_u);

[Phi_hat_x_1_grid, Phi_hat_x_2_grid] = meshgrid(Phi_hat_x_1, Phi_hat_x_2);
Phi_hat_x = prod([reshape(Phi_hat_x_1_grid,[],1), reshape(Phi_hat_x_2_grid,[],1)], 2);

% Gaussian Mixture distributions as next PDFs

% Define parameters for Gaussian Mixture
n_def = 3; % Number of Gaussian components

for i = 41 : -10 : 1

    Mu = L_i_l + (L_i_u - L_i_l).*rand(n_def, 2); % Random means for each component
    
    % Matriz de Covarianza
    Var_def = 0.0005*i;
    Cov_1 = [Var_def, 0;
             0, Var_def];
    Cov_2 = [Var_def, -Var_def/2;
            -Var_def/2, Var_def];
    Cov_3 = [Var_def, Var_def/2;
            Var_def/2, Var_def];
    
    Sigma = cat(3,Cov_1,Cov_2,Cov_3);
    
    % Distribución
    gm_dist = gmdistribution(Mu, Sigma);
    
    % PDF de referencia REAL
    Phi = pdf(gm_dist, Omega);
    
    Phi_hat_x = cat(2, Phi_hat_x, Phi);

end

%% Visualize distributions

% figure(1)
% contour(x_1_grid, x_2_grid, reshape(Phi_hat_x(:,6), length(x_2), length(x_1)))
% xlim([L_1_l, L_1_u])
% ylim([L_2_l, L_2_u])
% xlabel('$x_1$ [m]','Interpreter','latex')
% ylabel('$x_2$ [m]','Interpreter','latex')
% axis equal
% grid on
% legend('$\hat{\Phi}(\mathbf{x})$','Interpreter','latex','Location','northeastoutside')


%% Cálculo de los coeficientes de Fourier para la PDF de referencia

% Coeficientes por dimensión
K = 12;

% Conjunto de valores para k_i
k_1 = (0:K-1)';
k_2 = (0:K-1)';

[k_1_grid, k_2_grid] = meshgrid(k_1, k_2);

% Conjunto de vectores índice
K_cal = [reshape(k_1_grid,1,[]); reshape(k_2_grid,1,[])];

Par_struct.K = K;
Par_struct.n = n;
Par_struct.K_cal = K_cal;
Par_struct.Omega = Omega;
Par_struct.dx_1 = dx_1;
Par_struct.dx_2 = dx_2;
Par_struct.L_i_l = L_i_l;
Par_struct.L_i_u = L_i_u;

[phi_k_reg, f_k_reg, h_k_reg] = FourierCoef_RefPDF(Phi_hat_x(:,1), Par_struct);

%% Optimal problem parameters

N = 200;             % Número de puntos de trayectoria óptima
t_f = 10;           %Tiempo final por iteración
T_s = t_f/N;                  % Tiempo de muestreo
t = (0:T_s:t_f)';   %Vector de tiempo por iteración

R = [3e-04, 0;
     0, 3e-04]*(1/T_s); %N = 200; 7e-5

% Peso sobre métrica ergódica
gamma = 1;

% Estado inicial z = [z_1; z_2; z_3; z_4] = [x_1; x_1_dot; x_2; x_2_dot]
z_0 = [1.0; 0; 1.0; 0]; 

%Pre-cálculo de Lambda
p = 2; %norma 2
Lambda_k = (1 + vecnorm(K_cal, p, 1)').^(-(n + 1)/2);

%% %%%%%%%%%%%%% Casadi Problem Setup FATROP %%%%%%%%%%%%%%%%%%

% Ecuaciones x_1_ddot = u_1;    x_2_ddot = u_2;
z = MX.sym('z', 4); %states z = [z(1); z(2); z(3); z(4)] = [x_1; x_1_dot; x_2; x_2_dot]
u = MX.sym('u', 2); %controls u = [u(1); u(2)]

%ODE construction: z_dot = [z(2); u(1); z(4); u(2)]
z_dot = [z(2); u(1); z(4); u(2)];
f = Function('f', {z,u}, {z_dot}, {'z', 'u'}, {'z_dot'});

%DAE problem structure

intg_options = struct;
intg_options.tf = T_s;   %Integration time (one step ahead)
intg_options.simplify = true;
intg_options.number_of_finite_elements = 4; %intermediate steps on the integration horizon

dae = struct;

dae.x = z;  % states (formalized)
dae.p = u;  % parameter, fixed during integration horizon (just one step ahead)
dae.ode = f(z,u); % symbolic dynamics

intg = integrator('intg', 'rk', dae, intg_options); %RK4 integration method

%One step integration (numerical)
% res = intg('x0', z_0, 'p', [0.5; 0]);  %z_0 initial condition, p = u controls
% z_next = res.xf;

%One step integration (symbolic)
res = intg('x0', z, 'p', u);
z_next = res.xf;

F = Function('F', {z,u}, {z_next}, {'z','u'}, {'z_next'});

%% Multiple Shooting for one prediction horizon with N+1 samples
opti = casadi.Opti();

% z = opti.variable(4, N+1);
% u = opti.variable(2, N);

nz = numel(z);
nu = numel(u);
nc = K^n;
nx = nz + nc;    % 148 estados dinámicos en total

% 1. VARIABLES COMO CELDAS INDIVIDUALES
X = cell(1, N+1);
U = cell(1, N);

for k = 1:N
    X{k} = opti.variable(nx, 1);
    U{k} = opti.variable(nu, 1);
end
X{N+1} = opti.variable(nx, 1);

% Parámetro de referencia de Fourier (phi_k)
phi_k_sym = opti.parameter(nc, 1);
opti.set_value(phi_k_sym, phi_k_reg);

% 2. FUNCIÓN DE TRANSICIÓN DE ESTADO COMPLETA (X_k, U_k) -> X_{k+1}
X_sym = MX.sym('X_sym', nx, 1);
U_sym = MX.sym('U_sym', nu, 1);

z_sym = X_sym(1:4);
c_sym = X_sym(5:end);

z_next_sym = F(z_sym, U_sym);

temp_1 = cos(K_cal(1,:)'.*pi.*(z_sym(1) - L_1_l)./(L_1_u - L_1_l));
temp_2 = cos(K_cal(2,:)'.*pi.*(z_sym(3) - L_2_l)./(L_2_u - L_2_l));
f_k_traj_sym = (temp_1 .* temp_2) ./ h_k_reg;
c_next_sym = c_sym + (f_k_traj_sym * T_s) / t_f;

X_next_sym = [z_next_sym; c_next_sym];
F_full = Function('F_full', {X_sym, U_sym}, {X_next_sym});

% 3. ASIGNACIÓN DIRECTA DE COTAS DE VARIABLES (lbx / ubx)
% Separar en '>=' y '<=' con vectores numéricos evita entrar al vector g(x)

X_init = [z_0; zeros(nc, 1)];

% Condición inicial para k = 1 (z_0 y c_1 = 0)
opti.subject_to( X{1} >= X_init );
opti.subject_to( X{1} <= X_init );

% Cotas para las demás etapas (k = 2 ... N+1)
inf_c = inf(nc, 1);
X_lb_base = [L_1_l; -inf; L_2_l; -inf; -inf_c];
X_ub_base = [L_1_u;  inf; L_2_u;  inf;  inf_c];

for k = 2:N+1
    lb_k = X_lb_base;
    ub_k = X_ub_base;
    
    % Velocidades fijas a cero en muestras clave
    if k == 2 || k == N || k == N+1
        lb_k(2) = 0; ub_k(2) = 0;
        lb_k(4) = 0; ub_k(4) = 0;
    end
    
    opti.subject_to( X{k} >= lb_k );
    opti.subject_to( X{k} <= ub_k );
end

% Cotas de control
for k = 1:N
    opti.subject_to( U{k} >= [-50; -50] );
    opti.subject_to( U{k} <= [50; 50] );
end

% 4. ÚNICAS IGUALDADES '==' DEL PROBLEMA: TRANSICIONES DINÁMICAS
% Esto garantiza que g[0] sea la primera ecuación de dinámica (gap-closing)
J = 0;
for k = 1:N
    % Transición de estado (Única igualdad en g(x))
    opti.subject_to( X{k+1} == F_full(X{k}, U{k}) );

    % Función de costo de la etapa
    c_next = X{k+1}(5:end);
    Varepsilon_sym = sum( Lambda_k .* (c_next - phi_k_sym).^2 );
    J = J + gamma * Varepsilon_sym + U{k}' * R * U{k} * T_s;
end

opti.minimize( J );

% 5. CONFIGURACIÓN Y RESOLUCIÓN CON FATROP
solver_opts = struct;
solver_opts.expand = true;
solver_opts.structure_detection = 'auto';
solver_opts.debug = true;

opti.solver('fatrop', solver_opts);

% Inicialización recomendada
for k = 1:N+1
    opti.set_initial(X{k}, X_init);
end

% Resolver
solution = opti.solve();

% 6. EXTRACCIÓN DE RESULTADOS
X_opt = zeros(nx, N+1);
for k = 1:N+1
    X_opt(:, k) = solution.value(X{k});
end

U_opt = zeros(nu, N);
for k = 1:N
    U_opt(:, k) = solution.value(U{k});
end

z_opt = X_opt(1:4, :);     % Trayectoria física (x1, v1, x2, v2)
c_opt = X_opt(5:end, :);   % Coeficientes ergódicos

% M_fatrop = opti.to_function('M_fatrop', {z_0_sym, phi_k_sym}, {z, u}, {'z_0','phi_k'}, {'z','u'});

