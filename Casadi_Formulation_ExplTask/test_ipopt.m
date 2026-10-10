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

%% %%%%%%%%%%%%% Casadi Problem Setup IPOPT %%%%%%%%%%%%%%%%%%

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

z = opti.variable(4, N+1);
u = opti.variable(2, N);
z_0_sym = opti.parameter(4, 1);  % parameter (not optimized over): initial condition
%u_0_sym = opti.parameter(2, 1);
phi_k_sym = opti.parameter(K^n, 1);
%u_d_sym = opti.parameter(2, 1);

% Symbolic Fourier functions, coefficients and ergodic metric with casadi
X_e_sym = [z(1,:)', z(3,:)'];     %Position [x_1, x_2] for all N samples

c_k_sym = 0;
f_k_traj_sym = opti.variable(K^n,1);
J = 0;
for i = 1:N
    for j = 1:K^n
        %problems using(.*) with casadi when K_cal is a matrix
        temp = cos(K_cal(:,j)'.*pi.*(X_e_sym(i,:) - L_i_l)./(L_i_u - L_i_l));   
        %problems using prod() function
        f_k_traj_sym(j,1) = temp(1)*temp(2)/h_k_reg(j);
    end
    c_k_sym = c_k_sym + (f_k_traj_sym*T_s)/(t_f);
    Varepsilon_sym = sum( Lambda_k.*(c_k_sym - phi_k_sym).^2 );

    % Objetive function
    J = J + gamma*Varepsilon_sym + u(:,i)'*R*u(:,i)*T_s; %(u(:,i) - u_d)   u(:,i)'*R*u(:,i)
end

% f_k_traj_sym = prod(cos( K_cal'.*pi.*(X_e_sym(i,:) - L_i_l)./(L_i_u - L_i_l) ), 2) ./ h_k_reg ;
% c_k_sym = c_k_sym + (f_k_traj_act*T_s)/t_f ;
% Varepsilon_act = sum( Lambda_k .* (c_k_act - phi_k_reg).^2 ); 

opti.minimize( J );

% Equality Constraints
for k = 1:N
    opti.subject_to( z(:,k+1) == F( z(:,k),u(:,k) ) );
end
opti.subject_to( z(1,1) == z_0_sym(1) ); % posiciones iniciales
opti.subject_to( z(3,1) == z_0_sym(3) );
opti.subject_to( z(2,1:2) == 0 ); % velocidades iniciales cero
opti.subject_to( z(4,1:2) == 0 );
opti.subject_to( z(2,end-1:end) == 0 ); % velocidades finales cero
opti.subject_to( z(4,end-1:end) == 0 );
% opti.subject_to( u(:,1) == u_0_sym ); % controles (aceleraciones) iniciales

% Inequality Constraints
opti.subject_to( -50 <= u <= 50 );
opti.subject_to( L_1_l <= z(1,:) <= L_1_u );    % x_1 boundaries
opti.subject_to( L_2_l <= z(3,:) <= L_2_u );    % x_2 boundaries

% constrain the derivate of control inputs
% for k = 1:N-1
%     opti.subject_to( abs(u(:,k+1) - u(:,k))/T_s <= 300 ); % 15 N
% end

% Solver definition
% p_opts = struct('expand',true);
% s_opts = struct('max_iter',100);
solver_opts = struct;
solver_opts.expand = true;
solver_opts.ipopt.max_iter = 2000;
opts.ipopt.print_level = 0;
%opts.ipopt.acceptable_tol = 1e-8;
%opts.ipopt.acceptable_obj_change_tol = 1e-6;

opti.solver('ipopt', solver_opts);

% opti.solver('sqpmethod', struct('qpsol','qrqp'));

%setting the initial conditions
opti.set_value(z_0_sym, z_0);
%opti.set_value(u_0_sym, u_0);
opti.set_value(phi_k_sym, phi_k_reg);
%opti.set_value(u_d_sym, u_d);

%Solving the optimization problem over the horizon N
solution = opti.solve();

% Function mapping from initial condition z0 to the optimal control action
%(first u of the N control actions). M contains IPOPT method embedded and
%integration method RK4

M_ipopt = opti.to_function('Erg_traj_ipopt',...
            {z_0_sym, phi_k_sym}, {z, u},...
            {'z_0','phi_k'}, {'z','u'});