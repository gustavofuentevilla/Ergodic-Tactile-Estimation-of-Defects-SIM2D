close all
clear
clc

import casadi.*

%% Search space \Omega = [L_1_l, L_1_u] \times [L_2_l, L_2_u]

n = 2;

L_1_l = 1.0;
L_1_u = 1.5;
dx_1 = (L_1_u - L_1_l)/100;

L_2_l = 1.0;
L_2_u = 1.5;
dx_2 = (L_2_u - L_2_l)/100;

% Axes \mathbf{x} = [x_1 x_2]^T
x_1 = (L_1_l:dx_1:L_1_u)';
x_2 = (L_2_l:dx_2:L_2_u)';

% Spatial Bounds (L_i_l = Lower bounds, L_i_u = Upper bounds)
L_i_l = [L_1_l, L_2_l];
L_i_u = [L_1_u, L_2_u];

[x_1_grid, x_2_grid] = meshgrid(x_1, x_2);

% Discretized Search space
Omega = [reshape(x_1_grid,[],1), reshape(x_2_grid,[],1)];

%% Uniform PDF as some testing distribution

Phi_hat_x_1 = unifpdf(x_1, L_1_l, L_1_u);
Phi_hat_x_2 = unifpdf(x_2, L_2_l, L_2_u);

[Phi_hat_x_1_grid, Phi_hat_x_2_grid] = meshgrid(Phi_hat_x_1, Phi_hat_x_2);
Phi_hat_x = prod([reshape(Phi_hat_x_1_grid,[],1),...
                  reshape(Phi_hat_x_2_grid,[],1)], 2);

%% Fourier coefficients computation

% Number of coefficients per dimension
K = 12;

% Build the index vectors set
k_1 = (0:K-1)';
k_2 = (0:K-1)';

[k_1_grid, k_2_grid] = meshgrid(k_1, k_2);

K_cal = [reshape(k_1_grid,1,[]); reshape(k_2_grid,1,[])];

% Parameters struct
Par_struct.K = K;
Par_struct.n = n;
Par_struct.K_cal = K_cal;
Par_struct.Omega = Omega;
Par_struct.dx_1 = dx_1;
Par_struct.dx_2 = dx_2;
Par_struct.L_i_l = L_i_l;
Par_struct.L_i_u = L_i_u;

% Coefficients for uniform distribution
[phi_k_reg, f_k_reg, h_k_reg] = FourierCoef_RefPDF(Phi_hat_x, Par_struct);

%% Optimal problem parameters

N = 200;             % Number of points in trajectory
t_f = 10;            % Time window
T_s = t_f/N;         % Sampling time
t = (0:T_s:t_f)';    % Time vector

% Control effort weights
R = [3e-04, 0;
    0, 3e-04]*(1/T_s);

% Ergodic metric weight
gamma = 1;

% Initial state z = [z_1; z_2; z_3; z_4] = [x_1; x_1_dot; x_2; x_2_dot]
z_0 = [1.0; 0; 1.0; 0]; 

% Lambda pre-computation
p = 2;
Lambda_k = (1 + vecnorm(K_cal, p, 1)').^(-(n + 1)/2);

%% %%%%%%%%%%%%% Casadi Problem Setup %%%%%%%%%%%%%%%%%%

% Dynamics x_1_ddot = u_1;    x_2_ddot = u_2;
% States z = [z(1); z(2); z(3); z(4)] = [x_1; x_1_dot; x_2; x_2_dot]
% Controls u = [u(1); u(2)]
z = MX.sym('z', 4); 
u = MX.sym('u', 2); 

% ODE construction: z_dot = [z(2); u(1); z(4); u(2)]
z_dot = [z(2); u(1); z(4); u(2)];
f = Function('f', {z,u}, {z_dot}, {'z', 'u'}, {'z_dot'});

% DAE problem structure
intg_options = struct;
intg_options.tf = T_s;   % Integration time (one step ahead)
intg_options.simplify = true;
intg_options.number_of_finite_elements = 4; % steps per integration horizon

dae = struct;

% formalization of dae problem
dae.x = z;          % states
dae.p = u;          % parameters, fixed during integration horizon
dae.ode = f(z,u);   % symbolic dynamics

intg = integrator('intg', 'rk', dae, intg_options); %RK4 integration method

% One step integration (numerical, for debbuging)
% res = intg('x0', z_0, 'p', [0.5; 0]);
% z_next = res.xf;

% One step integration (symbolic)
res = intg('x0', z, 'p', u);
z_next = res.xf;

% build the object with the dynamics
F = Function('F', {z,u}, {z_next}, {'z','u'}, {'z_next'});

%% Multiple Shooting for one prediction horizon with N+1 samples
opti = casadi.Opti();

% Variables (states and controls)
z = opti.variable(4, N+1);
u = opti.variable(2, N);
% Parameters (not optimized over: 
% initial conditions and the coefficients of the reference PDF)
z_0_sym = opti.parameter(4, 1);
phi_k_sym = opti.parameter(K^n, 1);

% Symbolic Fourier functions, coefficients and ergodic metric
X_e_sym = [z(1,:)', z(3,:)'];     % Positions [x_1, x_2] for all N samples

c_k_sym = 0;
f_k_traj_sym = opti.variable(K^n,1);
J = 0;
for i = 1:N
    for j = 1:K^n
        temp = cos(K_cal(:,j)'.*pi.*(X_e_sym(i,:) - L_i_l)./...
               (L_i_u - L_i_l));
        f_k_traj_sym(j,1) = temp(1)*temp(2)/h_k_reg(j);
    end
    c_k_sym = c_k_sym + (f_k_traj_sym*T_s)/(t_f);
    Varepsilon_sym = sum( Lambda_k.*(c_k_sym - phi_k_sym).^2 );

    % Objetive function
    J = J + gamma*Varepsilon_sym + u(:,i)'*R*u(:,i)*T_s;
end

opti.minimize( J );

% Equality Constraints
for k = 1:N
    % Dynamics
    opti.subject_to( z(:,k+1) == F( z(:,k),u(:,k) ) );
end
% Initial positions
opti.subject_to( z(1,1) == z_0_sym(1) ); 
opti.subject_to( z(3,1) == z_0_sym(3) );
% Initial velocities
opti.subject_to( z(2,1:2) == 0 ); 
opti.subject_to( z(4,1:2) == 0 );
% Final velocities
opti.subject_to( z(2,end-1:end) == 0 ); 
opti.subject_to( z(4,end-1:end) == 0 );

% Inequality Constraints
opti.subject_to( -50 <= u <= 50 );
opti.subject_to( L_1_l <= z(1,:) <= L_1_u );    % x_1 boundaries
opti.subject_to( L_2_l <= z(3,:) <= L_2_u );    % x_2 boundaries

% Solver definition
solver_opts = struct;
solver_opts.expand = true;
solver_opts.ipopt.max_iter = 2000;
opts.ipopt.print_level = 0;

opti.solver('ipopt', solver_opts);

% Setting the initial conditions
opti.set_value(z_0_sym, z_0);
opti.set_value(phi_k_sym, phi_k_reg);

% Solving the optimization problem over the horizon N
solution = opti.solve();

% Function mapping
Erg_traj_ipopt = opti.to_function('Erg_traj_ipopt',...
                                   {z_0_sym, phi_k_sym}, {z, u},...
                                   {'z_0','phi_k'}, {'z','u'});

% Save Casadi object as a file
Erg_traj_ipopt.save('Erg_traj_ipopt.casadi');
