% Primera ejecución
% import casadi.*
% M = Function.load('/home/gustavo-fuentevilla/MATLAB/Tactile_Defects_Localization/Casadi_Formulation_ExplTask/M_N200.casadi');

% for run_idx = 1:100

% Resto de ejecuciones
close all
clearvars -except M %run_idx
clc

%% Parámetros del espacio de búsqueda U = [L_1_l, L_1_u] \times [L_2_l, L_2_u]

%idx_tmp = repelem( (1:10)' ,repelem(10, 10));

n = 2; % Número de dimensiones espaciales

L_1_l = 0.5;
dx_1 = 0.01;
L_1_u = 1.5;

L_2_l = 0.5;
dx_2 = 0.01;
L_2_u = 1.5;

% Dimensiones \mathbf{x} = [x_1 x_2]^T
x_1 = (L_1_l:dx_1:L_1_u)';
x_2 = (L_2_l:dx_2:L_2_u)';

%vector de límites inferior y superiores de las dimensiones
L_i_l = [L_1_l, L_2_l];
L_i_u = [L_1_u, L_2_u];

L_1 = (L_1_u - L_1_l);
L_2 = (L_2_u - L_2_l);

[x_1_grid, x_2_grid] = meshgrid(x_1, x_2);

%Espacio de búsqueda discretizado
Omega = [reshape(x_1_grid,[],1), reshape(x_2_grid,[],1)]; 

%% Defects definition

n_def = 10; %idx_tmp(run_idx); % Greater than 0

% [Mu, Sigma, r_elips_Phi] = DefectsGen(n_def, L_i_l, L_i_u);

% run_idx = 91;
% % Para replicar los resultados con defectos ya existentes en /results/
% load(sprintf("Results/output_%d.mat",run_idx), "Mu", "Sigma", "r_elips_Phi")

gm_dist = gmdistribution(Mu, Sigma);% , proporciones);

%PDF de referencia REAL
Phi_x = pdf(gm_dist, Omega);

%% Uniform PDF as an initial guess

Phi_hat_x_1 = unifpdf(x_1, L_1_l, L_1_u);
Phi_hat_x_2 = unifpdf(x_2, L_2_l, L_2_u);

[Phi_hat_x_1_grid, Phi_hat_x_2_grid] = meshgrid(Phi_hat_x_1, Phi_hat_x_2);
Phi_hat_x = prod([reshape(Phi_hat_x_1_grid,[],1), reshape(Phi_hat_x_2_grid,[],1)], 2);

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

[phi_k_reg, f_k_reg, h_k_reg] = FourierCoef_RefPDF(Phi_hat_x, Par_struct);

%% Condiciones Iniciales y parámetros

N = 200; % Número de muestras por iteración
t_f = 10;           %Tiempo final por iteración
T_s = t_f/N;                  % Tiempo de muestreo
t = (0:T_s:t_f)';   %Vector de tiempo por iteración

% Peso sobre controles
R = [7e-5, 0;
     0, 7e-5]*(1/T_s); % N = 200

% Peso sobre métrica ergódica
gamma = 1;

% Estado inicial z = [z_1; z_2; z_3; z_4] = [x_1; x_1_dot; x_2; x_2_dot]
z_0 = [0.5; 0; 0.5; 0]; 

%Pre-cálculo de Lambda
p = 2; %norma 2
Lambda_k = (1 + vecnorm(K_cal, p, 1)').^(-(n + 1)/2);


%% %%%%%%%%%%%%% Casadi Problem Setup %%%%%%%%%%%%%%%%%%

% import casadi.*
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

%% vector to add more points on the trajectory and get more data from sensor

t_spline = (0:0.01:t_f)'; %Time vector por spline in one iteration

%% Loop for the Search task
n_iter_max = 10;

% Registers
z_reg = zeros(N+1, 4, n_iter_max);
u_reg = zeros(N, 2, n_iter_max);
X_e_reg = zeros(N+1, 2, n_iter_max);
X_e_dot_reg = zeros(N+1, 2, n_iter_max);
phi_k_REG = zeros(K^n, 1, n_iter_max);
Phi_hat_x_reg = zeros(height(Omega), 1, n_iter_max + 1);
X_e_spline_reg = zeros(length(t_spline), 2, n_iter_max);
X_e_dot_spline_reg = zeros(length(t_spline), 2, n_iter_max);
u_spline_reg = zeros(length(t_spline) - 1, 2, n_iter_max);
V_Xe_reg = zeros(length(t_spline), 1, n_iter_max);

% Initializations
z_act = z_0;
%u_act = u_0;
phi_k_act = phi_k_reg;
Phi_hat_x_act = Phi_hat_x;

% Measurement parameters
a = 5; %Reference force
b = 0.2;
c = 0.01; %Amplitud del ruido en la medición
n_points = length(t_spline);

% Parameters for PDF Estimator
Par_PDF.Omega = Omega;
Par_PDF.dx_1 = dx_1;
Par_PDF.dx_2 = dx_2;
Par_PDF.Meas_mean = a;
% Range of possible Number of defects to be found
Par_PDF.nbDef_range = [1, n_def + 2]; 

Par_PDF.Prev_Data = [];
% Par_PDF.Prev_Priors = [];
% Par_PDF.Prev_Mu = [0, 0; 0, 0; 0, 0];
% Par_PDF.Prev_Sigma = repmat(diag([realmax, realmax]), 1, 1, n_def);
% Par_PDF.Prev_Sigma_a = repmat(diag([realmax, realmax]), 1, 1, n_def);
Par_PDF.Prev_numComponents = [];

% Define the dimensions of the registers for the defects found with an
% initial value (these has to be removed at the end)
Par_PDF.Prev_Mu_found = [0, 0];
Par_PDF.Prev_Sigma_found = [0, 0; 0, 0];
NoDataIterCounter = 0;

n_iter = n_iter_max;

Par_PDF.DataEscFact = 1;
% total variation condition to find a defect
% Par_PDF.Thres_Variation = 0.135; % 13.5 cm 
Par_PDF.Thres_Variation = max(sum(r_elips_Phi)) + 0.001;
% Minimum axes lengths of gaussian elipses
Par_PDF.MinAxisLengths = 0; % 0 cm 
% Limit for One cluster
Par_PDF.OneClustDistLimit = 2*max(r_elips_Phi,[],"all") + 0.07;
Par_PDF.flag_ExplorationStage = true;

% Parameters definition for Variation constraint function

% Porcentage of max variation constraint, 
% porcentage of MaxVarCons to match with first D_KL value
nu_p = 0.35;
% Another way to compute the MaxVarCons is define the number of times of
% Variation Threshold we want to cover per defect, \eta times.
eta = 6;
% D_KL that matches MaxVarCons
Par_PDF.D_KL_bar_u = []; % computed in estimator 
% Little offset under Variation Threshold for numerical estability
Par_PDF.eps = 0.001;

% Maximum Variation Constraint Computation (2 ways)
Par_PDF.MaxVarCons = nu_p*(L_1 + L_2) + (1 - nu_p)*...
                     (Par_PDF.Thres_Variation - Par_PDF.eps);
% Par_PDF.MaxVarCons = eta*Par_PDF.Thres_Variation;

for i = 1:n_iter_max

    [Z, U] = M(z_act, phi_k_act); %Soluciones
    %[Z, U] = M(z_act, u_act, phi_k_act); 
    Z = full(Z)';
    U = full(U)';

    % Trajectory X_e(t) = [x_1, x_2]
    X_e = [Z(:, 1), Z(:, 3)];
    X_e_dot = [Z(:, 2), Z(:, 4)];

    % Spline: Adding points to the trajectory to get more data from the sensor
    % and pass it to the estimator
    x_1e_spline = spline(t, X_e(:,1), t_spline);
    x_2e_spline = spline(t, X_e(:,2), t_spline);
    X_e_spline = [x_1e_spline, x_2e_spline];

    x_1e_dot_spline = spline(t, X_e_dot(:,1), t_spline);
    x_2e_dot_spline = spline(t, X_e_dot(:,2), t_spline);
    X_e_dot_spline = [x_1e_dot_spline, x_2e_dot_spline];

    u_1_spline = spline(t(1:end-1), U(:,1), t_spline(1:end-1));
    u_2_spline = spline(t(1:end-1), U(:,2), t_spline(1:end-1));
    u_spline = [u_1_spline, u_2_spline];

    % Measurement along the trajectory, V_Xe
    Upsilon = a + b*pdf(gm_dist, X_e_spline); %Real PDF
    delta = c*randn(n_points, 1); %Gaussian Noise with Variance c^2
    V_Xe = Upsilon + delta;
    
    Par_PDF.thres_meas = a + max(delta);

    % Registers
    z_reg(:,:,i) = Z;
    u_reg(:,:,i) = U;
    X_e_reg(:,:,i) = X_e;
    X_e_dot_reg(:,:,i) = X_e_dot;
    phi_k_REG(:,:,i) = phi_k_reg;
    Phi_hat_x_reg(:,:,i) = Phi_hat_x_act;
    X_e_spline_reg(:,:,i) = X_e_spline;
    X_e_dot_spline_reg(:,:,i) = X_e_dot_spline;
    u_spline_reg(:,:,i) = u_spline;
    V_Xe_reg(:,:,i) = V_Xe;

    % PDF Estimation
    Par_PDF.iteration = i;
    Par_PDF.Prev_Phi_hat_x = Phi_hat_x_act;
    [Phi_hat_x_next, Estim_sol(i)] = PDF_Estimator(X_e_spline, V_Xe, Par_PDF);

    % Update Iterations Counter where No data hav been found
    NoDataIterCounter = NoDataIterCounter + Estim_sol(i).flag_NoData;

    if NoDataIterCounter == 2
        n_iter = i;
        break;
    end

    % Save D_KL from first iteration to set the exploration function
    if i == 1
        Par_PDF.D_KL_bar_u = Estim_sol(i).D_KL;
    end

    % Detect the falling edge of Exploration flag termination
    if i > 1
        falling_edge = (Estim_sol(i).flag_ExplorationStage - ...
                        Estim_sol(i-1).flag_ExplorationStage) == -1;
    else
        falling_edge = false;
    end
    % Save the Exploration flag (turned off) and the number of components 
    if falling_edge
        Par_PDF.flag_ExplorationStage = Estim_sol(i).flag_ExplorationStage;
        Par_PDF.Prev_numComponents = Estim_sol(i).numComponents;
    end

    % Save found defects if any
    Par_PDF.Prev_Mu_found = cat(1, Par_PDF.Prev_Mu_found, Estim_sol(i).Mu_found);
    Par_PDF.Prev_Sigma_found = cat(3, Par_PDF.Prev_Sigma_found, Estim_sol(i).Sigma_found);

    if Estim_sol(i).flag_done
        n_iter = i; % Save number of iterations achieved
        break;
    end
    
    % Saving Data to use it as "Previous data" in next iterations
    Par_PDF.Prev_Data = Estim_sol(i).Data;
    % Par_PDF.Prev_Priors = Estim_sol(i).Priors;
    % Par_PDF.Prev_Mu = Estim_sol(i).Mu;
    % Par_PDF.Prev_Sigma = Estim_sol(i).Sigma;
    % Par_PDF.Prev_Sigma_a = Estim_sol(i).Sigma_a;

    % Compute new Fourier coefficients for \hat{Phi}(x)
    [phi_k_reg, ~, ~] = FourierCoef_RefPDF(Phi_hat_x_next, Par_struct);

    % Update parameters for next iteration
    z_act = Z(end,:)';           % Initial condition for state
    %u_act = U(end,:)';
    phi_k_act = phi_k_reg;      % New target coefficients
    Phi_hat_x_act = Phi_hat_x_next;

end

% Remove the initial value (zero values) for defects found
Mu_found = Par_PDF.Prev_Mu_found(2:end, :);
Sigma_found = Par_PDF.Prev_Sigma_found(:,:,2:end);


%% Reconstrucción de distribución empírica y métrica ergódica

Varepsilon_reg = zeros(length(t), 1, n_iter);
C_x_reg = zeros(height(Omega), length(t), n_iter);

for r = 1:n_iter
    c_k = zeros(K^n, 1);
    C_x = zeros(height(Omega), 1);
    for i = 1:length(t)

        % Compute Fourier Functions and coefficients on the new position
        f_k_traj = prod(cos( K_cal'.*pi.*(X_e_reg(i,:,r) - L_i_l)./(L_i_u - L_i_l) ), 2) ./ h_k_reg ;
        c_k = c_k + (f_k_traj*T_s)/(t_f) ; %i*T_s %t_f

        %Ergodic metric
        Varepsilon = sum( Lambda_k .* (c_k - phi_k_REG(:,:,r)).^2 );

        %Empirical distribution reconstruction
        C_x_i = zeros(height(Omega), 1);
        for j = 1:K^n
            C_x_i = C_x_i + c_k(j)*f_k_reg(:,j);
        end

        %Se suman todas las distribuciones generadas en cada muestra
        C_x = C_x + C_x_i;

        %Se registra
        C_x_reg(:,i,r) = C_x;
        Varepsilon_reg(i,:,r) = Varepsilon;

    end

end


%% Pre-Procesing for charts

t_total = zeros(n_iter*(length(t)-1) + 1, 1);
X_e_total = zeros(n_iter*(length(t)-1) + 1, 2);
X_e_dot_total = zeros(n_iter*(length(t)-1) + 1, 2);
Varepsilon_total = zeros(n_iter*(length(t)-1) + 1, 1);
u_total = zeros(n_iter*(length(t)-1) + 1, 2);

t_spline_total = zeros(n_iter*(length(t_spline)-1) + 1, 1);
X_e_spline_total = zeros(n_iter*(length(t_spline)-1) + 1, 2);
X_e_dot_spline_total = zeros(n_iter*(length(t_spline)-1) + 1, 2);
u_spline_total = zeros(n_iter*(length(t_spline)-1) + 1, 2);
V_Xe_total = zeros(n_iter*(length(t_spline)-1) + 1, 1);

for i = 1:n_iter

    id_init = ((i - 1)*(length(t)-1) + 1); %1,101,..
    id_last = (i*(length(t)-1) + 1);        %101, 201,...

    t_total( id_init:id_last ) = (i - 1)*t(end) + t;
    X_e_total( id_init:id_last, : ) = X_e_reg(:,:,i);
    X_e_dot_total( id_init:id_last, : ) = X_e_dot_reg(:,:,i);
    Varepsilon_total( id_init:id_last, : ) = Varepsilon_reg(:,:,i);
    u_total( id_init:id_last, : ) = [u_reg(:,:,i); [NaN, NaN]];
    
    id_init_spline = ((i - 1)*(length(t_spline)-1) + 1);
    id_last_spline = (i*(length(t_spline)-1) + 1);
    t_spline_total( id_init_spline:id_last_spline ) = (i - 1)*t_spline(end) + t_spline;
    X_e_spline_total( id_init_spline:id_last_spline, : ) = X_e_spline_reg(:,:,i);
    X_e_dot_spline_total( id_init_spline:id_last_spline, : ) = X_e_dot_spline_reg(:,:,i);
    u_spline_total( id_init_spline:id_last_spline, : ) = [u_spline_reg(:,:,i); [NaN, NaN]];
    V_Xe_total( id_init_spline:id_last_spline, : ) = V_Xe_reg(:,:,i);

end

%% Saving variables (Except M casadi function)

% Normal for-loop

% save(sprintf("Results/output_%d.mat",run_idx), "-regexp", "^(?!(M)$).");


% end


