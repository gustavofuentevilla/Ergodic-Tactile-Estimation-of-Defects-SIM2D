% % Ejecute primero estas dos lineas en consola
% import casadi.*
% M = Function.load('/home/gustavo-fuentevilla/MATLAB/Tactile_Defects_Localization/Casadi_Formulation_ExplTask/M_N200.casadi');

% for run_idx = 1:15

% Resto de ejecuciones
close all
clearvars -except M %run_idx
clc

%% Parámetros del espacio de búsqueda U = [L_1_l, L_1_u] \times [L_2_l, L_2_u]

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

[x_1_grid, x_2_grid] = meshgrid(x_1, x_2);

%Espacio de búsqueda discretizado
Omega = [reshape(x_1_grid,[],1), reshape(x_2_grid,[],1)]; 

%% Gaussian Mixture distribution (Real PDF)

% Media del Gaussiano
%load last run, enter on console  load("s.mat"); rng(s);
% s = rng;
% save("s.mat", "s")
offset = 0.05; %For random means
Mu_1 = (L_i_l + offset) + ((L_i_u - offset) - (L_i_l + offset)).*rand(1,2);
Mu_2 = (L_i_l + offset) + ((L_i_u - offset) - (L_i_l + offset)).*rand(1,2);
Mu_3 = (L_i_l + offset) + ((L_i_u - offset) - (L_i_l + offset)).*rand(1,2); 

Mu = [Mu_1; Mu_2; Mu_3];

%Testing 
% Mu = [1.4295    1.1986;	    0.9695    1.3885;	    0.7255    0.6193];

% Matriz de Covarianza
Var_def = 0.0005;
% Cov_1 = [Var_def, 0;
%          0, Var_def];
% Cov_2 = Cov_1;
% Cov_3 = Cov_1;

Cov_1 = [Var_def, Var_def/2;
         Var_def/2, Var_def];
Cov_2 = [Var_def, -Var_def/2;
         -Var_def/2, Var_def];
Cov_3 = Cov_1;

Sigma = cat(3,Cov_1,Cov_2,Cov_3);

% run_idx = 1;
% Para replicar los resultados con defectos ya existentes en /results/
% load(sprintf("results/output_%d.mat",run_idx), "Mu", "Sigma")

%Peso sobre el Gaussiano
n_def = size(Mu, 1);
%proporciones = [1/n_def, 1/n_def, 1/n_def];

gm_dist = gmdistribution(Mu, Sigma);% , proporciones);

%PDF de referencia REAL
Phi_x = pdf(gm_dist, Omega);

%% Geometric computations of real PDF

stdev_Phi = zeros(size(Sigma));
Sigma_ast_Phi = zeros(size(Sigma));
r_elips_Phi = zeros(n, n_def);
variation_Phi = zeros(1, n_def);
for j = 1:n_def
    stdev_Phi(:,:,j) = sqrtm(Sigma(:,:,j)); % Standard deviation
    Sigma_ast_Phi(:,:,j) = 3*stdev_Phi(:,:,j);
    % radios de ejes principales 
    r_elips_Phi(:,j) = eig(Sigma_ast_Phi(:,:,j));
    % Variación total = trace(Sigma_ast_phi) = sum(r_elips_phi)
    variation_Phi(:,j) = sum(r_elips_Phi(:,j));
end

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
     0, 7e-5]*(1/T_s);

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

n_iter = 5;

% Registers
z_reg = zeros(N+1, 4, n_iter);
u_reg = zeros(N, 2, n_iter);
X_e_reg = zeros(N+1, 2, n_iter);
X_e_dot_reg = zeros(N+1, 2, n_iter);
phi_k_REG = zeros(K^n, 1, n_iter);
Phi_hat_x_reg = zeros(height(Omega), 1, n_iter + 1);
Phi_hat_x_1_reg = zeros(length(x_1), n_def, n_iter + 1);
Phi_hat_x_2_reg = zeros(length(x_2), n_def, n_iter + 1);
X_e_spline_reg = zeros(length(t_spline), 2, n_iter);
X_e_dot_spline_reg = zeros(length(t_spline), 2, n_iter);
u_spline_reg = zeros(length(t_spline) - 1, 2, n_iter);
V_Xe_reg = zeros(length(t_spline), 1, n_iter);
V_hat_i_reg = zeros(length(t_spline), n_def, n_iter);

% Initializations
z_act = z_0;
%u_act = u_0;
phi_k_act = phi_k_reg;
Phi_hat_x_1_act = [Phi_hat_x_1, Phi_hat_x_1, Phi_hat_x_1];
Phi_hat_x_2_act = [Phi_hat_x_2, Phi_hat_x_2, Phi_hat_x_2];

% Measurement parameters
a = 5;
b = 0.2;
c = 0.01; %0.05 %Amplitud del ruido en la medición: Importante porque la estimación
% del siguiente PDF depende fuertemente de la amplitud de medición, si el
% ruido es muy grande el estimador puede fallar.
n_points = length(t_spline);

% Parameters for PDF Estimator
Par_PDF.x_1 = x_1;
Par_PDF.x_2 = x_2;
Par_PDF.Meas_mean = a;
Par_PDF.n_def = n_def;

%Parameters for Independent Measurements modeling
Par_MeasMdl.a = a;
Par_MeasMdl.beta = 0; % 0.01;
Par_MeasMdl.n_def = n_def;

Phi_hat_x_1_next = zeros(length(x_1), n_def);
Phi_hat_x_2_next = zeros(length(x_2), n_def);
Phi_hat_x_j = zeros(length(x_1)*length(x_2), n_def);

% V_bar = zeros(length(t_spline), n_def);
V_bar_i_reg = zeros(length(t_spline), n_def, n_iter);

PrevUse_Data = [];
% Par_MeasMdl.Prev_Priors = [];
% Par_MeasMdl.Prev_Mu = []; %[0, 0; 0, 0; 0, 0];
% Par_MeasMdl.Prev_Sigma = []; % repmat(diag([realmax, realmax]), 1, 1, n_def);

for i = 1:n_iter

    [Z, U] = M(z_act, phi_k_act);
    %[Z, U] = M(z_act, u_act, phi_k_act); %Soluciones
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
    
    Par_MeasMdl.thres_meas = a + max(delta);
    Par_MeasMdl.iter = i;

    [V_hat_i, MdlSol_struct(i)] = Indep_MeasMdl(V_Xe, X_e_spline, Par_MeasMdl, PrevUse_Data);
    
    %Saving sorted GMM for next iteration
    Par_MeasMdl.GMM_ordered_last = MdlSol_struct(i).GMM_ordered;
    
    %Update of preprocessed data from iteration 1 to actual to use in the
    %next iteration
    PrevUse_Data = MdlSol_struct(i).Preprocessed_Data;

    % Save gaussian parameters from last solution (feedback)
    % Par_MeasMdl.Prev_Priors = MdlSol_struct(i).GMModel.ComponentProportion;
    % Par_MeasMdl.Prev_Mu = MdlSol_struct(i).GMModel.mu;
    % Par_MeasMdl.Prev_Sigma = MdlSol_struct(i).GMModel.Sigma;

    % Sequence of V above the threshold
    V_bar_i = MdlSol_struct(i).V_bar_i;

    % Registers
    z_reg(:,:,i) = Z;
    u_reg(:,:,i) = U;
    X_e_reg(:,:,i) = X_e;
    X_e_dot_reg(:,:,i) = X_e_dot;
    phi_k_REG(:,:,i) = phi_k_reg;
    Phi_hat_x_reg(:,:,i) = Phi_hat_x;
    Phi_hat_x_1_reg(:,:,i) = Phi_hat_x_1_act;
    Phi_hat_x_2_reg(:,:,i) = Phi_hat_x_2_act;
    X_e_spline_reg(:,:,i) = X_e_spline;
    X_e_dot_spline_reg(:,:,i) = X_e_dot_spline;
    u_spline_reg(:,:,i) = u_spline;
    V_Xe_reg(:,:,i) = V_Xe;
    V_hat_i_reg(:,:,i) = V_hat_i;
    V_bar_i_reg(:,:,i) = V_bar_i;

    minlen_low = 0.085; % 8.5 cm
    minlen_up = 0.6;    % 70 cm

    Par_PDF.MinLength_x_j = -tanh((i-2)*5)*(minlen_up - minlen_low)/2 +...
                            (minlen_up + minlen_low)/2;

    % PDF Estimation
    for j = 1:n_def
        [Phi_hat_x_1_next(:,j), Phi_hat_x_2_next(:,j), Estim_sol(i,j)] = PDF_Estimator( ...
            Phi_hat_x_1_act(:,j), Phi_hat_x_2_act(:,j), X_e_spline, V_hat_i(:,j), Par_PDF );

        [Phi_hat_x_1_next_grid, Phi_hat_x_2_next_grid] = meshgrid( ...
            Phi_hat_x_1_next(:,j), Phi_hat_x_2_next(:,j) );

        Phi_hat_x_j(:,j) = prod([reshape(Phi_hat_x_1_next_grid,[],1), reshape(Phi_hat_x_2_next_grid,[],1)], 2);
    end

    Phi_hat_x = (1/n_def)*sum(Phi_hat_x_j, 2);

    % Compute new Fourier coefficients for \hat{Phi}(x)
    [phi_k_reg, ~, ~] = FourierCoef_RefPDF(Phi_hat_x, Par_struct);

    % Update parameters for next iteration
    z_act = Z(end,:)';           % Initial condition for state
    %u_act = U(end,:)';
    phi_k_act = phi_k_reg;      % New target coefficients
    Phi_hat_x_1_act = Phi_hat_x_1_next;
    Phi_hat_x_2_act = Phi_hat_x_2_next;   

end

% Save last estimation in registers
Phi_hat_x_reg(:,:,n_iter + 1) = Phi_hat_x;
Phi_hat_x_1_reg(:,:,n_iter + 1) = Phi_hat_x_1_act;
Phi_hat_x_2_reg(:,:,n_iter + 1) = Phi_hat_x_2_act;

% After the last iteration, the next PDF is still computed, so we have one
% more mu_hat and sigma_hat, even if it's not used for ergodic control

for i = 1:n_def
    Mu_hat(i,:) = [Estim_sol(n_iter,i).exp_x1_hat, Estim_sol(n_iter,i).exp_x2_hat];
    Sigma_hat(:,:,i) = diag([Estim_sol(n_iter,i).var_x1_hat, Estim_sol(n_iter,i).var_x2_hat]);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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


%% Pre-processing for ploting

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
% save(sprintf("results/output_%d.mat",run_idx), "-regexp", "^(?!(M)$).");
% 
% end

%% Charts

clearvars -except M
close all

load("results/output_13.mat")  %fails 5, 13

fig1h = figure(1);
layout1h = tiledlayout(fig1h, 3, 1);
nexttile
plot(t_total, X_e_total, 'LineWidth', 2)
title("Position States")
xlabel('Time [s]')
ylabel('Position [m]')
legend('$x_1$', '$x_2$')
grid on
nexttile
plot(t_total, X_e_dot_total, 'LineWidth', 2)
title("Velocity States")
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('$\dot{x}_1$', '$\dot{x}_2$')
grid on
nexttile
plot(t_total, u_total, 'LineWidth', 2) %stairs
title("Control actions")
xlabel('Time [s]')
ylabel('Force [N]')
legend('$u_1$', '$u_2$')
grid on

set(findall(fig1h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig1h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig1h, "-property", "FontSize"), "FontSize", 18)

fig2h = figure(2);
plot(t_total, Varepsilon_total, "LineWidth",2)
title("Ergodic Metric")
xlabel('Time [s]')
ylabel('$\varepsilon \left( \mathbf{X_e}(t), \Phi(\mathbf{x}) \right) $')
grid on

set(findall(fig2h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig2h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig2h, "-property", "FontSize"), "FontSize", 18)

%% 2D charts

if n_iter <= 2
    aux = 1;
elseif n_iter > 2 && n_iter <= 5
    aux = 2;
elseif n_iter > 5 && n_iter <= 8
    aux = 3;
elseif n_iter > 8 && n_iter <= 11
    aux = 4;
elseif n_iter > 11 && n_iter <= 14
    aux = 5;
end

nbDrawingSeg = 100;
tmp_vec = linspace(-pi, pi, nbDrawingSeg)';
Elipse_Phi = zeros(height(tmp_vec), 2, n_def); %Elipse
for j = 1:n_def
    Elipse_Phi(:,:,j) = [cos(tmp_vec), sin(tmp_vec)] * real(Sigma_ast_Phi(:,:,j)) + repmat(Mu(j,:),nbDrawingSeg,1);
end

contornos = 25;

fig3h = figure(3);
subplot(aux,3,1)
pcolor(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "EdgeColor","none","FaceColor","interp")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Real PDF")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
axis equal tight
grid on
colormap("default")
hold on
for j = 1:n_def
    plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "-.k", "LineWidth", 1.3)
end
plot(Mu(:,1), Mu(:,2), ".", 'MarkerSize', 8)
hold off
legend('$\Phi(\mathbf{x})$','Location', 'best') %'northeastoutside')
for i = 1:n_iter
    subplot(aux,3,i+1)
    contour(x_1_grid, x_2_grid, reshape(Phi_hat_x_reg(:,:,i), length(x_2), length(x_1)) )%, contornos)
    %pcolor(x_1_grid, x_2_grid, reshape(Phi_hat_x_reg(:,:,i), length(x_2), length(x_1)),...
        % "FaceColor","interp","EdgeColor","none")
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    title("Estimated PDF, iteration " + i)
    xlabel('$x_1$ [m]')
    ylabel('$x_2$ [m]')
    axis equal
    grid on
    hold on
    plot(X_e_reg(:,1,i), X_e_reg(:,2,i),'LineWidth',1.5)
    plot(X_e_reg(1,1,i), X_e_reg(1,2,i),'ksq','MarkerSize',7,'LineWidth',1.5)
    for j = 1:n_def
        plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "-.k", "LineWidth", 1)
    end
    plot(Mu(:,1),Mu(:,2),'.k','MarkerSize',8)
    legend('$\hat{\Phi}(\mathbf{x})$','$\mathbf{X_e}(t)$', '$\mathbf{X_e}(0)$',...
        "Real" + newline + "defects",'Location','best')%'northeastoutside')
    hold off
end

set(findall(fig3h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig3h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig3h, "-property", "FontSize"), "FontSize", 16)

%% Gráficas con adición de puntos Spline

figure(4)
subplot(3,1,1)
plot(t_spline_total, X_e_spline_total, 'LineWidth', 1.5)
title("Position States",'Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Position [m]','Interpreter','latex')
legend('$x_1$', '$x_2$','Interpreter','latex')
grid on
subplot(3,1,2)
plot(t_spline_total, X_e_dot_spline_total, 'LineWidth', 1.5)
title("Velocity States",'Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Velocity [m/s]','Interpreter','latex')
legend('$\dot{x}_1$', '$\dot{x}_2$','Interpreter','latex')
grid on
subplot(3,1,3)
plot(t_spline_total, u_spline_total, 'LineWidth', 1.5) %stairs
title("Control actions",'Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Force [N]','Interpreter','latex')
legend('$u_1$', '$u_2$','Interpreter','latex')
grid on

figure(5)
subplot(3,1,1)
plot(t_spline_total, V_Xe_total, 'LineWidth', 1.5)
title("Sensor measurement on time",'Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Force [N]','Interpreter','latex')
legend('$V(t)$','Interpreter','latex')
grid on

for i = 1:n_iter

    subplot(3,1,2)
    hold on
    scatter(X_e_spline_reg(:,1,i), V_Xe_reg(:,:,i), 14, 'filled',...
        'DisplayName', "$V(x_1)$, iteration " + i, "MarkerFaceAlpha", 0.2)
    title("Sensor measurement on spatial domain",'Interpreter','latex')
    xlabel('$x_1$ [m]','Interpreter','latex')
    ylabel('Force [N]','Interpreter','latex')
    grid on
    hold off
    legend('Interpreter','latex')
    subplot(3,1,3)
    hold on
    scatter(X_e_spline_reg(:,2,i), V_Xe_reg(:,:,i), 14, 'filled',...
        'DisplayName', "$V(x_2)$, iteration " + i, "MarkerFaceAlpha", 0.4)
    title("Sensor measurement on spatial domain",'Interpreter','latex')
    xlabel('$x_2$ [m]','Interpreter','latex')
    ylabel('Force [N]','Interpreter','latex')
    grid on
    hold off
    legend('Interpreter','latex')

end


% figure(6)
% subplot(3,3,1)
% contour(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)))
% xlim([L_1_l, L_1_u])
% ylim([L_2_l, L_2_u])
% title("Real PDF",'Interpreter','latex')
% xlabel('$x_1$ [m]','Interpreter','latex')
% ylabel('$x_2$ [m]','Interpreter','latex')
% axis equal
% grid on
% legend('$\Phi(\mathbf{x})$','Interpreter','latex','Location','northeastoutside')
% for i = 1:n_iter
%     subplot(3,3,i+1)
%     contour(x_1_grid, x_2_grid, reshape(Phi_hat_x_reg(:,:,i), length(x_2), length(x_1)))
%     xlim([L_1_l, L_1_u])
%     ylim([L_2_l, L_2_u])
%     title("Estimated PDF, iteration " + i,'Interpreter','latex')
%     xlabel('$x_1$ [m]','Interpreter','latex')
%     ylabel('$x_2$ [m]','Interpreter','latex')
%     axis equal
%     grid on
%     hold on
%     plot(X_e_spline_reg(:,1,i), X_e_spline_reg(:,2,i),'LineWidth',1.2)
%     plot(X_e_spline_reg(1,1,i), X_e_spline_reg(1,2,i),'ksq','MarkerSize',7,'LineWidth',1.5)
%     legend('$\hat{\Phi}(\mathbf{x})$','$\mathbf{X_e}(t)$', '$\mathbf{X_e}(0)$',...
%         'Interpreter','latex','Location','northeastoutside')
%     hold off
% end


%% Detailed charts about every iteration
%
% for i = 1:n_iter
% 
%     figure(i+10)
%     subplot(2,3,1);
%     contour(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)))
%     xlim([L_1_l, L_1_u])
%     ylim([L_2_l, L_2_u])
%     title("Real PDF",'Interpreter','latex')
%     xlabel('$x_1$ [m]','Interpreter','latex')
%     ylabel('$x_2$ [m]','Interpreter','latex')
%     axis equal
%     grid on
%     legend('$\Phi(\mathbf{x})$','Interpreter','latex','Location','northeastoutside')
%     subplot(2,3,2);
%     contour(x_1_grid, x_2_grid, reshape(Phi_hat_x_reg(:,:,i), length(x_2), length(x_1)))
%     xlim([L_1_l, L_1_u])
%     ylim([L_2_l, L_2_u])
%     title("Estimated PDF, iteration " + i,'Interpreter','latex')
%     xlabel('$x_1$ [m]','Interpreter','latex')
%     ylabel('$x_2$ [m]','Interpreter','latex')
%     axis equal
%     grid on
%     hold on
%     plot(X_e_reg(:,1,i), X_e_reg(:,2,i),'LineWidth',1.2)
%     plot(X_e_reg(1,1,i), X_e_reg(1,2,i),'ksq','MarkerSize',7,'LineWidth',1.5)
%     legend('$\hat{\Phi}(\mathbf{x})$','$\mathbf{X_e}(t)$', '$\mathbf{X_e}(0)$',...
%         'Interpreter','latex','Location','northeastoutside')
%     hold off
%     subplot(2,3,3);
%     contour(x_1_grid, x_2_grid, reshape(C_x_reg(:,end,i), length(x_2), length(x_1)))
%     xlim([L_1_l, L_1_u])
%     ylim([L_2_l, L_2_u])
%     title("Empirical distribution reconstruction",'Interpreter','latex')
%     xlabel('$x_1$ [m]','Interpreter','latex')
%     ylabel('$x_2$ [m]','Interpreter','latex')
%     axis equal
%     grid on
%     hold on
%     plot(X_e_reg(:,1,i), X_e_reg(:,2,i),'LineWidth',1.2)
%     plot(X_e_reg(1,1,i), X_e_reg(1,2,i),'ksq','MarkerSize',7,'LineWidth',1.5)
%     legend('$C(\mathbf{x})$','$\mathbf{X_e}(t)$', '$\mathbf{X_e}(0)$',...
%         'Interpreter','latex','Location','northeastoutside')
%     hold off
%     subplot(2,3,4);
%     surf(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)))
%     xlim([L_1_l, L_1_u])
%     ylim([L_2_l, L_2_u])
%     %title("Real PDF",'Interpreter','latex')
%     xlabel('$x_1$','Interpreter','latex')
%     ylabel('$x_2$','Interpreter','latex')
%     zlabel('$\Phi(\mathbf{x})$','Interpreter','latex')
%     grid on
%     subplot(2,3,5);
%     surf(x_1_grid, x_2_grid, reshape(Phi_hat_x_reg(:,:,i), length(x_2), length(x_1)))
%     xlim([L_1_l, L_1_u])
%     ylim([L_2_l, L_2_u])
%     %title("Reference PDF",'Interpreter','latex')
%     xlabel('$x_1$','Interpreter','latex')
%     ylabel('$x_2$','Interpreter','latex')
%     zlabel('$\hat{\Phi}(\mathbf{x})$','Interpreter','latex')
%     grid on
%     subplot(2,3,6);
%     surf(x_1_grid, x_2_grid, reshape(C_x_reg(:,end,i), length(x_2), length(x_1)))
%     xlim([L_1_l, L_1_u])
%     ylim([L_2_l, L_2_u])
%     %title("Empirical distribution reconstruction",'Interpreter','latex')
%     xlabel('$x_1$','Interpreter','latex')
%     ylabel('$x_2$','Interpreter','latex')
%     zlabel('$C(\mathbf{x})$','Interpreter','latex')
%     grid on
% 
% end

%% Charts about the independent measurement modeling

if n_iter <= 3
    aux_2 = 1;
elseif n_iter > 3 && n_iter <= 6
    aux_2 = 2;
elseif n_iter > 6 && n_iter <= 9
    aux_2 = 3;
elseif n_iter > 9 && n_iter <= 12
    aux_2 = 4;
elseif n_iter > 12 && n_iter <= 15
    aux_2 = 5;
end

for i = 1:n_iter

    % Model = MdlSol_struct(i).GMModel;
    GMM_dist = MdlSol_struct(i).GMM_ordered;

    TMP = pdf(GMM_dist, Omega);

    grps_plot = "$^" + MdlSol_struct(i).grps + "\hat{V}_k$";

    fig6h = figure(6);
    subplot(aux_2,3,i)
    gscatter(MdlSol_struct(i).Preprocessed_X_e(:,1),...
             MdlSol_struct(i).Preprocessed_X_e(:,2),...
             grps_plot)
    hold on
    % scatter(MdlSol_struct(i).Preprocessed_Data(:,1),...
    %          MdlSol_struct(i).Preprocessed_Data(:,2), ...
    %          5, '+', "black", "MarkerEdgeAlpha", 0.1, 'DisplayName', "Previous"+ newline +"Data")
    contour(x_1_grid, x_2_grid, reshape(TMP, length(x_2), length(x_1)), 'DisplayName', 'GMM')
    plot(GMM_dist.mu(:,1), GMM_dist.mu(:,2), ...
        "rs", "LineWidth", 1.5, 'DisplayName', '$\mu$')
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    title("GMM clustering, iteration " + i,'Interpreter','latex')
    xlabel('$x_1$ [m]','Interpreter','latex')
    ylabel('$x_2$ [m]','Interpreter','latex')
    axis equal
    grid on
    hold off
    legend('Interpreter','latex','Location','best')

    fig7h = figure(7);
    subplot(aux_2,3,i)
    scatter3(X_e_spline_reg(:,1,i), X_e_spline_reg(:,2,i), V_Xe_reg(:,:,i), 10, "black")
    hold on
    scatter3(MdlSol_struct(i).Preprocessed_X_e(:,1),...
        MdlSol_struct(i).Preprocessed_X_e(:,2),...
        MdlSol_struct(i).Preprocessed_V_Xe,...
        15, MdlSol_struct(i).grps, "filled")
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    title("GMM clustering, iteration " + i,'Interpreter','latex')
    xlabel('$x_1$ [m]','Interpreter','latex')
    ylabel('$x_2$ [m]','Interpreter','latex')
    zlabel('$V_k$ [N]','Interpreter','latex')
    grid on
    hold off

    fig8h = figure(8);
    subplot(aux_2,3,i)
    hist3(MdlSol_struct(i).Data_Xe_hist_V,'CDataMode','auto','FaceColor','interp',...
        'Nbins',[length(x_1)-1, length(x_2)-1], "EdgeColor", "none")
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    title("Training Data, iteration " + i,'Interpreter','latex')
    xlabel('$x_1$ [m]','Interpreter','latex')
    ylabel('$x_2$ [m]','Interpreter','latex')
    zlabel('Measurement int','Interpreter','latex')
    grid on

    fig9h = figure();
    subplot(aux_2,3,i)
    surf(x_1_grid, x_2_grid, reshape(TMP, length(x_2), length(x_1)),...
        'EdgeAlpha', 0.4, 'FaceAlpha', 0.7)
    colormap default
    %colormap sky
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    title("GMM fit, iteration " + i,'Interpreter','latex')
    xlabel('$x_1$ [m]','Interpreter','latex')
    ylabel('$x_2$ [m]','Interpreter','latex')
    zlabel('Measurement int','Interpreter','latex')
    grid on

end

set(findall(fig6h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig6h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig6h, "-property", "FontSize"), "FontSize", 16)

set(findall(fig7h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig7h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig7h, "-property", "FontSize"), "FontSize", 16)

set(findall(fig8h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig8h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig8h, "-property", "FontSize"), "FontSize", 16)

set(findall(fig9h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig9h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig9h, "-property", "FontSize"), "FontSize", 16)


%% probabilities on x_j given V_bar plots

fig10h = figure(10);
layout25h = tiledlayout(2,1);
title(layout25h,"PDF estimated for iteration " + (i+1) + " (Results)",...
    'Interpreter', 'latex', "FontSize", 16)

color1d = [0 0.4470 0.7410];
color2d = [0.4940 0.1840 0.5560];
color3d = [0.9290 0.6940 0.1250];

% Tiles Axis x_1
nexttile
yyaxis left
scatter(X_e_spline_reg(:,1,n_iter), V_hat_i_reg(:,1,n_iter), 16, "o", "filled",...
    "MarkerFaceAlpha", 0.6, "MarkerFaceColor", color1d)
hold on
scatter(X_e_spline_reg(:,1,n_iter), V_hat_i_reg(:,2,n_iter), 16, "s", "filled",...
    "MarkerFaceAlpha", 0.6, "MarkerFaceColor", color2d)
scatter(X_e_spline_reg(:,1,n_iter), V_hat_i_reg(:,3,n_iter), 16, "d", "filled",...
    "MarkerFaceAlpha", 0.6, "MarkerFaceColor", color3d)
ylabel('Force [N]','Interpreter','latex')
yyaxis right
plot(x_1,Phi_hat_x_1_reg(:,1,n_iter+1), "Color", color1d, "LineWidth", 2)
plot(x_1,Phi_hat_x_1_reg(:,2,n_iter+1), "Color", color2d, "LineWidth", 2)
plot(x_1,Phi_hat_x_1_reg(:,3,n_iter+1), "Color", color3d, "LineWidth", 2)
if i > 1
xline(Estim_sol(n_iter,1).exp_x1_hat, "-", {"$\hat{\mu} =$ " + Estim_sol(n_iter,1).exp_x1_hat},...
    "Color", color1d, "LineWidth", 2,'Interpreter', 'latex',...
    "LabelOrientation","horizontal", "LabelVerticalAlignment", "top");
xline(Estim_sol(n_iter,2).exp_x1_hat, "-", {"$\hat{\mu} =$ " + Estim_sol(n_iter,2).exp_x1_hat},...
    "Color", color2d, "LineWidth", 2,'Interpreter', 'latex',...
    "LabelOrientation","horizontal", "LabelVerticalAlignment", "middle");
xline(Estim_sol(n_iter,3).exp_x1_hat, "-", {"$\hat{\mu} =$ " + Estim_sol(n_iter,3).exp_x1_hat},...
    "Color", color3d, "LineWidth", 2,'Interpreter', 'latex',...
    "LabelOrientation","horizontal", "LabelVerticalAlignment", "bottom");
end
hold off
xlabel('$x_1$','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
lgd = legend("$^{1}\hat{V}_k(x_1)$", "$^{2}\hat{V}_k(x_1)$", "$^{3}\hat{V}_k(x_1)$", ...
    "$P(\theta_1 ^{x_1})$", "$P(\theta_2 ^{x_1})$", ...
    "$P(\theta_3 ^{x_1})$", 'Location','best','NumColumns',2,'Interpreter','latex');
title(lgd,'Data | Estimated PDF')
grid on
yyaxis left
ylim([0 40])

% Axis x_2
nexttile
yyaxis left
scatter(X_e_spline_reg(:,2,n_iter), V_hat_i_reg(:,1,n_iter), 16, "o", "filled",...
    "MarkerFaceAlpha", 0.6, "MarkerFaceColor", color1d)
hold on
scatter(X_e_spline_reg(:,2,n_iter), V_hat_i_reg(:,2,n_iter), 16, "s", "filled",...
    "MarkerFaceAlpha", 0.6, "MarkerFaceColor", color2d)
scatter(X_e_spline_reg(:,2,n_iter), V_hat_i_reg(:,3,n_iter), 16, "d", "filled",...
    "MarkerFaceAlpha", 0.6, "MarkerFaceColor", color3d)
ylabel('Force [N]','Interpreter','latex')
yyaxis right
plot(x_2,Phi_hat_x_2_reg(:,1,n_iter+1), "Color", color1d, "LineWidth", 2)
plot(x_2,Phi_hat_x_2_reg(:,2,n_iter+1), "Color", color2d, "LineWidth", 2)
plot(x_2,Phi_hat_x_2_reg(:,3,n_iter+1), "Color", color3d, "LineWidth", 2)
if i > 1
xline(Estim_sol(n_iter,1).exp_x2_hat, "-", {"$\hat{\mu} =$ " + Estim_sol(n_iter,1).exp_x2_hat},...
    "Color", color1d, "LineWidth", 2,'Interpreter', 'latex',...
    "LabelOrientation","horizontal", "LabelVerticalAlignment", "top");
xline(Estim_sol(n_iter,2).exp_x2_hat, "-", {"$\hat{\mu} =$ " + Estim_sol(n_iter,2).exp_x2_hat},...
    "Color", color2d, "LineWidth", 2,'Interpreter', 'latex',...
    "LabelOrientation","horizontal", "LabelVerticalAlignment", "middle");
xline(Estim_sol(n_iter,3).exp_x2_hat, "-", {"$\hat{\mu} =$ " + Estim_sol(n_iter,3).exp_x2_hat},...
    "Color", color3d, "LineWidth", 2,'Interpreter', 'latex',...
    "LabelOrientation","horizontal", "LabelVerticalAlignment", "bottom");
end
hold off
xlabel('$x_2$','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
lgd = legend("$^{1}\hat{V}_k(x_2)$", "$^{2}\hat{V}_k(x_2)$", "$^{3}\hat{V}_k(x_2)$", ...
    "$P(\theta_1 ^{x_2})$", "$P(\theta_2 ^{x_2})$", ...
    "$P(\theta_3 ^{x_2})$", 'Location','best','NumColumns',2,'Interpreter','latex');
title(lgd,'Data | Estimated PDF')
grid on
yyaxis left
ylim([0 40])

set(findall(fig10h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig10h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig10h, "-property", "FontSize"), "FontSize", 16)

%% Real PDF  vs  Estimated PDF: Resultados

% GMM_hat = gmdistribution(Mu_found, Sigma_found);
% Phi_hat_final = pdf(GMM_hat, Omega);
% n_def_found = size(Sigma_found, 3);

nbDrawingSeg = 100;
tmp_vec = linspace(-pi, pi, nbDrawingSeg)';
stdev_Phi_hat = zeros(size(Sigma_hat));
Sigma_ast_Phi_hat = zeros(size(Sigma_hat));
Elipse_Phi_hat = zeros(height(tmp_vec), 2, n_def); %Elipse
for j = 1:n_def
    stdev_Phi_hat(:,:,j) = sqrtm(Sigma_hat(:,:,j));
    Sigma_ast_Phi_hat(:,:,j) = 3*stdev_Phi_hat(:,:,j);
    Elipse_Phi_hat(:,:,j) = [cos(tmp_vec), sin(tmp_vec)] * real(Sigma_ast_Phi_hat(:,:,j)) +...
                            repmat(Mu_hat(j,:), nbDrawingSeg, 1);
end

Est_color = [0.8500 0.3250 0.0980];

fig24h = figure(24);
tiledlayout(2,3)

nexttile([2, 2])
contour(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)), contornos)
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Real PDF vs Estimated PDF",'Interpreter','latex')
xlabel('$x_1$ [m]','Interpreter','latex')
ylabel('$x_2$ [m]','Interpreter','latex')
axis equal
grid on
hold on
for j = 1:n_def
    patch(Elipse_Phi_hat(:,1,j), Elipse_Phi_hat(:,2,j), Est_color, 'LineWidth', 1.5,...
        'EdgeColor', Est_color, "FaceAlpha",0.2);
end
for j = 1:n_def
    plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "-.k", "LineWidth", 1)
end
plot(Mu(:,1),Mu(:,2),'.k','MarkerSize',8)
plot(Mu_hat(:,1), Mu_hat(:,2), '+', 'LineWidth', 1.5, 'color', Est_color);
hold off
legend('Real defects, $\Phi(\mathbf{x})$',...
    'Estimation, $\hat{\Phi}(\mathbf{x})$','Location','best')

nexttile(3)
surf(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "EdgeColor", "none", "FaceColor", "interp")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Real PDF")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
grid on

nexttile(6)
surf(x_1_grid, x_2_grid, reshape(Phi_hat_x, length(x_2), length(x_1)),...
    "EdgeColor", "none", "FaceColor", "interp")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Estimated PDF")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
grid on

set(findall(fig24h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig24h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig24h, "-property", "FontSize"), "FontSize", 20)

%% Datos para animación


% save("data_for_animation.mat","t", "n_iter","x_1_grid","x_2_grid",...
%     "x_1", "x_2", "Phi_x", "L_1_l", "L_1_u", "L_2_l", "L_2_u", "X_e_reg",...
%     "Phi_hat_x_reg", "C_x_reg")


%% MPC

% % Registers
% f_k_traj_reg = zeros(K^n, length(t));
% c_k_reg = zeros(K^n, length(t));
% Varepsilon_reg = zeros(length(t), 1);
% Z_reg = zeros(length(t), size(z,1));
% U_reg = zeros(length(t), size(u,1));
% 
% % Initializations
% z_act = z_0;
% c_k_act = c_k_0;
% f_k_traj_act = f_k_traj_0;
% Varepsilon_act = Varepsilon_0;
% 
% for i = 1:length(t)
% 
%     % Optimal control input computation given actual state z
%     u = full( M(z_act) );
% 
%     % Saving in registers
%     f_k_traj_reg(:,i) = f_k_traj_act;
%     c_k_reg(:,i) = c_k_act;
%     Varepsilon_reg(i) =  Varepsilon_act;
%     Z_reg(i,:) = z_act;
%     U_reg(i,:) = u;
% 
%     % Apply optimal control u and simulate system
%     z_act = full( F(z_act,u) );
% 
%     % Compute Fourier Functions and coefficients on the new position
%     X_e = [z_act(1), z_act(3)];   % Position [x_1, x_2]
%     f_k_traj_act = prod(cos( K_cal'.*pi.*(X_e - L_i_l)./(L_i_u - L_i_l) ), 2) ./ h_k_reg ;
%     c_k_act = c_k_act + (f_k_traj_act*T_s)/t_f ;
% 
%     %Ergodic metric
%     Varepsilon_act = sum( Lambda_k .* (c_k_act - phi_k_reg).^2 ); 
% 
% end
% 
% %% Charts MPC results
% 
% figure(3)
% subplot(1,2,1)
% plot(t, Z_reg, 'LineWidth', 1.5)
% title("States",'Interpreter','latex')
% xlabel('Time [s]','Interpreter','latex')
% ylabel('Position [m], velocity[m/s]','Interpreter','latex')
% legend('$x_1$', '$\dot{x}_1$', '$x_2$', '$\dot{x}_2$','Interpreter','latex')
% grid on
% subplot(1,2,2)
% stairs(t, U_reg, 'LineWidth', 1.5)
% title("Control actions",'Interpreter','latex')
% xlabel('Time [s]','Interpreter','latex')
% ylabel('Force [N]','Interpreter','latex')
% legend('$u_1$', '$u_2$','Interpreter','latex')
% grid on
% 
% figure(4)
% plot(t, Varepsilon_reg, "LineWidth",1.5)
% title("Ergodic Metric",'Interpreter','latex')
% xlabel('Time [s]','Interpreter','latex')
% ylabel('$\varepsilon \left( \mathbf{X_e}(t), \Phi(\mathbf{x}) \right) $','Interpreter','latex')
% grid on

