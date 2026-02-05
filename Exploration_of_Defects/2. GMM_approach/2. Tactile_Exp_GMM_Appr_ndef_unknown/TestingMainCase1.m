% Primera ejecución
% import casadi.*
% Erg_traj_ipopt = Function.load('/home/gustavo-fuentevilla/MATLAB/Ergodic_Tactile_Estimation_of_Defects_SIM2D/Casadi_Formulation_ExplTask/Erg_traj_2.casadi');

for idx_X0 = 1:20

for run_idx = 1:5

% Resto de ejecuciones
close all
clearvars -except Erg_traj_ipopt run_idx idx_X0
clc

import casadi.*

%% Parámetros del espacio de búsqueda U = [L_1_l, L_1_u] \times [L_2_l, L_2_u]

n = 2; % Número de dimensiones espaciales

L_1_l = 1;
L_1_u = 1.5;

L_2_l = 1;
L_2_u = 1.5;

L_1 = (L_1_u - L_1_l);
L_2 = (L_2_u - L_2_l);

dx_1 = L_1/100;
dx_2 = L_2/100;

% Dimensiones \mathbf{x} = [x_1 x_2]^T
x_1 = (L_1_l:dx_1:L_1_u)';
x_2 = (L_2_l:dx_2:L_2_u)';

%vector de límites inferior y superiores de las dimensiones
L_i_l = [L_1_l, L_2_l];
L_i_u = [L_1_u, L_2_u];

[x_1_grid, x_2_grid] = meshgrid(x_1, x_2);

%Espacio de búsqueda discretizado
Omega = [reshape(x_1_grid,[],1), reshape(x_2_grid,[],1)]; 

%% Initial positions vector

X0 = [L_1_l, L_2_l;
      L_1_l + 1*L_1/5, L_2_l;
      L_1_l + 2*L_1/5, L_2_l;
      L_1_l + 3*L_1/5, L_2_l;
      L_1_l + 4*L_1/5, L_2_l;
      L_1_u, L_2_l;
      L_1_u, L_2_l + 1*L_2/5;
      L_1_u, L_2_l + 2*L_2/5;
      L_1_u, L_2_l + 3*L_2/5;
      L_1_u, L_2_l + 4*L_2/5;
      L_1_u, L_2_u;
      L_1_u - 1*L_1/5, L_2_u;
      L_1_u - 2*L_1/5, L_2_u;
      L_1_u - 3*L_1/5, L_2_u;
      L_1_u - 4*L_1/5, L_2_u;
      L_1_l, L_2_u;
      L_1_l, L_2_u - 1*L_2/5;
      L_1_l, L_2_u - 2*L_2/5;
      L_1_l, L_2_u - 3*L_2/5;
      L_1_l, L_2_u - 4*L_2/5];

%% Defects definition

n_def = 5;

% Save random number generator to reproduce Tests
rng_saved = rng;

% Defects generator
[Mu, Sigma, r_elips_Phi] = DefectsGen(n_def, L_i_l, L_i_u);

proportions = [1/(1+1+2+2+3),...
               1/(1+1+2+2+3),...
               2/(1+1+2+2+3),...
               2/(1+1+2+2+3),...
               3/(1+1+2+2+3)];
gm_dist = gmdistribution(Mu, Sigma, proportions);

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
% R = [7e-5, 0;
%      0, 7e-5]*(1/T_s); % N = 200

% Peso sobre métrica ergódica
% gamma = 1;

% Estado inicial z = [z_1; z_2; z_3; z_4] = [x_1; x_1_dot; x_2; x_2_dot]
z_0 = [X0(idx_X0,1); 0; X0(idx_X0,2); 0]; 

%Pre-cálculo de Lambda
p = 2; %norma 2
Lambda_k = (1 + vecnorm(K_cal, p, 1)').^(-(n + 1)/2);

%% vector to add more points on the trajectory and get more data from sensor

freq_spline = 200;
t_spline = (0:1/freq_spline:t_f)';

%% Loop for the Search task
n_iter_max = 20;

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
c = 0.05; %Amplitud del ruido en la medición
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
Par_PDF.Thres_Variation = max(sum(r_elips_Phi)) + 0.001;
% Minimum axes lengths of gaussian elipses
Par_PDF.MinAxisLengths = 0; % 0 cm (Not using this constrain)
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

% Times
T_ErgC_i = zeros(n_iter_max, 1);
T_PDF_i = zeros(n_iter_max, 1);

% Time for the whole Simulation
T_sim_tic = tic;

for i = 1:n_iter_max
    
    t_erg_init = tic;
    [Z, U] = Erg_traj_ipopt(z_act, phi_k_act); % Trayectoria Ergodica
    T_ErgC_i(i) = toc(t_erg_init);

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

    t_pdf_init = tic;
    [Phi_hat_x_next, Estim_sol(i)] = PDF_EstimatorCase1(X_e_spline, V_Xe, Par_PDF);
    T_PDF_i(i) = toc(t_pdf_init);

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

T_sim_toc = toc(T_sim_tic);


%% Saving variables (Except the casadi function)

% Normal for-loop

save(sprintf("Test3/Case1/%dDef/X0_%d/output_%d.mat",...
             n_def, idx_X0, run_idx),...
     "-regexp", "^(?!(Erg_traj_ipopt)$).");


end

end


