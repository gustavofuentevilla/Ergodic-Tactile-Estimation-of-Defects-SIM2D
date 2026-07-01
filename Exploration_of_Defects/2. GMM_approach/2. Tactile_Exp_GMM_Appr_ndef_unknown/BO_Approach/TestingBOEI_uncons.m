% Script for testing the BO approach through several simulations

% Test2 : Equal Severity levels
% Test3 : Different severity levels
% Case1 : Full algorithm
% Case2 : Without repeating points

for idx_X0 = 1:20

for run_idx = 1:5

close all
clearvars -except run_idx idx_X0
clc

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

% vector de límites inferior y superiores de las dimensiones
L_i_l = [L_1_l, L_2_l];
L_i_u = [L_1_u, L_2_u];

xrange = [L_i_l', L_i_u'];

[x_1_grid, x_2_grid] = meshgrid(x_1, x_2);

% Espacio de búsqueda discretizado
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

% Load random number generator to reproduce Tests (from Case 1)
load(sprintf("/home/gustavo-fuentevilla/MATLAB/"...
    + "Ergodic_Tactile_Estimation_of_Defects_SIM2D/"...
    + "Exploration_of_Defects/2. GMM_approach/"...
    + "2. Tactile_Exp_GMM_Appr_ndef_unknown/"...
    + "Test3/Case1/%dDef/X0_%d/output_%d.mat",...
             n_def, idx_X0, run_idx), "rng_saved")
rng(rng_saved)

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

%% Tiempos

N_MaxSamples = 200;

% t_f = 20;           %Tiempo final por iteración
% T_s = 0.1;                  % Tiempo de muestreo
% t = (0:T_s:t_f)';   %Vector de tiempo por iteración

% t for interpolation
freq_interp = 200;
T_s_interp = 1/freq_interp;
% t_interp = (0:T_s_interp:T_s)';

% constant velocity between samples
vel = 0.5;
% Max displacement
displacement = vel*T_s_interp;

%% Measurement parameters
a = 5; %Reference force
b = 0.2;
c = 0.05; %Amplitud del ruido en la medición 0.05

V_real = a + b*Phi_x;

% Initial position
X_e = X0(idx_X0,:);

% Generate a random first displacement within the velocity constraint
% dmax = v_max*T_s; % Maximum displacement
% validPoint = false;
% 
% while ~validPoint
% 
% Xrand = [(X_e(1,1) - dmax) + 2*dmax*rand, ...
%           (X_e(1,2) - dmax) + 2*dmax*rand];
% 
% v_init2 = sum(((Xrand - X_e)/T_s).^2, 2);
% 
% IsVelValid = v_init2 <= v_max^2;
% IsInOmega = all(isbetween(Xrand, L_i_l, L_i_u));
% 
% validPoint = IsVelValid && IsInOmega;
% 
% end

% Generate a random first displacement within the maximum velocity
% theta = linspace(0, 2*pi, 100)';
% P_tmp = [X_e(1,1) + dmax*cos(theta), X_e(1,2) + dmax*sin(theta)];
% IsInOmega = false;
% while ~IsInOmega
%     Xrand = P_tmp(randi(height(P_tmp)),:);
%     IsInOmega = all(isbetween(Xrand, L_i_l, L_i_u));
% end

% Initialize trajectory wit zero velocity and a random first displacement
% X_e = [X_e; X_e]; % Xrand];

% Interpolate to a 200Hz trajectory
% x_e_sp = spline(t(1:height(X_e)), X_e(:,1), t_interp);
% y_e_sp = spline(t(1:height(X_e)), X_e(:,2), t_interp);
%
% X_e_sp = [x_e_sp, y_e_sp];

X_e_sp = X_e;

% Take the measurements
delta = c*randn(height(X_e_sp), 1);
V = a + b*pdf(gm_dist, X_e_sp) + delta;

t_exec = 0;

%% Loop

timer_iter = zeros(N_MaxSamples, 1);
loglikelihood = [];

for i = height(X_e):N_MaxSamples

    % avgV = mean(V);
    % stdV = std(V);
    % 
    % if stdV <= 1e-8 || isnan(stdV)
    %     stdV = 1;
    % end
    % 
    % V_normalized = (V - avgV)/stdV;
    % mdl = fitrgp(X_e_sp, V_normalized);

    timer_init = tic;

    mdl = fitrgp(X_e_sp, V,...
                'KernelFunction', 'squaredexponential',...
                'Standardize',true,...
                'SigmaLowerBound',0.1,...
                'BasisFunction','constant',...
                'PredictMethod','exact');
    [V_pred, sd] = predict(mdl, Omega);

    %% Expected Improvement
    % This EI is from http://krasserm.github.io/2018/03/21/bayesian-optimization/
    
    xi = 0;  % Exploration-exploitation parameter (greek letter, xi)
             % High xi = more exploration
             % Low xi = more exploitation (can be < 0)

    d = V_pred - max(V) - xi; % (y - f*) if maximization

    EI = (sd ~= 0).*(d.*normcdf(d./sd) + sd.*normpdf(d./sd));

    % sort EI (first element is the maximum)
    [eisorted, idx_ei] = sort(EI, 1, "descend");
    % re-arrange Omega
    OmegaSortedEI = Omega(idx_ei,:);

    % Algorithm for selecting the next point X
    for j = 1:height(eisorted)
        % save the maximum EI (first element)
        eimax = eisorted(j);
        % check for equal EI values
        idxEqEI = eimax == eisorted;
        % Extract the set of posible x's in which EI is maximum
        posibleX = OmegaSortedEI(idxEqEI, :);
        % Select a random X from the set
        idx_rand = randi(height(posibleX));
        xEI = posibleX(idx_rand, :);
        % check for already visited point, if already visited, iterate
        % again
        flag_visitedP = ismember(xEI, X_e, "rows");
        if ~flag_visitedP
            break;
        end
    end

    X_e(end + 1,:) = xEI;

    % distance from current BO sample to the next
    dist = pdist(X_e(end-1:end,:));
    % Estimated travel time with constant velocity
    t_traj = dist/vel;
    % Execution time vector
    t_exec = [t_exec; t_exec(end) + t_traj];
    % Make interpolation
    % t_interp = (t_exec(end-1):T_s_interp:t_exec(end)-T_s_interp)';
    t_interp = (t_exec(end-1):T_s_interp:t_exec(end))';
    x_e_sp = spline(t_exec(end-1:end), X_e(end-1:end,1), t_interp);
    y_e_sp = spline(t_exec(end-1:end), X_e(end-1:end,2), t_interp);
    X_e_sp_new = [x_e_sp, y_e_sp];
    X_e_sp(end+1:end+height(X_e_sp_new), :) = X_e_sp_new;
    % Update the measurement vector
    delta_new = c*randn(height(X_e_sp_new), 1);
    V_new = a + b*pdf(gm_dist, X_e_sp_new) + delta_new;
    delta = [delta; delta_new];
    V = [V; 
        V_new];

    % Update the log-likelihood

    loglikelihood = [loglikelihood;...
                    mdl.LogLikelihood];

    timer_iter(i) = toc(timer_init);

    % Exit the loop when the LogLikelihood is low enough
    if mdl.LogLikelihood <= -1000
        disp("LogLikelihood reached the condition")
        samples = i;
        break;
    end

    samples = i

    mdl.LogLikelihood

end

%% Post-processing

timer_PDF_init = tic;

Par_PDF.thres_meas = a + max(delta);
Par_PDF.D_max = 10;
Par_PDF.Thres_Variation = max(sum(r_elips_Phi));
Par_PDF.OneClustDistLimit = 2*max(r_elips_Phi,[],"all") + 0.07;

Estim_sol = BO_Postprocessing(X_e_sp, V, Par_PDF);

Priors_found = Estim_sol.Priors_found;
Mu_found = Estim_sol.Mu_found;
Sigma_found = Estim_sol.Sigma_found_a;

timer_PDF = toc(timer_PDF_init);

%% Save test

save(sprintf("/home/gustavo-fuentevilla/MATLAB/"...
    + "Ergodic_Tactile_Estimation_of_Defects_SIM2D/"...
    + "Exploration_of_Defects/2. GMM_approach/"...
    + "2. Tactile_Exp_GMM_Appr_ndef_unknown/"...
    + "BO_Approach/BO_unconsTest/Test3/Case1/%dDef/X0_%d/output_%d.mat",...
             n_def, idx_X0, run_idx));

end
end
