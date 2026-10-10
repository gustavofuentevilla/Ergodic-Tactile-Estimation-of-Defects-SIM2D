close all
clear
clc

idx_X0 = 1;
run_idx = 1;

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

xrange = [L_i_l', L_i_u'];

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

n_def = 3;

% Save random number generator to reproduce Tests
rng_saved = rng;

% Defects generator
[Mu, Sigma, r_elips_Phi] = DefectsGen(n_def, L_i_l, L_i_u);

% proportions = [1/(1+2+3),...
%                2/(1+2+3),...
%                3/(1+2+3)];
gm_dist = gmdistribution(Mu, Sigma);%, proportions);

%PDF de referencia REAL
Phi_x = pdf(gm_dist, Omega);

%% Parametros

% interpolation
% freq_interp = 200;
% T_s_interp = 1/freq_interp;

% constant velocity
vel = 0.05;

%% Measurement parameters
a = 5; %Reference force
b = 0.2;
c = 0.05; % noise amplitud

% V_real = a + b*Phi_x;

% Initial position
X_e = X0(idx_X0,:);

% X_e_sp = X_e;

% Take the measurements
delta = c*randn(height(X_e), 1);
V = a + b*pdf(gm_dist, X_e) + delta;

% Execution time (during motion)
t_exec = 0;

%% Load traveled distanced from TEDD

Case_i = 1; % TEDD
Test_i = 3; % 2: equal severities, 3: different severities

% Compute TEDD's traveled distance
filename = "/home/gustavo-fuentevilla/MATLAB/" + ...
           "Ergodic_Tactile_Estimation_of_Defects_SIM2D/" + ...
           "Exploration_of_Defects/2. GMM_approach/" + ...
           "2. Tactile_Exp_GMM_Appr_ndef_unknown/" + ...
           "Test" + Test_i + "/Case" + Case_i +...
           "/" + n_def + "Def/X0_" + idx_X0 + ...
           "/output_" + run_idx + ".mat";
load(filename, "n_iter", "X_e_spline_reg")

TEDD_TravDist = 0;
for j = 1:n_iter
    P_diff = diff(X_e_spline_reg(:,:,j));
    TEDD_TravDist = TEDD_TravDist + ...
                        sum( sqrt( P_diff(:,1).^2 + P_diff(:,2).^2 ) );

end

%% Loop

timer_iter = [];
samples = 0;
Rand_TravDist = 0;

while Rand_TravDist <= TEDD_TravDist

    samples = samples + 1;
    disp(samples)

    timer_init = tic;

    % Next random point in Omega
    X_e(end + 1,:) = L_i_l + (L_i_u - L_i_l).*rand(1,2);

    % distance from last sample to current
    dist = pdist(X_e(end-1:end,:));

    % Estimated travel time with constant velocity
    t_traj = dist/vel;

    % Execution time vector
    t_exec = [t_exec; t_exec(end) + t_traj];

    % Make interpolation
    % t_interp = (t_exec(end-1):T_s_interp:t_exec(end)-T_s_interp)';
    % t_interp = (t_exec(end-1):T_s_interp:t_exec(end))';
    % x_e_sp = spline(t_exec(end-1:end), X_e(end-1:end,1), t_interp);
    % y_e_sp = spline(t_exec(end-1:end), X_e(end-1:end,2), t_interp);
    % X_e_sp_new = [x_e_sp, y_e_sp];
    % X_e_sp(end+1:end+height(X_e_sp_new), :) = X_e_sp_new;
    
    % Update the measurement vector
    delta_new = c*randn(height(X_e(end,:)), 1);
    V_new = a + b*pdf(gm_dist, X_e(end,:)) + delta_new;
    delta = [delta; delta_new];
    V = [V; V_new];

    % Update the traveled distance after adding the new sample
    Rand_TravDist = Rand_TravDist + dist;

    timer_iter(end + 1, :) = toc(timer_init);
    

end

%% More plots xd

t = linspace(0, t_exec(end), height(X_e))';
% v_normplot = sqrt(sum((diff(X_e_sp)/T_s_interp).^2, 2));
% v_normplot(end + 1) = v_normplot(end);

fig2h = figure(2);
tiledlayout(3,5)

nexttile(1, [1, 3])
plot(t, X_e, "LineWidth", 3)
% plot(t_exec, X_e, "LineWidth", 3)
xlabel('Time [s]')
ylabel('Position [m]')
title("Position")
legend("$x_{e_1}$", "$x_{e_2}$")

% nexttile(2)
% plot(t, v_normplot, "LineWidth", 3)
% xlabel('Time [s]')
% ylabel('Velocity [m/s]')
% title('Velocity norm')
% grid on;

nexttile(4, [2,2])
plot(X_e(:,1), X_e(:,2), "k-*", "LineWidth", 2, "MarkerSize", 8)
hold on
plot(X_e(1,1), X_e(1,2), 'rsq', 'LineWidth',4, 'MarkerSize',12)
plot(X_e(end,1), X_e(end,2), "o", "MarkerFaceColor", "r", "MarkerSize", 12)
hold off
axis equal tight
xlabel("$x_1$")
ylabel("$x_2$")
title("Exploration trajectory")

set(findall(fig2h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig2h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig2h, "-property", "FontSize"), "FontSize", 20)

%% Post-processing

Par_PDF.thres_meas = a + max(delta);
Par_PDF.D_max = 10;
Par_PDF.Thres_Variation = max(sum(r_elips_Phi));
Par_PDF.OneClustDistLimit = 2*max(r_elips_Phi,[],"all") + 0.07;

Estim_sol = Postprocessing(X_e, V, Par_PDF);

Priors_found = Estim_sol.Priors_found;
Mu_found = Estim_sol.Mu_found;
Sigma_found = Estim_sol.Sigma_found_a;
