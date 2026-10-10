% Script for testing the BO-TD (Traveled distance) approach through several simulations

% Test2 : Equal Severity levels
% Test3 : Different severity levels
% Case1 : Full algorithm (TEDD)

for n_def = 3:2:7

for idx_X0 = 1:20

for run_idx = 1:5

close all
clearvars -except run_idx idx_X0 n_def
clc

Case_i = 1; % TEDD
Test_i = 3; % 2: equal severities, 3: different severities

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

% n_def = 3;

% Load variables from TEDD to reproduce
filename = "/home/gustavo-fuentevilla/MATLAB/" + ...
           "Ergodic_Tactile_Estimation_of_Defects_SIM2D/" + ...
           "Exploration_of_Defects/2. GMM_approach/" + ...
           "2. Tactile_Exp_GMM_Appr_ndef_unknown/" + ...
           "Test" + Test_i + "/Case" + Case_i +...
           "/" + n_def + "Def/X0_" + idx_X0 + ...
           "/output_" + run_idx + ".mat";
load(filename, "n_iter", "X_e_spline_reg", "rng_saved")

rng(rng_saved)

% Defects generator
[Mu, Sigma, r_elips_Phi] = DefectsGen(n_def, L_i_l, L_i_u);

if n_def == 3
    proportions = [1/(1+2+3),...
                   2/(1+2+3),...
                   3/(1+2+3)];
elseif n_def == 5
    proportions = [1/(1+1+2+2+3),...
                   1/(1+1+2+2+3),...
                   2/(1+1+2+2+3),...
                   2/(1+1+2+2+3),...
                   3/(1+1+2+2+3)];
elseif n_def == 7
    proportions = [1/(1+1+1+2+2+2+3),...
                   1/(1+1+1+2+2+2+3),...
                   1/(1+1+1+2+2+2+3),...
                   2/(1+1+1+2+2+2+3),...
                   2/(1+1+1+2+2+2+3),...
                   2/(1+1+1+2+2+2+3),...
                   3/(1+1+1+2+2+2+3)];
end

gm_dist = gmdistribution(Mu, Sigma, proportions);

% gm_dist = gmdistribution(Mu, Sigma);

%PDF de referencia REAL
Phi_x = pdf(gm_dist, Omega);

%% Parametros

% interpolation
freq_interp = 200;
T_s_interp = 1/freq_interp;

% constant velocity
vel = 0.1;

%% Measurement parameters
a = 5; %Reference force
b = 0.2;
c = 0.05; % noise amplitud

% Initial position
X_e = X0(idx_X0,:);

X_e_sp = X_e;

% Take the measurements
delta = c*randn(height(X_e), 1);
V = a + b*pdf(gm_dist, X_e) + delta;

% Execution time (during motion)
t_exec = 0;

% Compute traveled distanced from TEDD
TEDD_TravDist = 0;
for j = 1:n_iter
    P_diff = diff(X_e_spline_reg(:,:,j));
    TEDD_TravDist = TEDD_TravDist + ...
                        sum( sqrt( P_diff(:,1).^2 + P_diff(:,2).^2 ) );
end

%% Approximate h for serpentine trajectory given the amount of distance

% Offset on the borders of Omega
alpha = 0.0;

% Corners starting point detection
C_DownLeft = isequal(X_e, [L_1_l, L_2_l]);
C_DownRight = isequal(X_e, [L_1_u, L_2_l]);
C_TopRight = isequal(X_e, [L_1_u, L_2_u]);
C_TopLeft = isequal(X_e, [L_1_l, L_2_u]);

Corners = C_DownLeft || C_DownRight || C_TopRight || C_TopLeft;

% Sides starting point detection
Bottom = X_e(2) == L_2_l;
Right = X_e(1) == L_1_u;
Top = X_e(2) == L_2_u;
Left = X_e(1) == L_1_l;

if Corners
    d = 0;
elseif Bottom
    d = pdist2(X_e, [L_1_l, L_2_l]);
elseif Right
    d = pdist2(X_e, [L_1_u, L_2_l]);
elseif Top
    d = pdist2(X_e, [L_1_u, L_2_u]);
elseif Left
    d = pdist2(X_e, [L_1_l, L_2_u]);
end

h = (L_1 - 2*alpha)*(L_2 - 2*alpha) / (TEDD_TravDist + d - L_2 + alpha);

% Approximatete ciclos number
Nc = (L_2 - 2*alpha)/h;

% Assure the trajectory won't go outside the domain
h = h - 1e-6;

%% Loop

flag_init = true;
timer_iter = [];
samples = 0;
Serp_TravDist = 0;

while Serp_TravDist <= TEDD_TravDist

    samples = samples + 1;
    disp(samples)

    timer_init = tic;

    %% Next point in Omega logic

    % logic starting from the bottom
    if Bottom && ~C_DownRight
        if flag_init
            % First step: Displacement X+, Y remains
            X_e(end + 1,:) = [L_1_u, X_e(end,2)]; 
            flag_init = false;
        else
            % If last 2 Y coordinates didn't change
            if X_e(end,2) == X_e(end-1,2)
                % Displacement upwards
                X_e(end + 1,:) = X_e(end,:) + [0, h]; 
            else
                % If they did, then check X coordinate to move
                if X_e(end,1) == L_1_l
                    % Displacement X+, Y remains
                    X_e(end + 1,:) = [L_1_u, X_e(end,2)]; 
                else
                    % Displacement X-, Y remains
                    X_e(end + 1,:) = [L_1_l, X_e(end,2)]; 
                end
            end
        end
        

    % logic starting from the right
    elseif Right && ~C_TopRight
        if flag_init
            % First step: Displacement Y+, X remains
            X_e(end + 1,:) = [X_e(end,1), L_2_u];
            flag_init = false;
        else
            % If last 2 X coordinates didn't change
            if X_e(end,1) == X_e(end-1,1)
                % Displacement to the left
                X_e(end + 1,:) = X_e(end,:) - [h, 0]; 
            else
                % If they did, then check Y coordinate to move
                if X_e(end,2) == L_2_l
                    % Displacement Y+, X remains
                    X_e(end + 1,:) = [X_e(end,1), L_2_u]; 
                else
                    % Displacement Y-, X remains
                    X_e(end + 1,:) = [X_e(end,1), L_2_l]; 
                end
            end
        end

    % logic starting from top
    elseif Top && ~C_TopLeft
        if flag_init
            % First step: Displacement X-, Y remains
            X_e(end + 1,:) = [L_1_l, X_e(end,2)]; 
            flag_init = false;
        else
            % If last 2 Y coordinates didn't change
            if X_e(end,2) == X_e(end-1,2)
                % Displacement downwards
                X_e(end + 1,:) = X_e(end,:) - [0, h]; 
            else
                % If they did, then check X coordinate to move
                if X_e(end,1) == L_1_l
                    % Displacement X+, Y remains
                    X_e(end + 1,:) = [L_1_u, X_e(end,2)]; 
                else
                    % Displacement X-, Y remains
                    X_e(end + 1,:) = [L_1_l, X_e(end,2)]; 
                end
            end
        end


    % logic starting from the left
    elseif Left && ~C_DownLeft
        if flag_init
            % First step: Displacement Y-, X remains
            X_e(end + 1,:) = [X_e(end,1), L_2_l];
            flag_init = false;
        else
            % If last 2 X coordinates didn't change
            if X_e(end,1) == X_e(end-1,1)
                % Displacement to the right
                X_e(end + 1,:) = X_e(end,:) + [h, 0]; 
            else
                % If they did, then check Y coordinate to move
                if X_e(end,2) == L_2_l
                    % Displacement Y+, X remains
                    X_e(end + 1,:) = [X_e(end,1), L_2_u]; 
                else
                    % Displacement Y-, X remains
                    X_e(end + 1,:) = [X_e(end,1), L_2_l]; 
                end
            end
        end



    end

    % distance from last sample to current
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
    V = [V; V_new];

    % Update the traveled distance after adding the new sample
    Serp_TravDist = Serp_TravDist + dist;

    timer_iter(end + 1, :) = toc(timer_init);
    

end


%% Post-processing

timer_PDF_init = tic;

Par_PDF.thres_meas = a + max(delta);
Par_PDF.D_max = 10;
% Par_PDF.Thres_Variation = max(sum(r_elips_Phi));
Par_PDF.OneClustDistLimit = 2*max(r_elips_Phi,[],"all") + 0.07;

Estim_sol = Postprocessing(X_e_sp, V, Par_PDF);

Priors_found = Estim_sol.Priors_found;
Mu_found = Estim_sol.Mu_found;
Sigma_found = Estim_sol.Sigma_found;

timer_PDF = toc(timer_PDF_init);

%% Save test

filename = "/home/gustavo-fuentevilla/MATLAB/" + ...
           "Ergodic_Tactile_Estimation_of_Defects_SIM2D/" + ...
           "Exploration_of_Defects/2. GMM_approach/" + ...
           "2. Tactile_Exp_GMM_Appr_ndef_unknown/" + ...
           "Serpentine_Approach/Serp_TD_Test/Test" + Test_i +...
           "/" + n_def + "Def/X0_" + idx_X0 + ...
           "/output_" + run_idx + ".mat";
save(filename)

end
end
end