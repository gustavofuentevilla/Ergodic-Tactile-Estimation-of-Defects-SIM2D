close all
clearvars -except Erg_traj_ipopt
clc

import casadi.*

% Load the Casadi object if it's not already loaded
% Erg_traj_ipopt = Function.load('Erg_traj_ipopt.casadi');

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

%% Building of random testing distributions (6 distributions in total)

Phi_hat_x_1 = unifpdf(x_1, L_1_l, L_1_u);
Phi_hat_x_2 = unifpdf(x_2, L_2_l, L_2_u);

[Phi_hat_x_1_grid, Phi_hat_x_2_grid] = meshgrid(Phi_hat_x_1, Phi_hat_x_2);
Phi_hat_x = prod([reshape(Phi_hat_x_1_grid,[],1),...
    reshape(Phi_hat_x_2_grid,[],1)], 2);

% Number of Gaussian components per distribution
n_def = 3;

for i = 41 : -10 : 1

    % Random means for each component
    Mu = L_i_l + (L_i_u - L_i_l).*rand(n_def, 2); 

    % Covariance Matrices
    Var_def = 0.0005*i;
    Cov_1 = [Var_def, 0;
        0, Var_def];
    Cov_2 = [Var_def, -Var_def/2;
        -Var_def/2, Var_def];
    Cov_3 = [Var_def, Var_def/2;
        Var_def/2, Var_def];

    Sigma = cat(3,Cov_1,Cov_2,Cov_3);

    % Distributions set
    gm_dist = gmdistribution(Mu, Sigma);
    Phi = pdf(gm_dist, Omega);
    Phi_hat_x = cat(2, Phi_hat_x, Phi);

end

%% Visualize distributions
% 
% dist = 6; % change this to visualize different distributions within the set
% fig1h = figure(1);
% pcolor(x_1_grid, x_2_grid,...
%        reshape(Phi_hat_x(:,dist), length(x_2), length(x_1)),...
%        "EdgeColor","interp")
% colorbar
% xlim([L_1_l, L_1_u])
% ylim([L_2_l, L_2_u])
% xlabel('$x_1$ [m]','Interpreter','latex')
% ylabel('$x_2$ [m]','Interpreter','latex')
% axis equal tight
% grid on
% legend('$\hat{\Phi}(\mathbf{x})$','Location','northeastoutside')
% 
% set(findall(fig1h,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(fig1h,'-property','TickLabelInterpreter'), ...
%     'TickLabelInterpreter','latex')
% set(findall(fig1h, "-property", "FontSize"), "FontSize", 18)

%% Parameters for coefficients computation

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

%% Optimal problem parameters

N = 200;             % Number of points in trajectory
t_f = 10;            % Time window
T_s = t_f/N;         % Sampling time
t = (0:T_s:t_f)';    % Time vector

% Initial state z = [z_1; z_2; z_3; z_4] = [x_1; x_1_dot; x_2; x_2_dot]
z_0 = [1.0; 0; 1.0; 0]; 

% Lambda pre-computation
p = 2;
Lambda_k = (1 + vecnorm(K_cal, p, 1)').^(-(n + 1)/2);

%% Loop for Testing the Ergodic trajectory Casadi Object

n_iter = size(Phi_hat_x, 2);

% Registers
z_reg = zeros(N+1, 4, n_iter);
u_reg = zeros(N, 2, n_iter);
X_e_reg = zeros(N+1, 2, n_iter);
X_e_dot_reg = zeros(N+1, 2, n_iter);
phi_k_REG = zeros(K^n, 1, n_iter);

% Initial condition
z_act = z_0;

for i = 1:n_iter

    % Fourier Coefficients
    [phi_k_reg, f_k_reg, h_k_reg] = FourierCoef_RefPDF(Phi_hat_x(:,i), Par_struct);

    % Current coefficients
    phi_k_act = phi_k_reg;

    [Z, U] = Erg_traj_ipopt(z_act, phi_k_act);
    Z = full(Z)';
    U = full(U)';

    % Update states
    z_act = Z(end,:)';

    % Trajectory X_e(t) = [x_1, x_2]
    X_e = [Z(:, 1), Z(:, 3)];
    X_e_dot = [Z(:, 2), Z(:, 4)];

    % Registers
    z_reg(:,:,i) = Z;
    u_reg(:,:,i) = U;
    X_e_reg(:,:,i) = X_e;
    X_e_dot_reg(:,:,i) = X_e_dot;
    phi_k_REG(:,:,i) = phi_k_reg; 

end


%% Empirical distribution and ergodic metric reconstruction

Varepsilon_reg = zeros(length(t), 1, n_iter);
C_x_reg = zeros(height(Omega), length(t), n_iter);

for r = 1:n_iter
    c_k = zeros(K^n, 1);
    C_x = zeros(height(Omega), 1);
    for i = 1:length(t)

        % Compute Fourier Functions and coefficients on the new position
        f_k_traj = prod(cos( K_cal'.*pi.*(X_e_reg(i,:,r) - L_i_l)./...
                        (L_i_u - L_i_l) ), 2) ./ h_k_reg ;
        c_k = c_k + (f_k_traj*T_s)/(t_f) ;

        % Ergodic metric
        Varepsilon = sum( Lambda_k .* (c_k - phi_k_REG(:,:,r)).^2 );

        % Empirical distribution
        C_x_i = zeros(height(Omega), 1);
        for j = 1:K^n
            C_x_i = C_x_i + c_k(j)*f_k_reg(:,j);
        end

        C_x = C_x + C_x_i;

        % Registers
        C_x_reg(:,i,r) = C_x;
        Varepsilon_reg(i,:,r) = Varepsilon;

    end

end


%% Pre-processing for Plotting

t_total = zeros(n_iter*(length(t)-1) + 1, 1);
X_e_total = zeros(n_iter*(length(t)-1) + 1, 2);
X_e_dot_total = zeros(n_iter*(length(t)-1) + 1, 2);
Varepsilon_total = zeros(n_iter*(length(t)-1) + 1, 1);
u_total = zeros(n_iter*(length(t)-1) + 1, 2);

for i = 1:n_iter

    id_init = ((i - 1)*(length(t)-1) + 1);
    id_last = (i*(length(t)-1) + 1);       

    t_total( id_init:id_last ) = (i - 1)*t(end) + t;
    X_e_total( id_init:id_last, : ) = X_e_reg(:,:,i);
    X_e_dot_total( id_init:id_last, : ) = X_e_dot_reg(:,:,i);
    Varepsilon_total( id_init:id_last, : ) = Varepsilon_reg(:,:,i);
    u_total( id_init:id_last, : ) = [u_reg(:,:,i); [NaN, NaN]];

end

%% Charts

fig2h = figure(2);
subplot(3,1,1)
plot(t_total, X_e_total, 'LineWidth', 1.5)
title("Position States")
xlabel('Time [s]')
ylabel('Position [m]')
legend('$x_1$', '$x_2$')
grid on
subplot(3,1,2)
plot(t_total, X_e_dot_total, 'LineWidth', 1.5)
title("Velocity States")
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('$\dot{x}_1$', '$\dot{x}_2$')
grid on
subplot(3,1,3)
plot(t_total, u_total, 'LineWidth', 1.5)
title("Control actions")
xlabel('Time [s]')
ylabel('Force [N]')
legend('$u_1$', '$u_2$')
grid on

fig3h = figure(3);
plot(t_total, Varepsilon_total, "LineWidth",1.5)
title("Ergodic Metric")
xlabel('Time [s]')
ylabel('$\varepsilon \left( \mathbf{X_e}(t), \Phi(\mathbf{x}) \right) $')
grid on

fig4h = figure(4);
for i = 1:n_iter
    subplot(2,3,i)
    pcolor(x_1_grid, x_2_grid,...
        reshape(Phi_hat_x(:,i), length(x_2), length(x_1)),...
        "EdgeColor","none","FaceColor","interp")
    xlim([L_1_l, L_1_u])
    ylim([L_2_l, L_2_u])
    xlabel('$x_1$ [m]')
    ylabel('$x_2$ [m]')
    axis tight equal
    grid on
    hold on
    plot(X_e_reg(:,1,i), X_e_reg(:,2,i),'LineWidth',1.2, 'Color','white')
    plot(X_e_reg(1,1,i), X_e_reg(1,2,i),'wsq','MarkerSize',7,'LineWidth',1.5)
    legend('$\hat{\Phi}(\mathbf{x})$','$\mathbf{X_e}(t)$', '$\mathbf{X_e}(0)$',...
        'Location','northeastoutside')
    hold off
end

set(findall(fig2h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig2h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig2h, "-property", "FontSize"), "FontSize", 18)

set(findall(fig3h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig3h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig3h, "-property", "FontSize"), "FontSize", 18)

set(findall(fig4h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig4h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig4h, "-property", "FontSize"), "FontSize", 18)