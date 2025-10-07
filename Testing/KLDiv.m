close all
clear
clc

%% 
n_dim = 2; % Número de dimensiones espaciales
dx = 0.01;

L_1_l = 0.5;
dx_1 = dx;
L_1_u = 1.5;

L_2_l = 0.5;
dx_2 = dx;
L_2_u = 1.5;

% Dimensiones \mathbf{x} = [x_1 x_2]^T
x_1 = (L_1_l:dx_1:L_1_u)';
x_2 = (L_2_l:dx_2:L_2_u)';

% vector de límites inferior y superiores de las dimensiones
L_i_l = [L_1_l, L_2_l];
L_i_u = [L_1_u, L_2_u];

[x_1_grid, x_2_grid] = meshgrid(x_1, x_2);

% Espacio de búsqueda discretizado
Omega = [reshape(x_1_grid,[],1), reshape(x_2_grid,[],1)]; 

%% Uniform PDF as Prior

Phi_hat_x_1 = unifpdf(x_1, L_1_l, L_1_u);
Phi_hat_x_2 = unifpdf(x_2, L_2_l, L_2_u);

[Phi_hat_x_1_grid, Phi_hat_x_2_grid] = meshgrid(Phi_hat_x_1, Phi_hat_x_2);
Q = prod([reshape(Phi_hat_x_1_grid,[],1), reshape(Phi_hat_x_2_grid,[],1)], 2);

%% Gaussian Mixture distribution as posterior

% Means 
Mu = [1.3938    0.7218;	    
      0.7558    0.8584;	    
      1.2323    1.2161];

% Matriz de Covarianza
Var_def = 0.01;
% Var_def = 0.0005;

Cov_1 = [Var_def, Var_def/2;
         Var_def/2, Var_def];
Cov_2 = [Var_def, -Var_def/2;
         -Var_def/2, Var_def];
Cov_3 = Cov_1;

Sigma_P1 = cat(3,Cov_1,Cov_2,Cov_3);

n_components = size(Mu, 1);

gm_dist_P1 = gmdistribution(Mu, Sigma_P1);% , proporciones);

%PDF de referencia REAL
P_1 = pdf(gm_dist_P1, Omega);

%% Gaussian Mixture distribution as Thiner posterior

% Matriz de Covarianza
Var_def = 0.0005;

Cov_1 = [Var_def, Var_def/2;
         Var_def/2, Var_def];
Cov_2 = [Var_def, -Var_def/2;
         -Var_def/2, Var_def];
Cov_3 = Cov_1;

Sigma_P2 = cat(3,Cov_1,Cov_2,Cov_3);

gm_dist_P2 = gmdistribution(Mu, Sigma_P2);% , proporciones);

%PDF de referencia REAL
P_2 = pdf(gm_dist_P2, Omega);

%% KL Divergence Test on Overall Distribution

% The KL Divergence is computed with the Riemann numerical integral
% \Int(f(x))dx \approx sum( f(x_i) )*dx

% Calculate KL Divergence from Q to P_1
idx_P1_Q = (Q ~= 0) & (P_1 ~= 0);
D_KL_P1_Q = sum(P_1(idx_P1_Q) .* log(P_1(idx_P1_Q) ./ ...
                Q(idx_P1_Q)))*dx_1*dx_2;

% Calculate KL Divergence from P_1 to P_2
idx_P2_P1 = (P_1 ~= 0) & (P_2 ~= 0);
D_KL_P2_P1 = sum(P_2(idx_P2_P1) .* log(P_2(idx_P2_P1) ./ ...
                 P_1(idx_P2_P1)))*dx_1*dx_2;

% Calculate KL Divergence from Q to P_2
idx_P2_Q = (Q ~= 0) & (P_2 ~= 0);
D_KL_P2_Q = sum(P_2(idx_P2_Q) .* log(P_2(idx_P2_Q) ./ ...
                Q(idx_P2_Q)))*dx_1*dx_2;

%% KL Divergence Independently (General Form)

D_KL_P1_Q_ind = zeros(1, n_components);
for i = 1:n_components
    P_1_i = pdf(gmdistribution(Mu(i,:), Sigma_P1(:,:,i)), Omega);
    idx_P1_Q_i = (Q ~= 0) & (P_1_i ~= 0);
    D_KL_P1_Q_ind(i) = sum(P_1_i(idx_P1_Q_i) .* ...
                       log(P_1_i(idx_P1_Q_i) ./ ...
                       Q(idx_P1_Q_i)))*dx_1*dx_2;
end

D_KL_P2_P1_ind = zeros(1, n_components);
for i = 1:n_components
    P_1_i = pdf(gmdistribution(Mu(i,:), Sigma_P1(:,:,i)), Omega);
    P_2_i = pdf(gmdistribution(Mu(i,:), Sigma_P2(:,:,i)), Omega);
    idx_P2_P1_i = (P_1_i ~= 0) & (P_2_i ~= 0);
    D_KL_P2_P1_ind(i) = sum(P_2_i(idx_P2_P1_i) .* ...
                        log(P_2_i(idx_P2_P1_i) ./ ...
                        P_1_i(idx_P2_P1_i)))*dx_1*dx_2;
end

D_KL_P2_Q_ind = zeros(1, n_components);
for i = 1:n_components
    P_2_i = pdf(gmdistribution(Mu(i,:), Sigma_P2(:,:,i)), Omega);
    idx_P2_Q_i = (Q ~= 0) & (P_2_i ~= 0);
    D_KL_P2_Q_ind(i) = sum(P_2_i(idx_P2_Q_i) .* ...
                       log(P_2_i(idx_P2_Q_i) ./ ...
                       Q(idx_P2_Q_i)))*dx_1*dx_2;
end

%% KL Divergence Independently (Closed Form for Gaussian Distributions)

D_KL_P2_P1_ind_Gaus = zeros(1, n_components);
Prev_Mu = Mu;
for i = 1:n_components
    D_KL_P2_P1_ind_Gaus(i) = (1/2)*( trace( Sigma_P1(:,:,i)\Sigma_P2(:,:,i) ) ...
                        - n_dim +...
            (Prev_Mu(i,:) - Mu(i,:))*(Sigma_P1(:,:,i)\(Prev_Mu(i,:) - Mu(i,:))') + ...
            log(det(Sigma_P1(:,:,i)) / det(Sigma_P2(:,:,i))) );
end



%% Graficar

fig1h = figure(1);
layouth = tiledlayout(fig1h, 1, 3);

nexttile(layouth)
surfc(x_1_grid, x_2_grid, reshape(Q, length(x_2), length(x_1)),...
    "EdgeColor","none","FaceColor","interp", "FaceAlpha",0.3)
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Prior Q")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
axis tight
grid on
legend('$Q(\mathbf{x})$','Location','northeastoutside')

nexttile(layouth)
surfc(x_1_grid, x_2_grid, reshape(P_1, length(x_2), length(x_1)),...
    "EdgeColor","none","FaceColor","interp","FaceAlpha",0.3)
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Posterior $P_1$")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
axis tight
grid on
legend('$P_1(\mathbf{x})$','Location','northeastoutside')

nexttile(layouth)
surfc(x_1_grid, x_2_grid, reshape(P_2, length(x_2), length(x_1)),...
    "EdgeColor","none","FaceColor","interp","FaceAlpha",0.3)
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Thinner Posterior $P_2$")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
axis tight
grid on
legend('$P_2(\mathbf{x})$','Location','northeastoutside')

set(findall(fig1h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig1h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig1h, "-property", "FontSize"), "FontSize", 18)


