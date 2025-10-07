close all
clear
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
Mu_1 = [1.2, 0.7];    
Mu_2 = [0.7, 1];    
Mu_3 = [1.2, 1.3];

Mu = [Mu_1; Mu_2; Mu_3];

% Matriz de Covarianza
Var_def = 0.0005;

Cov_1 = [Var_def, Var_def/2;
         Var_def/2, Var_def];
Cov_2 = [Var_def, 0
        0, Var_def + 0.001];
Cov_3 = [Var_def + 0.001, -Var_def/2;
         -Var_def/2, Var_def];

Sigma = cat(3,Cov_1,Cov_2,Cov_3);

gm_dist = gmdistribution(Mu, Sigma);

%PDF
Phi_x = pdf(gm_dist, Omega);

n_def = size(Sigma, 3); %número de gaussianos


%% Geometric computations of Gaussian Mixture Distribution Phi

stdev_Phi = zeros(size(Sigma));
Sigma_ast_Phi = zeros(size(Sigma));
r_elips_Phi = zeros(n, n_def);
variation_Phi = zeros(1, n_def);
isInElipse = false([height(Omega), n_def]);
for j = 1:n_def
    % Standard deviation
    stdev_Phi(:,:,j) = sqrtm(Sigma(:,:,j)); 
    % 3*Standard deviation that represents 99% of data
    Sigma_ast_Phi(:,:,j) = 3*stdev_Phi(:,:,j); 

    [V, D] = eig(Sigma_ast_Phi(:,:,j));
    [lambdas_s, ind] = sort(diag(D));
    D_sorted = D(ind, ind); %El ultimo eigenvalor es el mayor
    V_sorted = V(:,ind);

    P_ast = V_sorted'*( Omega - Mu(j,:) )'; %Omega = vector de datos
    
    Elips_eq = sum( (P_ast.^2) ./ (lambdas_s.^2), 1 )';

    isInElipse(:,j) = Elips_eq <= 1; %Indice de los puntos que están dentro de las elipses
end

%índice de todos los datos a eliminar (los que están dentro de las elipses)
idx_erase = any(isInElipse, 2); 

%Eliminación de datos dentro de las elipses
Omega(idx_erase,:) = [];




%%

% For plotting

nbDrawingSeg = 1000;
tmp_vec = linspace(-pi, pi, nbDrawingSeg)';
Elipse_Phi = zeros(height(tmp_vec), 2, n_def); % Elipse
for j = 1:n_def  
    Elipse_Phi(:,:,j) = [cos(tmp_vec), sin(tmp_vec)] * real(Sigma_ast_Phi(:,:,j)) + repmat(Mu(j,:),nbDrawingSeg,1);
end

figure(1)
pcolor(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)), ...
    "FaceColor","interp","EdgeColor","none")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Real PDF",'Interpreter','latex')
xlabel('$x_1$ [m]','Interpreter','latex')
ylabel('$x_2$ [m]','Interpreter','latex')
axis equal tight
colormap sky
grid on
hold on
for j = 1:n_def
    plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "k", "LineWidth",2)
end
plot(Mu(:,1), Mu(:,2), ".")
scatter(Omega(:,1), Omega(:,2), 18, "filled", "o")
hold off

figure(2)
surf(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
        'EdgeColor','interp','FaceColor','interp', "FaceAlpha",0.2)
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
hold on
for j = 1:n_def
    plot3(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), zeros(size(Elipse_Phi(:,1,j))), "k", "LineWidth", 2)
end
hold off


