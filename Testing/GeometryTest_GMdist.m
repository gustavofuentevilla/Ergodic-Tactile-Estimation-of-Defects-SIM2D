
close all
clearvars -except M
% clearvars -except M %run_idx
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

Cov_1 = [Var_def, 0;
         0, Var_def];
Cov_2 = [Var_def, 0
        0, Var_def + 0.001];
Cov_3 = [Var_def + 0.001, -Var_def/2;
         -Var_def/2, Var_def];

Sigma = cat(3,Cov_1,Cov_2,Cov_3);

gm_dist = gmdistribution(Mu, Sigma);

%PDF
Phi_x = pdf(gm_dist, Omega);

n_def = size(Sigma, 3);

%% Geometric computations of Gaussian Mixture Distribution Phi

stdev_Phi = zeros(size(Sigma));
Sigma_ast_Phi = zeros(size(Sigma));
r_elips_Phi = zeros(n, n_def);
variation_Phi = zeros(1, n_def);
for j = 1:n_def
    % Standard deviation
    stdev_Phi(:,:,j) = sqrtm(Sigma(:,:,j)); 
    % 3*Standard deviation that represents 99% of data
    Sigma_ast_Phi(:,:,j) = 3*stdev_Phi(:,:,j); 
    % radios de ejes principales 
    r_elips_Phi(:,j) = eig(Sigma_ast_Phi(:,:,j));
    % Variación total = trace(Sigma_ast_phi) = sum(r_elips_phi)
    variation_Phi(:,j) = sum(r_elips_Phi(:,j));
end

%% Adding variance and compute geometrics
offset_axis = 0.1; %10 cm

% Augmented covariance matrix (the one we want to compute)
Sigma_a = zeros(size(Sigma)); 
for i = 1:n_def
    % Eigenvectors = V, Eigenvalues diagonal matriz = D
    [V, D] = eig(Sigma_ast_Phi(:,:,i));
    % Sorting with Max eigenvalue at last spot
    [r_j, ind] = sort(diag(D)); %r_j is the vector of elipse radius = eigenvalues
    D_sorted = D(ind, ind); 
    V_sorted = V(:,ind);

    % Respect aspect ratio (r_2 = Ratio_12*r_1)
    r_extension = zeros(size(r_j));
    Ratio_12 = r_j(end) / r_j(end-1);
    % Defining offset over minimum axis and computing the rest
    r_extension(1) = offset_axis; 
    r_extension(end) = Ratio_12*r_extension(end-1);

    % New diagonal eigenvalues matriz with entended radius
    D_a = D_sorted + diag(r_extension);

    % Augmented covariance matrix
    Sd_a = V_sorted*D_a*V_sorted' / 3;
    Sigma_a(:,:,i) = Sd_a * Sd_a;

end

% Sin sentido físico, no respeta el aspect ratio
% Sigma_a = Sigma + diag([offset_axis, offset_axis]);

%% Geometrics of augmented covariance matrix Sigma_a
stdev_a = zeros(size(Sigma_a));
Sigma_ast_a = zeros(size(Sigma_a));
r_elips_a = zeros(n, n_def);
variation_a = zeros(1, n_def);
for j = 1:n_def
    stdev_a(:,:,j) = sqrtm(Sigma_a(:,:,j)); % Standard deviation
    Sigma_ast_a(:,:,j) = 3*stdev_a(:,:,j);
    % radios de ejes principales 
    r_elips_a(:,j) = eig(Sigma_ast_a(:,:,j));
    % Variación total = trace(Sigma_ast_phi) = sum(r_elips_phi)
    variation_a(:,j) = sum(r_elips_a(:,j));
end

%%

% For plotting
%
nbDrawingSeg = 1000;
tmp_vec = linspace(-pi, pi, nbDrawingSeg)';
Elipse_Phi = zeros(height(tmp_vec), 2, n_def); % Elipse
Elipse_a = zeros(height(tmp_vec), 2, n_def);
for j = 1:n_def  
    Elipse_Phi(:,:,j) = [cos(tmp_vec), sin(tmp_vec)] * real(Sigma_ast_Phi(:,:,j)) + repmat(Mu(j,:),nbDrawingSeg,1);
    Elipse_a(:,:,j) = [cos(tmp_vec), sin(tmp_vec)] * real(Sigma_ast_a(:,:,j)) + repmat(Mu(j,:),nbDrawingSeg,1);
end

figure(1)
contour(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)), 25)
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Real PDF",'Interpreter','latex')
xlabel('$x_1$ [m]','Interpreter','latex')
ylabel('$x_2$ [m]','Interpreter','latex')
axis equal
grid on
hold on
for j = 1:n_def
    plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "k", "LineWidth",1.3)
    plot(Elipse_a(:,1,j), Elipse_a(:,2,j), "r", "LineWidth",1.3)
end
plot(Mu(:,1), Mu(:,2), ".")
hold off

figure(2)
surf(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
        'EdgeColor','interp','FaceColor','interp', "FaceAlpha",0.2)
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
hold on
for j = 1:n_def
    plot3(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), zeros(size(Elipse_Phi(:,1,j))), "k", "LineWidth", 2)
    plot3(Elipse_a(:,1,j), Elipse_a(:,2,j), zeros(size(Elipse_a(:,1,j))), "r", "LineWidth", 2)
end
hold off


