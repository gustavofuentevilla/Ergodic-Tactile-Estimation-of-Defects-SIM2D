close all
clearvars -except M %run_idx
clc

%% Parámetros del espacio de búsqueda U = [L_1_l, L_1_u] \times [L_2_l, L_2_u]

for k = 1:10

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

%% Gaussian Mixture distribution (Real PDF)

% Número de defectos a generar
n_def = 1;

% Registro de Covarianzas genéricas (3)
Var_def = 0.0005;

Gen_Cov = zeros(n,n,3);
Gen_Cov(:,:,1) = [Var_def, Var_def/2;
                  Var_def/2, Var_def];
Gen_Cov(:,:,2) = [Var_def, -Var_def/2;
                  -Var_def/2, Var_def];
Gen_Cov(:,:,3) = [Var_def, 0;
                  0, Var_def];

idx_gencov = (1:size(Gen_Cov,3))'; 
idx_gencov = repmat(idx_gencov, 5*size(Gen_Cov,3), 1);

% Covarianzas para los defectos
Sigma = zeros(n,n,n_def);
for i = 1:n_def
    Sigma(:,:,i) = Gen_Cov(:,:,idx_gencov(i));
end

% Geometric computations of Covariance Ellipses
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
r_elips_sorted = sort(reshape(r_elips_Phi,[], 1));

% Generate Means randomly
Verify_Mu_dist = true;
while(Verify_Mu_dist)

    offset = 0.05; 
    Mu = [];
    for i = 1:n_def
        Mu_tmp = (L_i_l + offset) + ((L_i_u - offset) - ...
                 (L_i_l + offset)).*rand(1,2);
        Mu = cat(1, Mu, Mu_tmp); 
    end

    if n_def == 1
        break;
    end
    
    % Calculate the pairwise distances between the means
    Mu_dist = pdist(Mu);
    
    % Distance beetween Mu's should be less than the two largest ellipse
    % axes to ensure that defects doesn't overlap
    Verify_Mu_dist = any(Mu_dist < ...
                        (r_elips_sorted(end) + r_elips_sorted(end-1)));

end

% Generate distribution
gm_dist = gmdistribution(Mu, Sigma);

% PDF de referencia REAL
Phi_x = pdf(gm_dist, Omega);

%% Charts

nbDrawingSeg = 100;
tmp_vec = linspace(-pi, pi, nbDrawingSeg)';
Elipse_Phi = zeros(height(tmp_vec), 2, n_def); %Elipse
for j = 1:n_def
    Elipse_Phi(:,:,j) = [cos(tmp_vec), sin(tmp_vec)] * real(3*sqrtm(Sigma(:,:,j))) + repmat(Mu(j,:),nbDrawingSeg,1);
end

figh = figure;
pcolor(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
       "FaceColor","interp","EdgeColor","none")
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
    plot(Elipse_Phi(:,1,j), Elipse_Phi(:,2,j), "-.", "LineWidth", 1.3)
end
plot(Mu(:,1), Mu(:,2), ".", 'MarkerSize', 8)
hold off
legend('$\Phi(\mathbf{x})$','Location', 'best') %'northeastoutside')

set(findall(figh,'-property','Interpreter'),'Interpreter','latex') 
set(findall(figh,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(figh, "-property", "FontSize"), "FontSize", 16)

end