function [Mu, Sigma, r_elips_Phi] = DefectsGen(n_def, L_i_l, L_i_u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% n_def: (1 x 1) Número de defectos a generar
% L_i_l: (1 x 2) Límites inferiores del espacio de búsqueda [L_1_l, L_2_l]
% L_i_u: (1 x 2) Límites superiores del espacio de búsqueda [L_1_u, L_2_u]
% OUTPUTS
% Mu: (n_def x 2) Bidimensional means, each row represents a defect center
% Sigma: (2 x 2 x n_def) Defect Shapes Covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 2; % Dimension of Search Space


% Registro de Covarianzas genéricas (3 formas de defectos)
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

end