function [Sigma_a, flag_readjust] = ReadjustConstrainedSigma(Sigma,...
                                    MinVariation, MinAxisLengths)

% Sigma_a: Augmented Sigma, the one to compute
% flag_readjust: indicates if any component have been readjusted or not

% Number of components
nbEllipses = size(Sigma, 3);
% Augmented covariance matrix (the one we want to compute)
Sigma_a = zeros(size(Sigma));
% Flag set to True if that component has been readjusted
flag_readjust = false(1, nbEllipses);
% Computation of current axis lengths
stdev = zeros(size(Sigma));
Sigma_ast = zeros(size(Sigma));
Variation_Sigma = zeros(1, nbEllipses);

for j = 1:nbEllipses
    % Standard deviation
    stdev(:,:,j) = sqrtm(Sigma(:,:,j)); 
    % 3*Standard deviation that represents 99% of data
    Sigma_ast(:,:,j) = 3*stdev(:,:,j); 
    % Eigenvectors = V, Eigenvalues diagonal matriz = D
    [V, D] = eig(Sigma_ast(:,:,j));
    % Sorting with Max eigenvalue at last spot
    % r_j is the vector of ellipse radius of Sigma on the plane
    [r_j, ind] = sort(diag(D)); 
    D_sorted = D(ind, ind); 
    V_sorted = V(:,ind);
    % Compute Ratio of original Sigma if aspect ratio is wanted to be
    % mantained, otherwise set Ratio to 1
    % Ratio = r_j(2)/r_j(1);
    Ratio = 1;
    % Variación total = trace(Sigma_ast_phi) = sum(r_j)
    Variation_Sigma(:,j) = sum(r_j);

    % Total Variation Constraint:
    % If the sum of axis lengths of Sigma is less than the 
    % limit (MinVariation), then re-adjust to achieve the minimal

    if (Variation_Sigma(:,j) < MinVariation)
        flag_readjust(j) = true;
        % Compute the variation needed to achieve the minimal
        DeltaVariation = MinVariation - Variation_Sigma(:,j);
        % Compute Extension to be made over axes mantaining the aspect 
        % ratio (r_2_bar = Ratio*r_1_bar, r_1_bar < r_2_bar)
        r_bar = zeros(size(r_j));
        % Compute axis lengths offsets to achieve the MinVariation
        r_bar(1) = DeltaVariation/(1 + Ratio); 
        r_bar(2) = Ratio*r_bar(1);
        % New diagonal eigenvalues matriz with extended radius
        D_a = D_sorted + diag(r_bar);
        % Augmented covariance matrix
        Sd_a = V_sorted*D_a*V_sorted' / 3;
        Sigma_a(:,:,j) = Sd_a * Sd_a;
        % Update lengths (V_sorted is the same as original Sigma)
        D_sorted = D_a; 
        r_j = diag(D_a);
    end

    % Minimum axis lengths constraint:

    if r_j(1) < MinAxisLengths
        flag_readjust(j) = true;
        % Compute Extension to be made over axes mantaining the aspect 
        % ratio (r_2_bar = Ratio*r_1_bar, r_1_bar < r_2_bar)
        r_bar = zeros(size(r_j));
        % Compute axis lengths offsets to achieve the MinAxisLength
        r_bar(1) = MinAxisLengths - r_j(1); 
        r_bar(2) = Ratio*r_bar(1);
        % New diagonal eigenvalues matriz with extended radius
        D_a = D_sorted + diag(r_bar);
        % Augmented covariance matrix
        Sd_a = V_sorted*D_a*V_sorted' / 3;
        Sigma_a(:,:,j) = Sd_a * Sd_a;
    end
    
    % If there was no re-adjustment, then the updated Sigma_a 
    % is the original Sigma
    if ~flag_readjust(j)
        Sigma_a(:,:,j) = Sigma(:,:,j);
    end

end