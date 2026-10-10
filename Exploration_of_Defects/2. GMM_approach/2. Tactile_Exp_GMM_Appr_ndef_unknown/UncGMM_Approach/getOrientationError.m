function delta_theta = getOrientationError(Sigma, Sigma_prime)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the acute orientation error between n pair of 
% Ellipses defined by the i-th Gaussian covariance matrix:
% Orientation(Sigma(i)) -- Orientation(Sigma_prime(i)), i = 1,...,n.
%
% INPUTS:
% Sigma:        (2 x 2 x n)     Covariances for the first set of Gaussians
% Sigma_prime : (2 x 2 x n)     Covariances for the second set of Gaussians
%
% OUTPUT:
% delta_theta : (1 x 1)         Acute angle error in DEGREES (0 to 90)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of ellipses
n = size(Sigma, 3);

% Error vector initialization
delta_theta = zeros(n,1);

for i = 1:n

    % Extract components for Matrix 1
    var_x1  = Sigma(1,1,i);
    var_y1  = Sigma(2,2,i);
    cov_xy1 = Sigma(1,2,i);
    
    % Extract components for Matrix 2
    var_x2  = Sigma_prime(1,1,i);
    var_y2  = Sigma_prime(2,2,i);
    cov_xy2 = Sigma_prime(1,2,i);
    
    % Compute individual orientation angles (in radians)
    theta1 = 0.5 * atan2(2 * cov_xy1, var_x1 - var_y1);
    theta2 = 0.5 * atan2(2 * cov_xy2, var_x2 - var_y2);
    
    % Convert to degrees
    theta1_deg = rad2deg(theta1);
    theta2_deg = rad2deg(theta2);
    
    % Compute absolute difference
    raw_diff = abs(theta1_deg - theta2_deg);
    
    % Wrap to a 180-degree period (due to ellipse axial symmetry)
    delta_theta(i) = mod(raw_diff, 180);
    
    % Constrain to acute angle (0 to 90 degrees)
    if delta_theta(i) > 90
        delta_theta(i) = 180 - delta_theta(i);
    end

end



end