function defect = rectangular_def(Omega, C_r, theta_r, rect_lens)
% =========================================================================
% Inputs:
%   C_r         : (n_r x 2) : Centers for n_r rectangles
%   theta_r     : (n_r x 1) : Vector or rotations in degrees
%   rect_lens   : (n_r x 2) : Lengths for n_r rectangles (widths, heights)   
%   Omega       : (n_p x 2) : 2D Space discretized in n_p points
%
% Outputs:
%   defect.pdf  : (n_p x 1) : Defect mask on Omega, 1 values = defective
%                             points, 0 = non-defective
%   defect.Center  : (1 x 2): Geometric center of defective area 
%
% =========================================================================

% number of rectangles
n_rect = height(C_r);

% Initialize matrix of logic values
idx_buff = false(height(Omega), n_rect);

for i = 1:n_rect
    % Rotation matrix
    R_theta = [cosd(theta_r(i)), -sind(theta_r(i)); 
               sind(theta_r(i)), cosd(theta_r(i))];
    
    % Omega points expressed in rentangle local frame
    X_r = (R_theta'*(Omega - C_r(i,:))')';
    
    % Check if points are inside the rectangle
    idx_r = abs(X_r) <= rect_lens(i,:) / 2;
    idx_r = idx_r(:, 1) & idx_r(:, 2);
    
    % Save indexes in matrix
    idx_buff(:,i) = idx_r;
end

% Join the rectangles to have one single component (defect)
idx_joint = any(idx_buff, 2);

% Transform logic vector to double vector
defect.pdf = 1*idx_joint;

% Compute geometric center of defective area
if n_rect > 1
    % Area of each rectangle
    A = prod(rect_lens, 2); 
    % Weighted average for geometric center
    defect.Center = mean(C_r, "Weights", A); 
else
    defect.Center = C_r;
end

end