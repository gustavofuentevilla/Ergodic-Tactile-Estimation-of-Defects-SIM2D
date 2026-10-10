function score = getEllipseIoU(Mu, Sigma, Mu_prime, Sigma_prime, points)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Ellipse Intersection-over-Union metric of n
% pair of Ellipses defined by the i-th Gaussian parameters:
% Ellipse(Mu(i), Sigma(i)) -- Ellipse(Mu_prime(i), Sigma_prime(i)),
% i = 1,...,n.
%
% INPUTS
% Mu:           (n x 2)         2D Gaussian centroids (n Gaussians)
% Sigma:        (2 x 2 x n)     Covariances related to Mu
% Mu_prime:     (n x 2)         Centroids of Gaussians to compare with
% Sigma_prime:  (2 x 2 x n)     Covariances related to Mu_prime
% points:       (1 x 1)         Number of points to generate the ellipses
%
% OUTPUTS
% score:        (n x 1)         Ellipse IoU metrics of n Gaussians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of ellipses to match
n = height(Mu);

% Parametric vector
t = linspace(-pi, pi, points)';

% Coeficient for standard deviation (a*std)
a = 3;

% Initialize
Ellipse = zeros(length(t), 2, n);
Ellipse_prime = zeros(length(t), 2, n);
score = zeros(n,1);

for i = 1:n

    % Compute ellipse points
    Ellipse(:,:,i) = [cos(t), sin(t)] * ...
                     real(a*sqrtm(Sigma(:,:,i))) + ...
                     repmat(Mu(i,:),points,1);
    Ellipse_prime(:,:,i) = [cos(t), sin(t)] * ...
                           real(a*sqrtm(Sigma_prime(:,:,i))) + ...
                           repmat(Mu_prime(i,:),points,1);

    % Compute polygons
    ellipsePolygon = polyshape(Ellipse(:,:,i));
    ellipsePrimePolygon = polyshape(Ellipse_prime(:,:,i));

    % Compute intersection and union areas
    intersectionArea = area(intersect(ellipsePolygon, ellipsePrimePolygon));
    unionArea = area(union(ellipsePolygon, ellipsePrimePolygon));

    % Compute Ellipse IoU score
    score(i,1) = intersectionArea / unionArea;

end



end