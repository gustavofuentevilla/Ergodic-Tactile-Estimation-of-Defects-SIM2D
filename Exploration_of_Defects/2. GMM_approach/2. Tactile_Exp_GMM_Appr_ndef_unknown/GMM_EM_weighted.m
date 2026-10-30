function GMModel = GMM_EM_weighted(Data, Weights, K, Opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted Gaussian Mixture Model fitting with EM algorithm
%
% Inputs:
%   Data                (N x D)         Data with N observations and D 
%                                       dimensions.
% 
%   Weights             (N x 1)         Weights column vector.
%   K                   (1 x 1)         Number of clusters.
%
% Options (Opts):
%   Max_iter            (1 x 1)         Maximum number of iterations
%   loglik_threshold    (1 x 1)         Log-Likelihood threshold
%   Min_var             (1 x 1)         Minimum variation on Sigma
%
% Outputs:
%   GMModel:            (struct)
%       Mu              (K x D)         Cluster centroids
%       Sigma           (D x D x K)     Covariance matrixes
%       Priors          (1 x K)         Mixture coefficients
%       BIC             (1 x 1)         Bayesian Information Criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    Data (:,:) {mustBeNonempty, mustBeNumeric}
    Weights (:,:) {mustBeNonempty, mustBeNumeric}
    K    (1,1) {mustBeInteger, mustBePositive}
    Opts.Max_iter (1,1) {mustBeInteger, mustBePositive} = 100
    Opts.loglik_threshold (1,1) {mustBeNumeric, mustBePositive} = 1e-10
    Opts.Min_var (1,1) {mustBeNumeric, mustBePositive} = 1e-5
end

% Extract Data and normalize weights
X = Data;
W = Weights;
W = W(:) / sum(W);

[nbData, nbDim] = size(X);
loglik_old = -realmax;
nbStep = 0;

%% Initialization with weighted K-means
[Data_id, Mu0] = kmeans_weighted(X, K, W,"replicates",3);

Priors0 = zeros(1, K);
Sigma0 = zeros(nbDim, nbDim, K);

for i = 1:K
    idtmp = find(Data_id == i);
    if isempty(idtmp)
        Priors0(i) = 1 / K;
        Sigma0(:,:,i) = eye(nbDim);
        continue;
    end
    
    % Sum of weights in current cluster
    Priors0(i) = sum(W(idtmp)); 
    
    % Initial weighted covariance matrix for current cluster
    w_cluster = W(idtmp) / sum(W(idtmp));
    X_c = X(idtmp, :) - Mu0(i, :);
    Sigma0(:,:,i) = (X_c' * (X_c .* w_cluster));
end

% Add tiny variance to avoid numerical issues
Sigma0 = Sigma0 + eye(nbDim)*Opts.Min_var;
Priors0 = Priors0 / sum(Priors0);

Mu = Mu0;
Sigma = Sigma0;
Priors = Priors0;

%% Expectation-Maximization (EM) Algorithm
while true

    % ---- E-step ----
    Pxi = zeros(nbData, K);
    for i = 1:K
        % Compute probability p(x|i)
        % Probability of having x due to Gaussian i
        Pxi(:,i) = gaussPDF(X, Mu(i,:), Sigma(:,:,i));
    end
    % To avoid zeros or NaNs
    Pxi = Pxi + realmin;

    % Compute posterior probability p(i|x): bayes rule
    % Probability that a given x belongs to Gaussian i
    Pix_tmp = Pxi .* Priors; 
    Pix = Pix_tmp ./ sum(Pix_tmp, 2); % (nbData, K)

    % Weighted responsabilities
    Pix_w = Pix .* W; % (nbData, K)
    % Weighted cumulative probability
    E = sum(Pix_w, 1); % (1, K)

    % ---- M-step ----
    for i = 1:K
        % Priors update
        Priors(i) = E(i); 

        % Means update
        Mu(i,:) = (Pix_w(:,i)' * X) / E(i);

        % Covariance update
        Data_Mu = X - Mu(i,:);
        Sigma(:,:,i) = (Data_Mu' * (Data_Mu .* Pix_w(:,i))) / E(i);
    end
    
    % Adding tiny variance
    Sigma = Sigma + eye(nbDim)*Opts.Min_var;

    % Stopping criteria using Weighted Log-likelihood
    for i = 1:K
        Pxi(:,i) = gaussPDF(X, Mu(i,:), Sigma(:,:,i));
    end
    
    F = Pxi * Priors';
    F(F < realmin) = realmin;
    
    loglik = sum(W .* log(F)); 

    nbStep = nbStep + 1;
    if abs((loglik / loglik_old) - 1) < Opts.loglik_threshold
        disp("Iterations in E-M: " + nbStep);
        break;
    end
    loglik_old = loglik;

    if nbStep >= Opts.Max_iter
        warning("Maximum number of iterations in E-M reached");
        break;
    end

end

% Adding tiny variance
Sigma = Sigma + eye(nbDim)*Opts.Min_var;

%% Weighted BIC
p = K * (nbDim * (nbDim + 1) / 2 + nbDim + 1) - 1;
BIC = -2 * sum(W .* log(F)) + p * log(nbData);

% Model output
GMModel.Mu = Mu;
GMModel.Sigma = Sigma;
GMModel.Priors = Priors;
GMModel.BIC = BIC;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Local Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function prob = gaussPDF(Data, Mu, Sigma)
    [nbData, nbDim] = size(Data);
    Data_Mu = Data - Mu;
    prob = sum((Data_Mu / Sigma) .* Data_Mu, 2);
    prob = exp(-0.5 * prob) / sqrt((2*pi)^nbDim * (abs(det(Sigma)) + realmin));
end

