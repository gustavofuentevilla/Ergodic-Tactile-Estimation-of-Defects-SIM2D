function GMModel = GMM_EM_weighted(Data, Weights, K, Opts)
% GMM_EM_WEIGHTED Ajuste de Gaussian Mixture Model con pesos por muestra.
%
% Entradas:
%   Data -    (N x 2) 2D Data.
%   Weights - (N x 1) Weights vector
%   K    -    (1 x 1) number of clusters.
%   Opts - Options: Max_iter, loglik_threshold, Min_var.

arguments
    Data (:,:) {mustBeNonempty, mustBeNumeric}
    Weights (:,:) {mustBeNonempty, mustBeNumeric}
    K    (1,1) {mustBeInteger, mustBePositive}
    Opts.Max_iter (1,1) {mustBeInteger, mustBePositive} = 100
    Opts.loglik_threshold (1,1) {mustBeNumeric, mustBePositive} = 1e-10
    Opts.Min_var (1,1) {mustBeNumeric, mustBePositive} = 1e-5
end

% Separación de coordenadas espaciales (X) y pesos (W)
X = Data;
W = Weights;
W = W(:) / sum(W); % Normalization

[nbData, nbDim] = size(X);
loglik_old = -realmax;
nbStep = 0;

%% 1. Inicialización con K-Means Ponderado
[Data_id, Mu0, ~] = kmeans_weighted(X, K, W,"replicates",3);

Priors0 = zeros(1, K);
Sigma0 = zeros(nbDim, nbDim, K);

for i = 1:K
    idtmp = find(Data_id == i);
    if isempty(idtmp)
        Priors0(i) = 1 / K;
        Sigma0(:,:,i) = eye(nbDim);
        continue;
    end
    
    % Suma de pesos en el clúster
    Priors0(i) = sum(W(idtmp)); 
    
    % Covarianza ponderada inicial del clúster
    w_cluster = W(idtmp) / sum(W(idtmp));
    X_c = X(idtmp, :) - Mu0(i, :);
    Sigma0(:,:,i) = (X_c' * (X_c .* w_cluster));
end

% Regularización e igualación de Priors
Sigma0 = Sigma0 + repmat(eye(nbDim) * Opts.Min_var, [1, 1, K]);
Priors0 = Priors0 / sum(Priors0);

Mu = Mu0;
Sigma = Sigma0;
Priors = Priors0;

%% 2. Algoritmo Expectation-Maximization (EM)
while true

    % ---- E-step ----
    Pxi = zeros(nbData, K);
    for i = 1:K
        Pxi(:,i) = gaussPDF(X, Mu(i,:), Sigma(:,:,i));
    end
    Pxi = Pxi + realmin; % Prevenir ceros y NaN

    % Probabilidad a posteriori p(i|x) (Regla de Bayes)
    Pix_tmp = Pxi .* Priors; 
    Pix = Pix_tmp ./ sum(Pix_tmp, 2); %(nbData, K)

    % Responsabilidades ponderadas por los pesos de muestra
    Pix_w = Pix .* W; % (nbData, K)
    E = sum(Pix_w, 1); % (1, K) Probabilidad acumulada ponderada

    % ---- M-step ----
    for i = 1:K
        % Actualización de Priors
        Priors(i) = E(i); 

        % Actualización de Centros (Medias)
        Mu(i,:) = (Pix_w(:,i)' * X) / E(i);

        % Actualización de Matrices de Covarianza
        Data_Mu = X - Mu(i,:);
        Sigma(:,:,i) = (Data_Mu' * (Data_Mu .* Pix_w(:,i))) / E(i);
    end
    
    % Regularización de covarianza
    Sigma = Sigma + repmat(eye(nbDim) * Opts.Min_var, [1, 1, K]);

    % ---- Criterio de Parada (Log-Likelihood Ponderado) ----
    for i = 1:K
        Pxi(:,i) = gaussPDF(X, Mu(i,:), Sigma(:,:,i));
    end
    
    F = Pxi * Priors';
    F(F < realmin) = realmin;
    
    % Log-likelihood ponderado por W
    loglik = sum(W .* log(F)); 

    nbStep = nbStep + 1;
    if abs((loglik / loglik_old) - 1) < Opts.loglik_threshold
        disp("Iteraciones en E-M: " + nbStep);
        break;
    end
    loglik_old = loglik;

    if nbStep >= Opts.Max_iter
        warning("Máximo número de iteraciones alcanzado");
        break;
    end

end

% Regularización final
Sigma = Sigma + repmat(eye(nbDim) * Opts.Min_var, [1, 1, K]);

%% 3. Cálculo de BIC Ponderado
p = K * (nbDim * (nbDim + 1) / 2 + nbDim + 1) - 1;
BIC = -2 * sum(W .* log(F)) + p * log(nbData);

% Salida del modelo
GMModel.Mu = Mu;
GMModel.Sigma = Sigma;
GMModel.Priors = Priors;
GMModel.BIC = BIC;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% FUNCIONES LOCALES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function prob = gaussPDF(Data, Mu, Sigma)
    [nbData, nbDim] = size(Data);
    Data_Mu = Data - Mu;
    prob = sum((Data_Mu / Sigma) .* Data_Mu, 2);
    prob = exp(-0.5 * prob) / sqrt((2*pi)^nbDim * (abs(det(Sigma)) + realmin));
end

