function [labels, centroids] = kmeans_weighted(X, k, weights, Opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted K-means algorithm: the centroids update takes into account the
% weights specified, initializations doesn't.
%
% Inputs:
%   X                   (N x D)    Data with N observations and D 
%                                  dimensions
%   k                   (1 x 1)    Number of clusters
%   weights             (N x 1)    Weights vector
%
% Options (Opts):
%   replicates          (1 x 1)    Number of replicates
%   init                (string)   Initialization methods: "kmeans++",
%                                  "random" or "density_based".
%   tol                 (1 x 1)    Relative tolerance according to 
%                                  Frobenius norm
%   max_iter            (1 x 1)    Maximum number of iterations
%   percent             (1 x 1)    Percent covering parameter for
%                                  "density_based" initialization method
%
% Outputs:
%   labels              (N x 1)    Cluster Labels asigned to data points.
%   centroids           (k x D)    Cluster centroids
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    X (:,:) {mustBeNumeric, mustBeNonempty}
    k (1,1) {mustBeInteger, mustBePositive}
    weights (:,1) {mustBeNumeric} = []
    % Optional inputs
    Opts.replicates (1,1) {mustBeInteger, mustBePositive} = 1
    Opts.init (1,1) string {mustBeMember(Opts.init, ["kmeans++",...
                            "random", "density_based"])} = "kmeans++"
    Opts.tol (1,1) {mustBeNumeric, mustBePositive} = 1e-4
    Opts.max_iter (1,1) {mustBeInteger, mustBePositive} = 300
    Opts.percent (1,1) {mustBeNumeric, mustBePositive} = 0.02
end

[N, ~] = size(X);

% Assign unitary weights if empty
if isempty(weights)
    weights = ones(N, 1);
else
    weights = weights(:);
end

best_inertia = Inf;
best_labels = [];
best_centroids = [];

for rep = 1:Opts.replicates
    % Centroids initialization (unweighted)
    centroids = initialize_centroids(X, k, Opts.init, Opts.percent);

    for iter = 1:Opts.max_iter
        old_centroids = centroids;

        % Data assignment to the closest centroid
        D = pdist2(X, centroids, 'squaredeuclidean');
        [~, labels] = min(D, [], 2);

        % Centroids weighted update
        for j = 1:k
            mask = (labels == j);
            w_cluster = weights(mask);
            sum_w = sum(w_cluster);

            if sum_w > 0
                centroids(j, :) = sum(X(mask, :) .* w_cluster, 1) / sum_w;
            else
                % Random re-initialization if cluster is empty
                centroids(j, :) = X(randi(N), :);
            end
        end

        % Convergence check (frobenius norm)
        frob_shift = norm(centroids - old_centroids, 'fro');
        frob_old = norm(old_centroids, 'fro');
        
        if (frob_shift / (frob_old + eps)) < Opts.tol
            break;
        end
    end

    % Weighted inertia computation to select the best replicate
    D = pdist2(X, centroids, 'squaredeuclidean');
    [min_dist_sq, labels] = min(D, [], 2);
    inertia = sum(weights .* min_dist_sq);

    if inertia < best_inertia
        best_inertia = inertia;
        best_labels = labels;
        best_centroids = centroids;
    end
end

labels = best_labels;
centroids = best_centroids;
% inertia = best_inertia;

end

%% Función auxiliar para la inicialización de centroides
function centroids = initialize_centroids(X, k, method, percent)
    [N, D] = size(X);
    centroids = zeros(k, D);

    if method == "random"
        centroids = X(randperm(N, k), :);

    elseif method == "kmeans++"
        % 1. Primer centroide elegido aleatoriamente con distribución uniforme
        first_idx = randi(N);
        centroids(1, :) = X(first_idx, :);

        % 2. Selección probabilística de los k-1 centroides restantes
        for m = 2:k
            % Distancia mínima al cuadrado hacia los centroides ya existentes
            D2 = min(pdist2(X, centroids(1:m-1, :), 'squaredeuclidean'), [], 2);

            sum_D2 = sum(D2);
            if sum_D2 == 0
                prob = ones(N, 1) / N;
            else
                prob = D2 / sum_D2;
            end

            % Muestreo proporcional a D(x)^2 (sin considerar pesos)
            cumprob = cumsum(prob);
            r = rand();
            next_idx = find(cumprob >= r, 1, 'first');

            if isempty(next_idx)
                next_idx = randi(N);
            end
            centroids(m, :) = X(next_idx, :);
        end

    elseif method == "density_based"
        % Matriz de distancias euclídeas de todos los pares
        D_matrix = pdist2(X, X);

        % 1. Radio del vecindario R a partir del percentil
        triu_mask = triu(true(N), 1);
        dist_vec = sort(D_matrix(triu_mask));
        num_dist = length(dist_vec);
        pos = max(1, round(num_dist * percent));
        R = dist_vec(pos);
        if R == 0, R = eps; end

        % 2. Cálculo de densidad de cada objeto
        density = sum(exp(-(D_matrix / R).^2), 2);
        avg_density = mean(density);

        % 3. Cálculo del MDD (Minimum Density Distance)
        MDD = zeros(N, 1);
        for i = 1:N
            higher_dens_idx = find(density > density(i));
            if isempty(higher_dens_idx)
                % Punto de máxima densidad: distancia máxima al objeto más lejano
                MDD(i) = max(D_matrix(i, :));
            else
                % Distancia al punto más cercano con mayor densidad
                MDD(i) = min(D_matrix(i, higher_dens_idx));
            end
        end

        % 4. Filtrado por densidad superior a la media y ordenación por MDD descendente
        candidate_mask = (density > avg_density);
        candidate_indices = find(candidate_mask);

        if length(candidate_indices) >= k
            [~, sort_order] = sort(MDD(candidate_indices), 'descend');
            selected_indices = candidate_indices(sort_order(1:k));
        else
            % Resguardo si la cantidad de candidatos filtra demasiados puntos
            [~, sort_order] = sort(MDD, 'descend');
            selected_indices = sort_order(1:k);
        end

        centroids = X(selected_indices, :);
    end

end