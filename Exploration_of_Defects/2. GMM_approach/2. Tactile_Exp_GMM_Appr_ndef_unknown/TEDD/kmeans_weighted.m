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

%% Initialization methods
function centroids = initialize_centroids(X, k, method, percent)
    [N, D] = size(X);
    centroids = zeros(k, D);

    if method == "random"
        centroids = X(randperm(N, k), :);

    elseif method == "kmeans++"
        % Choosing the first centroid randomly
        first_idx = randi(N);
        centroids(1, :) = X(first_idx, :);

        % Choosing the rest of centroids
        for m = 2:k
            % Minimum squared distance to the existing centroids
            D2 = min(pdist2(X, centroids(1:m-1, :), 'squaredeuclidean'), [], 2);

            sum_D2 = sum(D2);
            if sum_D2 == 0
                prob = ones(N, 1) / N;
            else
                prob = D2 / sum_D2;
            end

            % Proportional sampling
            cumprob = cumsum(prob);
            r = rand();
            next_idx = find(cumprob >= r, 1, 'first');

            if isempty(next_idx)
                next_idx = randi(N);
            end
            centroids(m, :) = X(next_idx, :);
        end

    elseif method == "density_based"
        % Distance matrix
        D_matrix = pdist2(X, X);

        % Computing radius R
        triu_mask = triu(true(N), 1);
        dist_vec = sort(D_matrix(triu_mask));
        num_dist = length(dist_vec);
        pos = max(1, round(num_dist * percent));
        R = dist_vec(pos);
        if R == 0, R = eps; end

        % Density computation of data points
        density = sum(exp(-(D_matrix / R).^2), 2);
        avg_density = mean(density);

        % Minimum Density Distance (MDD) computation
        MDD = zeros(N, 1);
        for i = 1:N
            higher_dens_idx = find(density > density(i));
            if isempty(higher_dens_idx)
                % Maximum density point: maximum distance to the furthest point
                MDD(i) = max(D_matrix(i, :));
            else
                % Distance to the closest point with higher density
                MDD(i) = min(D_matrix(i, higher_dens_idx));
            end
        end

        % Density filtering by the average
        candidate_mask = (density > avg_density);
        candidate_indices = find(candidate_mask);

        % Select the centroids
        if length(candidate_indices) >= k
            [~, sort_order] = sort(MDD(candidate_indices), 'descend');
            selected_indices = candidate_indices(sort_order(1:k));
        else
            [~, sort_order] = sort(MDD, 'descend');
            selected_indices = sort_order(1:k);
        end

        centroids = X(selected_indices, :);
    end

end