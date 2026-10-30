function [S, s_i] = silhouettes_criterion(Data, Labels, Opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the global and independent continuos weighted Silhouette scores
%
% Inputs:
%   Data          (N x D)   Data matriz with N observations and D 
%                           dimensions.
%   Labels        (N x 1)   Cluster labels of data points.
%   Opts.Weights  (N x 1)   Positive weights vector.
%   Opts.Distance (string)  Distance metric (default: "squaredeuclidean")
%
% Outputs:
%   S             (1 x 1)   Global weighted Silhouette score.
%   s_i           (N x 1)   Indvidual Silhouettes score.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    Data (:,:) double {mustBeReal, mustBeNonempty}
    Labels (:,1) {mustBeNonempty}
    % Optional inputs
    Opts.Weights (:,1) double {mustBePositive} = ones(size(Data, 1), 1)
    Opts.Distance (1,1) string {mustBeMember(Opts.Distance,...
        ["euclidean", "squaredeuclidean", "seuclidean",...
        "fasteuclidean", "fastsquaredeuclidean", "fastseuclidean",...
        "mahalanobis", "cityblock", "minkowski", "chebychev",...
        "cosine", "correlation", "hamming", "jaccard", "spearman",...
        ])} = "squaredeuclidean"
end

N = size(Data, 1);

% Consistent lenghts validation
assert(size(Labels, 1) == N,...
    "Labels column vector must have the same number of elements" + ...
    "as rows in Data.");
assert(size(Opts.Weights, 1) == N,...
    "Weights column vector must have the same number of elements" + ...
    "as rows in Data.");

Weights = Opts.Weights;

% Get unique cluster identifiers
[unique_labels, ~, label_idx] = unique(Labels);
K = length(unique_labels);

% Initialize intra cluster distance (a) and inter cluster distance (b)
a = zeros(N, 1);
b = inf(N, 1);

% Define Silhouette as zero if there is only one cluster
if K == 1
    s_i = zeros(N, 1);
    S = 0;
    return;
end

% Pre-compute the total weights sum of every cluster
cluster_weights = zeros(K, 1);
for k = 1:K
    cluster_weights(k) = sum(Weights(label_idx == k));
end

% Compute a(x_i) and b(x_i)
for k = 1:K
    % Masks for data inside and outside the current cluster k
    in_k_mask = (label_idx == k);
    out_k_mask = ~in_k_mask;
    
    % Extract data and weight inside cluster k
    pts_k = Data(in_k_mask, :);
    w_k = Weights(in_k_mask);
    
    if isempty(pts_k)
        continue;
    end
    
    % Compute distances and weighted sum
    D_to_k = pdist2(Data, pts_k, Opts.Distance);
    sum_wd = D_to_k * w_k; 
    
    % Intra-cluster distance a(x_i)
    W_rem = cluster_weights(k) - w_k; 
    in_k_idx = find(in_k_mask);
    
    valid_a = W_rem > 0; 
    a(in_k_idx(valid_a)) = sum_wd(in_k_idx(valid_a)) ./ W_rem(valid_a);
    
    % Inter-cluster distance b(x_i)
    if cluster_weights(k) > 0
        delta_to_k = sum_wd(out_k_mask) / cluster_weights(k);
        b(out_k_mask) = min(b(out_k_mask), delta_to_k);
    end
end

% Individual Silhouettes score s_i
max_ab = max(a, b);
s_i = zeros(N, 1);

valid_s = max_ab > 0;
s_i(valid_s) = (b(valid_s) - a(valid_s)) ./ max_ab(valid_s);

% Weighted global Silhouette score S
S = sum(Weights .* s_i) / sum(Weights);

end