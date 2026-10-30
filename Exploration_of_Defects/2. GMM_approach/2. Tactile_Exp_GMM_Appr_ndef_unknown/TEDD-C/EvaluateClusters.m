function optimal_clust = EvaluateClusters(Data, W, k_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the optimal number of clusters based on silhouettes criterion
%
% Inputs:
%   Data            (N x 2)     2D Data points
%   W               (N x 1)     Weights vector
%   k_max           (1 x 1)     Maximum number of clusters to evaluate
%
% Outputs:
%   optimal_clust   (1 x 1)     Optimal number of clusters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

score = NaN(k_max, 1);

for k = 2:k_max

    [L,~] = kmeans_weighted(Data, k, W, "replicates", 3);

    % Alternative: Fast aproximation of weighted silhouette criterion
    %
    % s = silhouette(Data, L);
    % 
    % score(k) = mean(s, "Weights", W);

    % Exact weighted silhouette criterion
    score(k) = silhouettes_criterion(Data, L, "Weights", W);

end

% Select the number of clusters with the highest silhouettes score
[sil_val, idx] = max(score);

optimal_clust = idx;



end