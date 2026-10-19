function optimal_clust = EvaluateClusters(Data, W, k_max)
% Data  (N x 2)     Data points
% W     (N x 1)     Weights
% k_max (1 x 1)     Maximum number of clusters to evaluate

score = NaN(k_max, 1);

for k = 2:k_max

    [L,~,~] = kmeans_weighted(Data, k, W, "replicates", 3);

    % Aproximation of weighted silhouette criterion
    %
    % s = silhouette(Data, L);
    % 
    % score(k) = mean(s, "Weights", W);

    % Exact weighted silhouette criterion
    score(k) = silhouettes_criterion(Data, L, "Weights", W);

end

[sil_val, idx] = max(score);

optimal_clust = idx;



end