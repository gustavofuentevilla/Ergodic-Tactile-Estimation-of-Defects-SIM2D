function [S, s_i] = silhouettes_criterion(Data, Labels, Opts)
    % SILHOUETTES_CRITERION Calcula el coeficiente de Silhouette Ponderado Continuo.
    %
    % Entradas:
    %   Data          - Matriz (N x D) de características.
    %   Labels        - Vector columna (N x 1) con las etiquetas de los clusters.
    %   Opts.Weights  - (Opcional) Vector columna (N x 1) de pesos continuos 
    %                   positivos (w_i > 0). Por defecto: vector de unos.
    %   Opts.Distance - (Opcional) Métrica de distancia. Por defecto: "sqeuclidean".
    %
    % Salidas:
    %   S       - Escalar con el valor global de Silhouette ponderado.
    %   s_i     - Vector columna (N x 1) con los valores individuales S(x_i).

    arguments
        % Data debe ser una matriz numérica real y no vacía
        Data (:,:) double {mustBeReal, mustBeNonempty}
        
        % Labels debe ser un vector columna
        Labels (:,1) {mustBeNonempty}
        
        % Opts.Weights debe ser un vector columna de números positivos. 
        % Si no se especifica, toma un vector de unos del tamaño de Data.
        Opts.Weights (:,1) double {mustBePositive} = ones(size(Data, 1), 1)
        
        % Métrica de distancia por defecto
        Opts.Distance (1,1) string {mustBeMember(Opts.Distance,...
            ["euclidean", "squaredeuclidean", "seuclidean",...
            "fasteuclidean", "fastsquaredeuclidean", "fastseuclidean",...
            "mahalanobis", "cityblock", "minkowski", "chebychev",...
            "cosine", "correlation", "hamming", "jaccard", "spearman",...
            ])} = "squaredeuclidean"
    end
    
    N = size(Data, 1);
    
    % Validación de longitudes consistentes
    % (Se realiza aquí porque el bloque arguments nativo no permite vincular 
    % explícitamente tamaños entre distintas variables de forma directa)
    assert(size(Labels, 1) == N, 'Labels debe tener el mismo número de filas que Data.');
    assert(size(Opts.Weights, 1) == N, 'Opts.Weights debe tener el mismo número de filas que Data.');
    
    Weights = Opts.Weights;
    
    % Obtener identificadores únicos de los clusters
    [unique_labels, ~, label_idx] = unique(Labels);
    K = length(unique_labels);
    
    % Inicializar distancias intra (a) e inter (b)
    a = zeros(N, 1);
    b = inf(N, 1);
    
    % Si hay un solo cluster, el Silhouette se define como 0
    if K == 1
        s_i = zeros(N, 1);
        S = 0;
        return;
    end
    
    % Precalcular los pesos totales de cada cluster (W_k)
    cluster_weights = zeros(K, 1);
    for k = 1:K
        cluster_weights(k) = sum(Weights(label_idx == k));
    end
    
    % Calcular a(x_i) y b(x_i) iterando por cada cluster destino
    for k = 1:K
        % Máscaras lógicas para puntos dentro y fuera del cluster k
        in_k_mask = (label_idx == k);
        out_k_mask = ~in_k_mask;
        
        % Extraer datos y pesos específicos del cluster k
        pts_k = Data(in_k_mask, :);
        w_k = Weights(in_k_mask);
        
        if isempty(pts_k)
            continue;
        end
        
        % Calcular distancias (N x N_k) y suma ponderada (Vectorización)
        D_to_k = pdist2(Data, pts_k, Opts.Distance);
        sum_wd = D_to_k * w_k; 
        
        % --- 1. Distancia intra-cluster a(x_i) ---
        W_rem = cluster_weights(k) - w_k; 
        in_k_idx = find(in_k_mask);
        
        valid_a = W_rem > 0; 
        a(in_k_idx(valid_a)) = sum_wd(in_k_idx(valid_a)) ./ W_rem(valid_a);
        
        % --- 2. Distancia inter-cluster b(x_i) ---
        if cluster_weights(k) > 0
            delta_to_k = sum_wd(out_k_mask) / cluster_weights(k);
            b(out_k_mask) = min(b(out_k_mask), delta_to_k);
        end
    end
    
    % Calcular Silhouette individual s_i
    max_ab = max(a, b);
    s_i = zeros(N, 1);
    
    valid_s = max_ab > 0;
    s_i(valid_s) = (b(valid_s) - a(valid_s)) ./ max_ab(valid_s);
    
    % Calcular Silhouette global ponderado S
    S = sum(Weights .* s_i) / sum(Weights);
end