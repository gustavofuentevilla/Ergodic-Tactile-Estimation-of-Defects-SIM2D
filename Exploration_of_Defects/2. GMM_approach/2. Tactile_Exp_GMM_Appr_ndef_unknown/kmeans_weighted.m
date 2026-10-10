function [labels, centroids, inertia] = kmeans_weighted(X, k, weights, Opts)
% KMEANS_WEIGHTED Algoritmo K-Means ponderado con inicialización K-Means++.
%
% Entradas:
%   X          - Matriz de datos de tamaño (N x D).
%   k          - Número de clústeres.
%   weights    - Vector de pesos de las muestras de tamaño (N x 1).
%   Opts.replicates - Número de réplicas/ejecuciones independientes (defecto = 1).
%   Opts.init       - Método de inicio: "kmeans++" (defecto) o "random".
%   Opts.tol        - Tolerancia relativa según la norma de Frobenius (defecto = 1e-4).
%   Opts.max_iter   - Número máximo de iteraciones (defecto = 300).
%
% Salidas:
%   labels    - Vector de etiquetas asignadas a cada punto (N x 1).
%   centroids - Centroides finales calculados (k x D).
%   inertia   - Suma ponderada de las distancias al cuadrado más cercanas.

arguments
    X (:,:) {mustBeNumeric, mustBeNonempty}
    k (1,1) {mustBeInteger, mustBePositive}
    weights (:,1) {mustBeNumeric} = []
    Opts.replicates (1,1) {mustBeInteger, mustBePositive} = 1
    Opts.init (1,1) string {mustBeMember(Opts.init, ["kmeans++", "random"])} = "kmeans++"
    Opts.tol (1,1) {mustBeNumeric, mustBePositive} = 1e-4
    Opts.max_iter (1,1) {mustBeInteger, mustBePositive} = 300
end

[N, ~] = size(X);

% Asignar pesos unitarios si el vector viene vacío
if isempty(weights)
    weights = ones(N, 1);
else
    weights = weights(:);
end

best_inertia = Inf;
best_labels = [];
best_centroids = [];

for rep = 1:Opts.replicates
    % 1. Inicialización de centroides (sin utilizar pesos)
    centroids = initialize_centroids(X, k, Opts.init);
    labels = zeros(N, 1);

    for iter = 1:Opts.max_iter
        old_centroids = centroids;

        % Paso de asignación: distancia euclídea al cuadrado al centroide más cercano
        D = pdist2(X, centroids, 'squaredeuclidean');
        [~, labels] = min(D, [], 2);

        % Paso de actualización: media ponderada con los pesos
        for j = 1:k
            mask = (labels == j);
            w_cluster = weights(mask);
            sum_w = sum(w_cluster);

            if sum_w > 0
                centroids(j, :) = sum(X(mask, :) .* w_cluster, 1) / sum_w;
            else
                % Re-inicializar aleatoriamente si un clúster queda vacío
                centroids(j, :) = X(randi(N), :);
            end
        end

        % Criterio de convergencia: Tolerancia relativa basada en la norma de Frobenius
        % del cambio de los centroides respecto a la iteración previa
        frob_shift = norm(centroids - old_centroids, 'fro');
        frob_old = norm(old_centroids, 'fro');
        
        if (frob_shift / (frob_old + eps)) < Opts.tol
            break;
        end
    end

    % Cálculo de la inercia ponderada para seleccionar la mejor réplica
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
inertia = best_inertia;

end

%% Función auxiliar para la inicialización de centroides
function centroids = initialize_centroids(X, k, method)
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
    end
end