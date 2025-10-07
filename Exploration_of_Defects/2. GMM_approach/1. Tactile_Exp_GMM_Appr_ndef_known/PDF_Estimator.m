function [Phi_hat_x_next, Estim_sol] = PDF_Estimator(X_e, V, Par_PDF)

Omega = Par_PDF.Omega;
K = Par_PDF.K;
thres_meas = Par_PDF.thres_meas;
iteration = Par_PDF.iteration;
% Feedbacks
Prev_Data = Par_PDF.Prev_Data;
Prev_Priors = Par_PDF.Prev_Priors;
Prev_Mu = Par_PDF.Prev_Mu;
Prev_Sigma = Par_PDF.Prev_Sigma;
Prev_Sigma_a = Par_PDF.Prev_Sigma_a;
Prev_Mu_found = Par_PDF.Prev_Mu_found;
Prev_Sigma_found = Par_PDF.Prev_Sigma_found;
% Parameters
MinVariation = Par_PDF.MinVariation;
MinAxisLengths = Par_PDF.MinAxisLengths;
DataEscFact = Par_PDF.DataEscFact;
% Threshold about total variation for finding a defect
Thres_PostVariation = Par_PDF.Thres_Variation; %13 cm

Prev_Mu_found(1,:) = []; %Remove initial element (which was stablished outside only to define dimensions)
Prev_Sigma_found(:,:,1) = [];

%% %%%%%%%%%%%%%%%% Preprocesing data %%%%%%%%%%%%%%%%%%%%%%%%%
idx_V = V > thres_meas; %Locate values above threshold
Preprocessed_V = V;
Preprocessed_V(idx_V) = V(idx_V)*DataEscFact; %Scale those values to give it more importance %15
Preprocessed_V(~idx_V) = V(~idx_V)*0; %and less importance to the values under the threshold 

Data_current = [X_e, Preprocessed_V];

% adding the previous data (if any)
Data = [Prev_Data; 
        Data_current];

% Removing data related to already found defects
if ~isempty(Prev_Mu_found)        %If some defect has been found
    % nbDef_found = size(Prev_Mu_found, 1);
    % Prev_Sigma_found = Prev_Sigma_found + diag([0.0005, 0.0005]); % adding some variance
    % sd_found = zeros(size(Prev_Sigma_found));
    % Sigma_ast_found = zeros(size(Prev_Sigma_found));
    % isInElipse = false([height(Data), nbDef_found]);
    % 
    % for j = 1:nbDef_found
    %     % Standard deviation
    %     sd_found(:,:,j) = sqrtm(Prev_Sigma_found(:,:,j)); 
    %     % 3*Standard deviation
    %     Sigma_ast_found(:,:,j) = 3*sd_found(:,:,j);
    %     % Finding rotation matrix "V" and principal axis half-lenghts
    %     % "lambdas = [l_1; ...; l_max]" (sorted)
    %     [V, D] = eig(Sigma_ast_found(:,:,j));
    %     [lambdas_s, ind] = sort(diag(D));
    %     % D_sorted = D(ind, ind);
    %     V_sorted = V(:,ind);
    %     % Coordinates transformation
    %     P_ast = V_sorted'*( Data(:,1:2) - Prev_Mu_found(j,:) )';
    %     % Elipse equation
    %     Elips_eq = sum( (P_ast.^2) ./ (lambdas_s.^2), 1 )';
    %     % Indice de los datos que están dentro de las elipses
    %     isInElipse(:,j) = Elips_eq <= 1; 
    % end
    % 
    % %índice de todos los datos a eliminar (los que están dentro de las elipses)
    % idx_erase = any(isInElipse, 2); 
    % % Eliminación de todos los datos dentro de las elipses 
    % (incluso de los datos recopilados en iteraciones anteriores)
    
    % Remove data points that are inside the already found defects
    idx_remData = IsDataInEllipse(Data, Prev_Mu_found, Prev_Sigma_found);
    Data(idx_remData,:) = [];

end

V_int = round(Data(:,3));    %V Conversion to int

Data_Xe_hist_V = repelem(Data(:,1:2), V_int, 1); %Repeat elements on spatial domain (trajectory points)

%% %%%%%%%%%%%% Gaussian Mixture Model of Preprocessed data %%%%%%%%%%%%%%%

% Solution 1 for unknow number of clusters, but a range given, compute:
% cev = evalclusters(statsNorm,"gmdistribution","silhouette",...
%     KList=2:4);
% bestK = cev.OptimalK;

% GMM with incremental EM algorithm
[Priors, Mu, Sigma] = GMM_EM(Data_Xe_hist_V, K, Prev_Priors, Prev_Mu, Prev_Sigma, iteration, ...
                            "Max_iter", 500, "Min_var", 1e-8);

% Solution 2 for unknown number of clusters, but a range given, computing
% AIC criterion
% AIC = zeros(1,4);
% GMModels = cell(1,4);
% gmm_opt = statset('MaxIter',2000);
% for k = 1:4
%     GMModels{k} = fitgmdist(X,k,'Options',gmm_opt,'CovarianceType','diagonal',...
%                   'SharedCovariance', true, "Replicates", 5);
%     AIC(k)= GMModels{k}.AIC;
% end
% 
% [minAIC,numComponents] = min(AIC);
% BestModel = GMModels{numComponents}

%% Minimum variance and Minimum axis lengths constraints

% Definition of minimum lenght sum of elipse axis --> "MinVariation"
% Definition of minimum axis lengths --> "MinAxesLength"

% Augmented covariance matrix (the one we want to compute)
Sigma_a = zeros(size(Sigma));
flag_readjust = false(1, K);
% Computation of actual axis lengths
stdev = zeros(size(Sigma));
Sigma_ast = zeros(size(Sigma));
Variation_Sigma = zeros(1, K);

for j = 1:K
    % Standard deviation
    stdev(:,:,j) = sqrtm(Sigma(:,:,j)); 
    % 3*Standard deviation that represents 99% of data
    Sigma_ast(:,:,j) = 3*stdev(:,:,j); 
    % Eigenvectors = V, Eigenvalues diagonal matriz = D
    [V, D] = eig(Sigma_ast(:,:,j));
    % Sorting with Max eigenvalue at last spot
    [r_j, ind] = sort(diag(D)); %r_j is the vector of elipse radius of Sigma on the plane
    D_sorted = D(ind, ind); 
    V_sorted = V(:,ind);
    % Compute Ratio of original Sigma
    % Ratio = r_j(2)/r_j(1);
    Ratio = 1;
    % Variación total = trace(Sigma_ast_phi) = sum(r_j)
    Variation_Sigma(:,j) = sum(r_j);

    % Total Variation Constraint:
    % If the sum of axis lengths of Sigma is less than the 
    % limit (MinVariation), then re-adjust to achieve the minimal

    if (Variation_Sigma(:,j) < MinVariation) %|| (min(r_elips(:,j)) < MinAxesLength)
        flag_readjust(j) = true;
        % Compute the variation needed to achieve the minimal
        DeltaVariation = MinVariation - Variation_Sigma(:,j);
        % Compute Extension to be made over axes mantaining the aspect 
        % ratio (r_2_bar = Ratio*r_1_bar, r_1_bar < r_2_bar)
        r_bar = zeros(size(r_j));
        % Compute axis lengths offsets to achieve the MinVariation
        r_bar(1) = DeltaVariation/(1 + Ratio); 
        r_bar(2) = Ratio*r_bar(1);
        % New diagonal eigenvalues matriz with extended radius
        D_a = D_sorted + diag(r_bar);
        % Augmented covariance matrix
        Sd_a = V_sorted*D_a*V_sorted' / 3;
        Sigma_a(:,:,j) = Sd_a * Sd_a;
        % Update lengths (V_sorted is the same as original Sigma)
        D_sorted = D_a; 
        r_j = diag(D_a);
    end

    % Minimum axis lengths constraint:

    if r_j(1) < MinAxisLengths
        flag_readjust(j) = true;
        % Compute Extension to be made over axes mantaining the aspect 
        % ratio (r_2_bar = Ratio*r_1_bar, r_1_bar < r_2_bar)
        r_bar = zeros(size(r_j));
        % Compute axis lengths offsets to achieve the MinAxisLength
        r_bar(1) = MinAxisLengths - r_j(1); 
        r_bar(2) = Ratio*r_bar(1);
        % New diagonal eigenvalues matriz with extended radius
        D_a = D_sorted + diag(r_bar);
        % Augmented covariance matrix
        Sd_a = V_sorted*D_a*V_sorted' / 3;
        Sigma_a(:,:,j) = Sd_a * Sd_a;
    end
    
    % If there was no re-adjustment, the updated Sigma_a is the original
    % Sigma
    if ~flag_readjust(j)
        Sigma_a(:,:,j) = Sigma(:,:,j);
    end

end

%% Kullback–Leibler divergence (Relative entropy)
% Pos_p = zeros(height(Omega), n_def);
% for i = 1:n_def
%     Pos_p(:,i) = mvnpdf(Omega, Mu(i,:), Sigma(:,:,i));
% end
% 
% if iteration == 1
% 
%     Pre_p = Prev_Phi_hat_x;
%     D_KL = zeros(1, n_def);
%     for i = 1:n_def
%         D_KL(:,i) = sum( Pos_p(:,i).*log( Pos_p(:,i) ./ Pre_p ) );
%     end
% 
% else
% 
%     D_KL = zeros(1, n_def);
%     n_dim = size(Omega, 2);
%     for i = 1:n_def
%         D_KL(:,i) = (1/2)*( trace(Prev_Sigma(:,:,i)\Sigma(:,:,i)) - n_dim +...
%             (Prev_Mu(i,:) - Mu(i,:))*(Prev_Sigma(:,:,i)\(Prev_Mu(i,:) - Mu(i,:))') + ...
%             log(det(Prev_Sigma(:,:,i)) / det(Sigma(:,:,i))) );
%     end
% 
% end

%% Metrics

PriorSigma_SumDist = zeros(K, 1);
PostSigma_SumDist = zeros(K, 1);
for i = 1:K
    stdev_prev = sqrtm(Prev_Sigma_a(:,:,i));
    r_elips_Prev = eig(3*stdev_prev);
    PriorSigma_SumDist(i) = sum(r_elips_Prev); %Total variation for previous Sigma

    stdev_post = sqrtm(Sigma_a(:,:,i));
    r_elips_post = eig(3*stdev_post);
    PostSigma_SumDist(i) = sum(r_elips_post); %Total variation for posterior Sigma
end
% Change distance between Previous Mu and Posterior Mu 
Delta_SumDist = abs(PriorSigma_SumDist - PostSigma_SumDist);
Delta_Mu = sqrt(sum((Mu - Prev_Mu).^2, 2));

%% Finding defects and defining the next PDF

Thres_DeltaMu = 0.001; % 1 mm
Thres_DeltaVariation = 0.001; %1 mm
% Thres_PostVariation = Par_PDF.Thres_Variation; %13 cm

cond_DeltaMu = (Delta_Mu <= Thres_DeltaMu);
cond_DeltaVariation = (Delta_SumDist <= Thres_DeltaVariation);
cond_PostVariation = (PostSigma_SumDist <= Thres_PostVariation);

Cond = cond_PostVariation; % | cond_DeltaMu | cond_DeltaSumDist  ;

Def_found.Mu = [];
Def_found.Sigma = [];
Def_found.Priors = [];

% If true --> One or more Defects have been found
if any(Cond) 
    % Index of defects found ( idx_def = find(Cond) )
    idx_def = Cond;  
    % Saving Gaussian parameters related
    Def_found.Mu = Mu(idx_def,:);
    Def_found.Sigma = Sigma_a(:,:,idx_def);
    Def_found.Priors = Priors(idx_def);
    % Removing defect parameters from GMM solution
    Mu(idx_def,:) = [];
    Sigma(:,:,idx_def) = [];
    Sigma_a(:,:,idx_def) = [];
    Priors(:,idx_def) = [];
end

% Check if we've done
flag_done = isempty(Mu);

% If we've done, export a zero PDF 
if flag_done
    Phi_hat_x_next = zeros(height(Omega));
    Estim_sol.GMModel = [];
else % IF we've not done, compute and export next GMM
    % Next PDF
    GMModel = gmdistribution(Mu, Sigma_a, Priors); %This function normalizes Priors if these doesn't sum to 1.
    Phi_hat_x_next = pdf(GMModel, Omega);
    % Export
    Estim_sol.GMModel = GMModel;
end

%% Output structure
Estim_sol.Data = Data;
Estim_sol.Preprocessed_V = Preprocessed_V;
Estim_sol.Data_Xe_hist_V = Data_Xe_hist_V;
Estim_sol.Priors = Priors;
Estim_sol.Mu = Mu;
Estim_sol.Sigma = Sigma;
Estim_sol.Sigma_a = Sigma_a;
% Estim_sol.D_KL = D_KL;
Estim_sol.Priors_found = Def_found.Priors;
Estim_sol.Mu_found = Def_found.Mu;
Estim_sol.Sigma_found = Def_found.Sigma;
Estim_sol.Delta_Mu = Delta_Mu;
Estim_sol.Delta_SumDist = Delta_SumDist;
Estim_sol.PriorSigma_SumDist = PriorSigma_SumDist;
Estim_sol.PostSigma_SumDist = PostSigma_SumDist;
Estim_sol.MinVariation = MinVariation;
Estim_sol.cond_DeltaMu = cond_DeltaMu;
Estim_sol.cond_DeltaVariation = cond_DeltaVariation;
Estim_sol.cond_PostVariation = cond_PostVariation;
Estim_sol.flag_readjust = flag_readjust;
Estim_sol.flag_done = flag_done;

end

function idx_InEllipse = IsDataInEllipse(Data, Mu, Sigma)

% Evaluates and indentifies the Data points that are inside the 2D 
% ellipses defined by (Mu, Sigma)

nbEllipses = size(Mu, 1);
Sigma = Sigma + diag([0.0005, 0.0005]); % adding a little offset
S_D = zeros(size(Sigma));
Sigma_ast = zeros(size(Sigma));
isInEllipse = false([height(Data), nbEllipses]);

for j = 1:nbEllipses
    % Standard deviation
    S_D(:,:,j) = sqrtm(Sigma(:,:,j)); 
    % 3*Standard deviation
    Sigma_ast(:,:,j) = 3*S_D(:,:,j);
    % Finding rotation matrix "V" and principal axis half-lenghts
    % "lambdas = [l_1; ...; l_max]" (sorted)
    [V, D] = eig(Sigma_ast(:,:,j));
    [lambdas_s, ind] = sort(diag(D));
    % D_sorted = D(ind, ind);
    V_sorted = V(:,ind);
    % Coordinates transformation
    P_ast = V_sorted'*( Data(:,1:2) - Mu(j,:) )';
    % Ellipse equation
    Ellipse_eq = sum( (P_ast.^2) ./ (lambdas_s.^2), 1 )';
    % Indice de los datos que están dentro de las elipses
    isInEllipse(:,j) = Ellipse_eq <= 1; 
end

%índice de todos los datos a eliminar (los que están dentro de las elipses)
idx_InEllipse = any(isInEllipse, 2); 

% Eliminación de todos los datos dentro de las elipses
% Data(idx_erase,:) = [];

end