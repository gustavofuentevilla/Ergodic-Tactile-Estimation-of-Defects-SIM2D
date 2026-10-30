function [Phi_hat_x_next, Estim_sol] = PDF_Estimator(X_e, V, Par_PDF)

Omega = Par_PDF.Omega;
dx_1 = Par_PDF.dx_1;
dx_2 = Par_PDF.dx_2;
D_max = Par_PDF.D_max;
thres_meas = Par_PDF.thres_meas;
iteration = Par_PDF.iteration;

% Feedbacks
Prev_Data = Par_PDF.Prev_Data;
Prev_Mu_found = Par_PDF.Prev_Mu_found;
Prev_Sigma_found = Par_PDF.Prev_Sigma_found;
Prev_numComponents = Par_PDF.Prev_numComponents;
Prev_Phi_hat_x = Par_PDF.Prev_Phi_hat_x;

% Parameters
MinAxisLengths = Par_PDF.MinAxisLengths;
MaxVarCons = Par_PDF.MaxVarCons;
D_KL_bar_u = Par_PDF.D_KL_bar_u;
eps = Par_PDF.eps;
% Threshold for finding a defect
Thres_PostVariation = Par_PDF.Thres_Variation;
% Limit to consider one single cluster
OneClustDistLimit = Par_PDF.OneClustDistLimit;
flag_ExplorationStage = Par_PDF.flag_ExplorationStage;

% Remove initial element 
% (which was established outside only to define dimensions)
Prev_Mu_found(1,:) = []; 
Prev_Sigma_found(:,:,1) = [];

%% %%%%%%%%%%%%%%%% Preprocessing data %%%%%%%%%%%%%%%%%%%%%%%%%
% Locate sensor measurements above threshold
idx_V = V > thres_meas;

% Current iteration data (Position, Measurement)
Data_current = [X_e(idx_V,:), V(idx_V)];

% Set of all historical data
Data = [Prev_Data; 
        Data_current];

% Removing data related to already found defects 
if ~isempty(Prev_Mu_found)        
    % If some defect has been found
    % Remove data points inside the estimated ellipses
    idx_remData = IsDataInEllipse(Data(:,1:2), Prev_Mu_found, Prev_Sigma_found);
    Data(idx_remData,:) = [];
end

%% Check for No Data Case
% If there is no data above the threshold, return the same distribution for
% double check and terminate the execution of this function
flag_NoData = 0;
if isempty(Data_current)
    Phi_hat_x_next = Prev_Phi_hat_x;
    flag_NoData = 1;
    % Output structure
    Estim_sol.Data = Data;
    Estim_sol.Priors = [];
    Estim_sol.Mu = [];
    Estim_sol.Sigma = [];
    Estim_sol.Sigma_a = [];
    Estim_sol.numComponents = 0;
    Estim_sol.D_KL = [];
    Estim_sol.Priors_found = [];
    Estim_sol.Mu_found = [];
    Estim_sol.Sigma_found = [];
    Estim_sol.PostSigma_SumDist = [];
    Estim_sol.MinVariation = [];
    Estim_sol.cond_PostVariation = [];
    Estim_sol.flag_readjust = [];
    Estim_sol.flag_done = [];
    Estim_sol.flag_NoData = flag_NoData;
    Estim_sol.flag_OneCluster = false;
    Estim_sol.flag_ExplorationStage = true;
    return;
end

%% %%%%%%%%%%%% Gaussian Mixture Model of Preprocessed data %%%%%%%%%%%%%%%

% %%%%%%%%%%%%% Logic for one single cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check distance among X_e Points in the plane and test:
% If the maximum distance of the current batch is less than a certain 
% limit, then there is only one cluster (defect).
% *flag_OneCluster is used to avoid the computation of more than one
% clusters logic below

flag_OneCluster = false;

cHull = convhull(Data_current(:,1:2));
cHull(end) = [];
d_max = max(pdist(Data_current(cHull,1:2)));
if d_max <= OneClustDistLimit
    % One cluster
    numComponents = 1;
    Model = GMM_EM_weighted(Data(:,1:2), Data(:,3), numComponents,...
                   "Max_iter", 500, "Min_var", 1e-10);
    flag_OneCluster = true;
end


% %%%%%%%%%%%%% Logic for more than one cluster %%%%%%%%%%%%%%%%%%%%%%%%%%

if ~flag_OneCluster
    if flag_ExplorationStage
        % EXPLORATION PHASE
        % Define the number of gaussian components to estimate evaluating 
        % the clusters under the distance-based Silhouette criteria 

        numComponents = EvaluateClusters(Data(:,1:2), Data(:,3), D_max);
        
    else
        % EXPLOITATION PHASE
        % Use the Optimal number of components from the last iteration in 
        % the exploration phase minus the number of defects already found 
        % (if any)

        numComponents = Prev_numComponents - size(Prev_Mu_found, 1);

    end
    % Define the model
    Model = GMM_EM_weighted(Data(:,1:2), Data(:,3), numComponents,...
                            "Max_iter", 500, "Min_var", 1e-10);
end

% Extracting model parameters
Mu = Model.Mu;
Sigma = Model.Sigma;
Priors = Model.Priors;

%% Kullback–Leibler divergence
% Calculate KL Divergence from Q (Prior) to P (Posterior) PDFs
Q = Prev_Phi_hat_x + 1e-10;
P = pdf(gmdistribution(Mu, Sigma, Priors), Omega) + 1e-10;
% Riemann integral approximation
D_KL = sum( P .* log(P ./ Q) )*dx_1*dx_2;

%% Compute Constraint as a function of KL divergence

Variation_thres_eps = Thres_PostVariation - eps;

% Define the constraint upper bound based on the first KL measure
if iteration == 1
    D_KL_bar_u = D_KL;
end

% Check for Exploration Stage to compute the constraint
if flag_ExplorationStage
    MinVariation = D_KL*MaxVarCons/D_KL_bar_u;
    if MinVariation <= Variation_thres_eps
        MinVariation = Variation_thres_eps;
    elseif MinVariation >= MaxVarCons
        MinVariation = MaxVarCons;
    end
else
    MinVariation = Variation_thres_eps;
end

% Exploration stage flag set to false if Variation Constraint reaches the
% Threshold for finding a defect
if MinVariation <= Thres_PostVariation
    flag_ExplorationStage = false;
end

%% Minimum variance and Minimum axis lengths constraints

% Definition of minimum length sum of ellipse axis --> "MinVariation"
% Definition of minimum axis lengths --> "MinAxisLengths" (Not used)

% Update Sigma to Sigma_a to match the minimum constraint
[Sigma_a, flag_readjust] = ReadjustConstrainedSigma(Sigma,...
                                                    MinVariation,...
                                                    MinAxisLengths);

%% Metrics

% Compute posterior Sigma axis lengths sum
PostSigma_SumDist = zeros(numComponents, 1);
for i = 1:numComponents

    stdev_post = sqrtm(Sigma_a(:,:,i));
    r_elips_post = eig(3*stdev_post);

    PostSigma_SumDist(i) = sum(r_elips_post); 
end

%% Finding defects and defining the next PDF

% Condition to "find" a defect
cond_PostVariation = (PostSigma_SumDist <= Thres_PostVariation);

Cond = cond_PostVariation;

% Found defects struct
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
flag_done = isempty(Mu); %(comment for no removing case)

% If we've done, export a zero PDF 
if flag_done
    Phi_hat_x_next = zeros(height(Omega), 1);
    Estim_sol.GMModel = [];
else % IF we've not done, compute and export next GMM
    GMModel = gmdistribution(Mu, Sigma_a, Priors); 
    Phi_hat_x_next = pdf(GMModel, Omega);
    % Export
    Estim_sol.GMModel = GMModel;
end

%% Output structure
Estim_sol.Data = Data;
Estim_sol.Priors = Priors;
Estim_sol.Mu = Mu;
Estim_sol.Sigma = Sigma;
Estim_sol.Sigma_a = Sigma_a;
Estim_sol.numComponents = numComponents;
Estim_sol.D_KL = D_KL;
Estim_sol.Priors_found = Def_found.Priors;
Estim_sol.Mu_found = Def_found.Mu;
Estim_sol.Sigma_found = Def_found.Sigma;
Estim_sol.PostSigma_SumDist = PostSigma_SumDist;
Estim_sol.MinVariation = MinVariation;
Estim_sol.cond_PostVariation = cond_PostVariation;
Estim_sol.flag_readjust = flag_readjust;
Estim_sol.flag_done = flag_done;
Estim_sol.flag_NoData = flag_NoData;
Estim_sol.flag_OneCluster = flag_OneCluster;
Estim_sol.flag_ExplorationStage = flag_ExplorationStage;

end

function idx_InEllipse = IsDataInEllipse(Data, Mu, Sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluates and identifies the Data points that are inside the 2D 
% ellipses defined by (Mu, Sigma)
% Inputs:
% Data      N x 2       matrix, N data points, 2 dimensional
% Mu        n x D       n GMM components of dimension D
% Sigma     D x D x n   Covariance matrices
% Outputs:
% idx_InEllipse N x 1   Index vector (true = Data point inside ellipses) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    P_ast = V_sorted'*( Data - Mu(j,:) )';
    % Ellipse equation
    Ellipse_eq = sum( (P_ast.^2) ./ (lambdas_s.^2), 1 )';
    % Index matrix for data points inside the j-th ellipse
    isInEllipse(:,j) = Ellipse_eq <= 1; 
end

% Index vector considering all ellipses
idx_InEllipse = any(isInEllipse, 2); 

end

function [Sigma_a, flag_readjust] = ReadjustConstrainedSigma(Sigma,...
                                    MinVariation, MinAxisLengths)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-adjust covariance Matrix Sigma by checking to the constraint
% MinVariation, additionally one can set a minimum axis lengths value
% Inputs:
% Sigma     D x D x n   Covariance matrices
% MinVariation  1 x 1   Semi-axes sum constraint value
% MinAxisLengths    1 x 1   Optional minimum axis length 
%                           (set to 0 if not used)
% Outputs:
% Sigma_a   D x D x n   Augmented Covariance matrices
% flag_readjust 1 x n   Logic vector indicating for readjusted Sigma
%                       component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of components
nbEllipses = size(Sigma, 3);
% Augmented covariance matrix (the one we want to compute)
Sigma_a = zeros(size(Sigma));
% Flag set to True if that component has been readjusted
flag_readjust = false(1, nbEllipses);

% Computation of current axis lengths
stdev = zeros(size(Sigma));
Sigma_ast = zeros(size(Sigma));
Variation_Sigma = zeros(1, nbEllipses);

for j = 1:nbEllipses
    % Standard deviation
    stdev(:,:,j) = sqrtm(Sigma(:,:,j)); 
    % 3*Standard deviation that represents 99% of data
    Sigma_ast(:,:,j) = 3*stdev(:,:,j); 
    % Eigenvectors = V, Eigenvalues diagonal matrix = D
    [V, D] = eig(Sigma_ast(:,:,j));
    % Sorting with Max eigenvalue at last spot
    % r_j is the vector of ellipse radius of Sigma on the plane
    [r_j, ind] = sort(diag(D)); 
    D_sorted = D(ind, ind); 
    V_sorted = V(:,ind);
    % Compute Ratio of original Sigma if aspect ratio is wanted to be
    % maintained, otherwise set Ratio to 1
    % Ratio = r_j(2)/r_j(1);
    Ratio = 1;
    % Total variation = trace(Sigma_ast) = sum(r_j)
    Variation_Sigma(:,j) = sum(r_j);

    % Total Variation Constraint:
    % If the sum of axis lengths of Sigma is less than the 
    % limit (MinVariation), then re-adjust to achieve the minimal

    if (Variation_Sigma(:,j) < MinVariation)
        flag_readjust(j) = true;
        % Compute the variation needed to achieve the minimal
        DeltaVariation = MinVariation - Variation_Sigma(:,j);
        % Compute Extensions to be made over semi-axes maintaining the 
        % aspect ratio (r_2_bar = Ratio*r_1_bar, r_1_bar < r_2_bar)
        r_bar = zeros(size(r_j));
        % Compute semi-axis lengths offsets to achieve the MinVariation
        r_bar(1) = DeltaVariation/(1 + Ratio); 
        r_bar(2) = Ratio*r_bar(1);
        % New diagonal eigenvalues matrix with extended radius
        D_a = D_sorted + diag(r_bar);
        % Augmented covariance matrix computation
        Sd_a = V_sorted*D_a*V_sorted' / 3;
        Sigma_a(:,:,j) = Sd_a * Sd_a;
        % Update lengths (V_sorted is the same as original Sigma)
        D_sorted = D_a; 
        r_j = diag(D_a);
    end

    % Minimum axis lengths constraint:

    if r_j(1) < MinAxisLengths
        flag_readjust(j) = true;
        % Compute Extension to be made over axes maintaining the aspect 
        % ratio (r_2_bar = Ratio*r_1_bar, r_1_bar < r_2_bar)
        r_bar = zeros(size(r_j));
        % Compute axis lengths offsets to achieve the MinAxisLength
        r_bar(1) = MinAxisLengths - r_j(1); 
        r_bar(2) = Ratio*r_bar(1);
        % New diagonal eigenvalues matrix with extended radius
        D_a = D_sorted + diag(r_bar);
        % Augmented covariance matrix
        Sd_a = V_sorted*D_a*V_sorted' / 3;
        Sigma_a(:,:,j) = Sd_a * Sd_a;
    end
    
    % If there was no re-adjustment, then the updated Sigma_a 
    % is the original Sigma
    if ~flag_readjust(j)
        Sigma_a(:,:,j) = Sigma(:,:,j);
    end

end

end

