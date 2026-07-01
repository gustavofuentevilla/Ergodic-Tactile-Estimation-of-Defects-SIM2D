function Estim_sol = BO_Postprocessing(X_e_sp, V, Par_PDF)

thres_meas = Par_PDF.thres_meas;
D_max = Par_PDF.D_max;
Thres_Variation = Par_PDF.Thres_Variation;
OneClustDistLimit = Par_PDF.OneClustDistLimit;

% Locate values above threshold
idx_V = V > thres_meas; 

Data = [X_e_sp(idx_V,:), V(idx_V)];

% V Conversion to int
V_int = round(Data(:,3));    
% Repeat elements on spatial domain (trajectory points)
Data_Xe_hist_V = repelem(Data(:,1:2), V_int, 1); 
% Without repeating position elements:
% Data_Xe_hist_V = Data(:,1:2); 

flag_NoData = 0;
if isempty(Data)
    flag_NoData = 1;
    % Output structure
    Estim_sol.Data = Data;
    Estim_sol.Data_Xe_hist_V = Data_Xe_hist_V;
    Estim_sol.numComponents = 0;
    Estim_sol.Priors_found = [];
    Estim_sol.Mu_found = [];
    Estim_sol.Sigma_found = [];
    Estim_sol.Sigma_found_a = [];
    Estim_sol.flag_readjust = [];
    Estim_sol.flag_NoData = flag_NoData;
    Estim_sol.flag_OneCluster = false;
    return;
end

% %%%%%%%%%%%%% Logic for one single cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check distance among X_e Points in the plane and test:
% If the maximum distance of the CURRENT batch is less than a certain 
% limit, then we are sure there is only one cluster (defect).
% *flag_OneCluster is used to avoid the computation of more than one
% clusters logic below
flag_OneCluster = false;

d_max = max(pdist(Data(:,1:2)));
if d_max <= OneClustDistLimit
    % One cluster
    clust_eval = [];
    numComponents = 1;
    Model = GMM_EM(Data_Xe_hist_V, numComponents,...
                   "Max_iter", 500, "Min_var", 1e-10);
    flag_OneCluster = true;
end


% %%%%%%%%%%%%% Logic for multiple clusters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~flag_OneCluster
    % Evaluate for optimal number of clusters
    D_max = 10;
    clust_eval = evalclusters(Data_Xe_hist_V,"kmeans","silhouette",...
                                      KList=2:D_max);
    numComponents = clust_eval.OptimalK;
end

% Define the model
Model = GMM_EM(Data_Xe_hist_V, numComponents,...
                   "Max_iter", 500, "Min_var", 1e-10);

% Extracting model parameters
Mu_found = Model.Mu;
Sigma_found = Model.Sigma;
Priors_found = Model.Priors;

[Sigma_found_a, flag_readjust] = ReadjustConstrainedSigma(Sigma_found,...
                                    Thres_Variation, 0);


% %%%%%%%%%%%%%%%%%% Output structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Estim_sol.Data = Data;
Estim_sol.Data_Xe_hist_V = Data_Xe_hist_V;
Estim_sol.numComponents = numComponents;
Estim_sol.Priors_found = Priors_found;
Estim_sol.Mu_found = Mu_found;
Estim_sol.Sigma_found = Sigma_found;
Estim_sol.Sigma_found_a = Sigma_found_a;
Estim_sol.flag_readjust = flag_readjust;
Estim_sol.flag_NoData = flag_NoData;
Estim_sol.flag_OneCluster = flag_OneCluster;


end