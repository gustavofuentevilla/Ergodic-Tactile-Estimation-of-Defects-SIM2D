function [Phi_hat_x_1_next, Phi_hat_x_2_next, Estim_sol] = PDF_Estimator(Phi_hat_x_1_last, Phi_hat_x_2_last, X_e, V_hat, Par_PDF)

x_1 = Par_PDF.x_1;
x_2 = Par_PDF.x_2;
MinLength_x_j = Par_PDF.MinLength_x_j;
% a = Par_PDF.Meas_mean;

% Likelihood function for Theta given V_hat for one defect (functions projected on x_1 and x_2)

exp_x1 = sum( V_hat.*X_e(:,1) ) / sum(V_hat);
var_x1 = sum( V_hat.*(X_e(:,1).^2) ) / sum(V_hat) - exp_x1^2;
P_V_k_x1 = (1/(sqrt(2*pi*var_x1)))*exp( -(x_1 - exp_x1).^2/(2*var_x1));

exp_x2 = sum( V_hat.*X_e(:,2)) / sum(V_hat);
var_x2 = sum( V_hat.*(X_e(:,2).^2)) / sum(V_hat) - exp_x2^2;
P_V_k_x2 = (1/(sqrt(2*pi*var_x2)))*exp( -(x_2 - exp_x2).^2/(2*var_x2));

% Previous PDF
P_theta_x1 = Phi_hat_x_1_last;
P_theta_x2 = Phi_hat_x_2_last;

% Next PDF

P_theta_x1_next = P_V_k_x1.*P_theta_x1;
exp_x1_next = sum( P_theta_x1_next.*x_1) / sum(P_theta_x1_next);
var_x1_next = sum( P_theta_x1_next.*(x_1.^2)) / sum(P_theta_x1_next) - exp_x1_next^2;

P_theta_x2_next = P_V_k_x2.*P_theta_x2;
exp_x2_next = sum( P_theta_x2_next.*x_2) / sum(P_theta_x2_next);
var_x2_next = sum( P_theta_x2_next.*(x_2.^2)) / sum(P_theta_x2_next) - exp_x2_next^2;

if 3*sqrt(var_x1_next) < MinLength_x_j
    var_x1_next = (MinLength_x_j/3)^2;
end
if 3*sqrt(var_x2_next) < MinLength_x_j
    var_x2_next = (MinLength_x_j/3)^2;
end

x1_length = 3*sqrt(var_x1_next);
x2_length = 3*sqrt(var_x2_next);

Phi_hat_x_1_next = normpdf(x_1, exp_x1_next, sqrt(var_x1_next));
Phi_hat_x_2_next = normpdf(x_2, exp_x2_next, sqrt(var_x2_next));

%% Condition for a found defect



% variation_Phi_hat = sd_x1_next + sd_x2_next;
% 
% Cond = (variation_Phi_hat < Thres_PostVariation);
% %%%%%%%%%%%%%%%%%%%%%%%%% ME QUEDE AQUI XD
% Def_found.mu = [];
% Def_found.sigma = [];
% 
% If true --> One or more Defects have been found
% if any(Cond) 
%     % Index of defects found ( idx_def = find(Cond) )
%     idx_def = Cond;  
%     % Saving Gaussian parameters related
%     Def_found.mu = Mu(idx_def,:);
%     Def_found.sigma = Sigma(:,:,idx_def);
%     % Removing defect parameters from GMM solution
%     Mu(idx_def,:) = [];
%     Sigma(:,:,idx_def) = [];
% end
% 
% % Check if we've done
% flag_done = isempty(Mu);
% 
% % If we've done, export a zero PDF 
% if flag_done
%     Phi_hat_x_next = zeros(height(Omega));
%     Estim_sol.GMModel = [];
% else % IF we've not done, compute and export next GMM
%     % Next PDF
%     GMModel = gmdistribution(Mu, Sigma, Priors); %This function normalizes Priors if these doesn't sum to 1.
%     Phi_hat_x_next = pdf(GMModel, Omega);
%     % Export
%     Estim_sol.GMModel = GMModel;
% end

%% Output structure

Estim_sol.exp_x1_V_Xe = exp_x1;
Estim_sol.var_x1_V_Xe = var_x1;
Estim_sol.exp_x2_V_Xe = exp_x2;
Estim_sol.var_x2_V_Xe = var_x2;
Estim_sol.exp_x1_hat = exp_x1_next;
Estim_sol.var_x1_hat = var_x1_next;
Estim_sol.exp_x2_hat = exp_x2_next;
Estim_sol.var_x2_hat = var_x2_next;
Estim_sol.x1_length = x1_length;
Estim_sol.x2_length = x2_length;

end