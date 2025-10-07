function [Phi_hat_x_1_next, Phi_hat_x_2_next] = PDF_Estimator(Phi_hat_x_1_last, Phi_hat_x_2_last, X_e, V_Xe, Par_PDF)

% mu_1 = Par_PDF.mu_1;
% sigma_1 = Par_PDF.sigma_1;
% mu_2 = Par_PDF.mu_2;
% sigma_2 = Par_PDF.sigma_2;
% N = Par_PDF.N;
x_1 = Par_PDF.x_1;
x_2 = Par_PDF.x_2;
a = Par_PDF.Meas_mean;

% Likelihood function for Theta given measurements (functions projected on x_1 and x_2)
Meas = V_Xe - a;
Meas = Meas - min(Meas);

exp_x1 = sum( Meas.*X_e(:,1) ) / sum(Meas);
var_x1 = sum( Meas.*(X_e(:,1).^2) ) / sum(Meas) - exp_x1^2;
P_V_k_x1 = (1/(sqrt(2*pi*var_x1)))*exp( -(x_1 - exp_x1).^2/(2*var_x1));

exp_x2 = sum( Meas.*X_e(:,2)) / sum(Meas);
var_x2 = sum( Meas.*(X_e(:,2).^2)) / sum(Meas) - exp_x2^2;
P_V_k_x2 = (1/(sqrt(2*pi*var_x2)))*exp( -(x_2 - exp_x2).^2/(2*var_x2));

% Previous PDF
P_theta_x1 = Phi_hat_x_1_last;
P_theta_x2 = Phi_hat_x_2_last;

% Next PDF
P_theta_x1_next = P_V_k_x1.*P_theta_x1;
exp_x1_next = sum( P_theta_x1_next.*x_1) / sum(P_theta_x1_next);
var_x1_next = sum( P_theta_x1_next.*(x_1.^2)) / sum(P_theta_x1_next) - exp_x1_next^2;
Phi_hat_x_1_next = normpdf(x_1, exp_x1_next, sqrt(var_x1_next));

P_theta_x2_next = P_V_k_x2.*P_theta_x2;
exp_x2_next = sum( P_theta_x2_next.*x_2) / sum(P_theta_x2_next);
var_x2_next = sum( P_theta_x2_next.*(x_2.^2)) / sum(P_theta_x2_next) - exp_x2_next^2;
Phi_hat_x_2_next = normpdf(x_2, exp_x2_next, sqrt(var_x2_next));

%% Unrealistic update

% mu_1 = Par_PDF.mu_1;
% sigma_1 = Par_PDF.sigma_1;
% mu_2 = Par_PDF.mu_2;
% sigma_2 = Par_PDF.sigma_2;
% N = Par_PDF.N;
% x_1 = Par_PDF.x_1;
% x_2 = Par_PDF.x_2;
% 
% % Measurement model V
% a = 5;
% b = 0.5;
% c = 0.1;
% Upsilon_1 = a + b*normpdf(X_e(:,1), mu_1, sigma_1);    
% delta_1 = c*randn(N+1, 1);       % Gaussian noise on x_1
% Upsilon_2 = a + b*normpdf(X_e(:,2), mu_2, sigma_2);    
% delta_2 = c*randn(N+1, 1);       % Gaussian noise on x_2
% %noise_var = c^2;
% 
% % Mediciones sobre x_1 y x_2
% V_x_1 = Upsilon_1 + delta_1; %b*rand(N+1, 1).*
% V_x_2 = Upsilon_2 + delta_2;
% 
% % Likelihood function for Theta given measurements (functions projected on x_1 and x_2)
% temp_1 = V_x_1 - a;
% temp_1 = temp_1 - min(temp_1);
% exp_x1 = sum( temp_1.*X_e(:,1) ) / sum(temp_1);
% var_x1 = sum( temp_1.*(X_e(:,1).^2) ) / sum(temp_1) - exp_x1^2;
% P_V_k_x1 = (1/(sqrt(2*pi*var_x1)))*exp( -(x_1 - exp_x1).^2/(2*var_x1));
% 
% temp_2 = V_x_2 - a;
% temp_2 = temp_2 - min(temp_2);
% exp_x2 = sum( temp_2.*X_e(:,2)) / sum(temp_2);
% var_x2 = sum( temp_2.*(X_e(:,2).^2)) / sum(temp_2) - exp_x2^2;
% P_V_k_x2 = (1/(sqrt(2*pi*var_x2)))*exp( -(x_2 - exp_x2).^2/(2*var_x2));
% 
% % Previous PDF
% P_theta_x1 = Phi_hat_x_1_last;
% P_theta_x2 = Phi_hat_x_2_last;
% 
% % Next PDF
% P_theta_x1_next = P_V_k_x1.*P_theta_x1;
% exp_x1_next = sum( P_theta_x1_next.*x_1) / sum(P_theta_x1_next);
% var_x1_next = sum( P_theta_x1_next.*(x_1.^2)) / sum(P_theta_x1_next) - exp_x1_next^2;
% Phi_hat_x_1_next = normpdf(x_1, exp_x1_next, sqrt(var_x1_next));
% 
% P_theta_x2_next = P_V_k_x2.*P_theta_x2;
% exp_x2_next = sum( P_theta_x2_next.*x_2) / sum(P_theta_x2_next);
% var_x2_next = sum( P_theta_x2_next.*(x_2.^2)) / sum(P_theta_x2_next) - exp_x2_next^2;
% Phi_hat_x_2_next = normpdf(x_2, exp_x2_next, sqrt(var_x2_next));


end