function [Phi_hat_x_1_next, Phi_hat_x_2_next] = PDF_Estimator(Phi_hat_x_1_last, Phi_hat_x_2_last, X_e, V_Xe, Par_PDF)

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
var_x1_next = var_x1_next + 1e-6; % Add tiny variance (So it never reaches zero)
Phi_hat_x_1_next = normpdf(x_1, exp_x1_next, sqrt(var_x1_next));

P_theta_x2_next = P_V_k_x2.*P_theta_x2;
exp_x2_next = sum( P_theta_x2_next.*x_2) / sum(P_theta_x2_next);
var_x2_next = sum( P_theta_x2_next.*(x_2.^2)) / sum(P_theta_x2_next) - exp_x2_next^2;
var_x2_next = var_x2_next + 1e-6; % Add tiny variance (So it never reaches zero)
Phi_hat_x_2_next = normpdf(x_2, exp_x2_next, sqrt(var_x2_next));

end