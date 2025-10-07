function [phi_k_reg, f_k_reg, h_k_reg] = FourierCoef_RefPDF(Phi_hat_x, Par_struct)

%K, n, K_cal, Omega, dx_1, dx_2, L_i_l, L_i_u

K = Par_struct.K;
n = Par_struct.n;
K_cal = Par_struct.K_cal;
Omega = Par_struct.Omega;
dx_1 = Par_struct.dx_1;
dx_2 = Par_struct.dx_2;
L_i_l = Par_struct.L_i_l;
L_i_u = Par_struct.L_i_u;

% Registro para guardar los coeficientes
phi_k_reg = zeros(K^n, 1);

% Registro para guardar las funciones de Fourier ortonormales
f_k_reg = zeros(height(Omega), K^n);

% Registro para el término normalizador
h_k_reg = zeros(K^n,1);

% Cálculo de los coeficientes
for j = 1:K^n
    
    k_vect_j = K_cal(:,j)';

    % función ortogonal de Fourier
    f_tilde_k_j = prod(cos( k_vect_j.*pi.*(Omega - L_i_l)./(L_i_u - L_i_l) ), 2);

    % Término normalizador = escalar por cada k
    h_k_j = sqrt(sum(f_tilde_k_j .^2) * dx_1 * dx_2);

    % función ortonormal de Fourier
    f_k_j = f_tilde_k_j ./ h_k_j;

    % Coeficientes de Fourier, la integral se aproxima con Riemann
    phi_k_j = sum(Phi_hat_x .* f_k_j)*dx_1*dx_2;

    % Registros
    h_k_reg(j) = h_k_j;
    f_k_reg(:,j) = f_k_j;
    phi_k_reg(j) = phi_k_j;

end




end