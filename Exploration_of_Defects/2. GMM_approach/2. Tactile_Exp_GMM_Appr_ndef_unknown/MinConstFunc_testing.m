%% Funcion tanh para la definición de la restricción de minima variación
t_vect = (-1:0.1:10)';
n_iter_vect = (1:n_iter)';

D_KL = zeros(n_iter, 1);
for i = 1:n_iter
    D_KL(i) = Estim_sol(i).D_KL;
end

L_1 = (L_1_u - L_1_l);
L_2 = (L_2_u - L_2_l);

MaxVarCons = L_1 + L_2;
D_KL_bar_l = 1;         % D_KL value that matches the variation Threshold
D_KL_bar_u = 3.3;       % D_KL that matches certain porcentage (nu_p) 
                        % of maximum variation constraint (L1 + L2)
nu_p = 0.75;            % Porcentage of max variation constraint,
                        % namely, nu_p*(L_1 + L_2)
Variation_thres = Par_PDF.Thres_Variation;

% Resolver sistema de ecuaciones no lineales para alpha_i

F = @(alpha) [tanh(alpha(1)*D_KL_bar_l - alpha(2)) + 1 - 2*Variation_thres/MaxVarCons;
              tanh(alpha(1)*D_KL_bar_u - alpha(2)) + 1 - 2*nu_p];

alpha0 = [1; 3];

alpha = fsolve(F, alpha0);

f = 0.5*MaxVarCons*(tanh(alpha(1)*D_KL - alpha(2)) + 1);

figure(101)
plot(n_iter_vect, f, "o-", "LineWidth", 2, "MarkerSize", 8)
grid on

% Misc
func1 = 0.5*MaxVarCons*(tanh(alpha(1)*t_vect - alpha(2)) + 1);
figure(102)
plot(t_vect, func1, "LineWidth", 2)
grid on

%% Función saturación

nu_p = 0.9;

if D_KL <= D_KL_bar_l
    f_2 = Variation_thres;
elseif D_KL >= D_KL_bar_u
    f_2 = nu_p*MaxVarCons;
else
    f_2 = (nu_p*MaxVarCons - Variation_thres)*(D_KL - D_KL_bar_l)...
           /(D_KL_bar_u - D_KL_bar_l) + Variation_thres;
end

figure(201)
plot(n_iter_vect, f_2, "o-", "LineWidth", 2, "MarkerSize", 8)
grid on

% Misc
func2 = zeros(size(t_vect));
for i_tmp = 1:length(t_vect)
    if t_vect(i_tmp) <= D_KL_bar_l
        func2(i_tmp) = Variation_thres;
    elseif t_vect(i_tmp) >= D_KL_bar_u
        func2(i_tmp) = nu_p*MaxVarCons;
    else
        func2(i_tmp) = (nu_p*MaxVarCons - Variation_thres)*(t_vect(i_tmp) - D_KL_bar_l)...
               /(D_KL_bar_u - D_KL_bar_l) + Variation_thres;
    end
end

figure(202)
plot(t_vect, func2, "LineWidth", 2)
grid on

