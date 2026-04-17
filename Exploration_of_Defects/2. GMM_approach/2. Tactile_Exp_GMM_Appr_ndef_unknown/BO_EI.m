close all
clear
clc

%% Parámetros del espacio de búsqueda U = [L_1_l, L_1_u] \times [L_2_l, L_2_u]

n = 2; % Número de dimensiones espaciales

L_1_l = 1;
L_1_u = 1.5;

L_2_l = 1;
L_2_u = 1.5;

L_1 = (L_1_u - L_1_l);
L_2 = (L_2_u - L_2_l);

dx_1 = L_1/100;
dx_2 = L_2/100;

% Dimensiones \mathbf{x} = [x_1 x_2]^T
x_1 = (L_1_l:dx_1:L_1_u)';
x_2 = (L_2_l:dx_2:L_2_u)';

%vector de límites inferior y superiores de las dimensiones
L_i_l = [L_1_l, L_2_l];
L_i_u = [L_1_u, L_2_u];

xrange = [L_i_l', L_i_u'];

[x_1_grid, x_2_grid] = meshgrid(x_1, x_2);

%Espacio de búsqueda discretizado
Omega = [reshape(x_1_grid,[],1), reshape(x_2_grid,[],1)]; 

%% Initial positions vector

X0 = [L_1_l, L_2_l;
      L_1_l + 1*L_1/5, L_2_l;
      L_1_l + 2*L_1/5, L_2_l;
      L_1_l + 3*L_1/5, L_2_l;
      L_1_l + 4*L_1/5, L_2_l;
      L_1_u, L_2_l;
      L_1_u, L_2_l + 1*L_2/5;
      L_1_u, L_2_l + 2*L_2/5;
      L_1_u, L_2_l + 3*L_2/5;
      L_1_u, L_2_l + 4*L_2/5;
      L_1_u, L_2_u;
      L_1_u - 1*L_1/5, L_2_u;
      L_1_u - 2*L_1/5, L_2_u;
      L_1_u - 3*L_1/5, L_2_u;
      L_1_u - 4*L_1/5, L_2_u;
      L_1_l, L_2_u;
      L_1_l, L_2_u - 1*L_2/5;
      L_1_l, L_2_u - 2*L_2/5;
      L_1_l, L_2_u - 3*L_2/5;
      L_1_l, L_2_u - 4*L_2/5];

%% Defects definition

n_def = 3;

% Save random number generator to reproduce Tests
rng_saved = rng;

% Defects generator
[Mu, Sigma, r_elips_Phi] = DefectsGen(n_def, L_i_l, L_i_u);

proportions = [1/(1+2+3),...
               2/(1+2+3),...
               3/(1+2+3)];
gm_dist = gmdistribution(Mu, Sigma, proportions);

%PDF de referencia REAL
Phi_x = pdf(gm_dist, Omega);

%% Tiempos
t_f = 20;           %Tiempo final por iteración
T_s = 0.1;                  % Tiempo de muestreo
t = (0:T_s:t_f)';   %Vector de tiempo por iteración

% Initial state
z_0 = [X0(1,1); 0; X0(1,2); 0]; 

% t for interpolation
freq_interp = 200;
T_s_interp = 1/freq_interp;
t_interp = (0:T_s_interp:T_s-T_s_interp)';

% Max velocity
v_max = 1.0;
% Max displacement
dmax = v_max*T_s;

%% Measurement parameters
a = 5; %Reference force
b = 0.2;
c = 0.05; %Amplitud del ruido en la medición 0.05

V_real = a + b*Phi_x;

% Initial position
X_e = X0(1,:);

% Generate a random first displacement within the velocity constraint
% dmax = v_max*T_s; % Maximum displacement
% validPoint = false;
% 
% while ~validPoint
% 
% Xrand = [(X_e(1,1) - dmax) + 2*dmax*rand, ...
%           (X_e(1,2) - dmax) + 2*dmax*rand];
% 
% v_init2 = sum(((Xrand - X_e)/T_s).^2, 2);
% 
% IsVelValid = v_init2 <= v_max^2;
% IsInOmega = all(isbetween(Xrand, L_i_l, L_i_u));
% 
% validPoint = IsVelValid && IsInOmega;
% 
% end

% Generate a random first displacement within the maximum velocity
% theta = linspace(0, 2*pi, 100)';
% P_tmp = [X_e(1,1) + dmax*cos(theta), X_e(1,2) + dmax*sin(theta)];
% IsInOmega = false;
% while ~IsInOmega
%     Xrand = P_tmp(randi(height(P_tmp)),:);
%     IsInOmega = all(isbetween(Xrand, L_i_l, L_i_u));
% end

% Initialize trajectory wit zero velocity
X_e = [X_e; X_e]; % Xrand];

% Interpolate to a 200Hz trajectory
x_e_sp = spline(t(1:height(X_e)), X_e(:,1), t_interp);
y_e_sp = spline(t(1:height(X_e)), X_e(:,2), t_interp);

X_e_sp = [x_e_sp, y_e_sp];

% Take the measurements
V = a + b*pdf(gm_dist, X_e_sp) + c*randn(height(X_e_sp), 1);

%% Plots

fig1h = figure(1);
tiledlayout(3,3)

nexttile(1)
TrueF_plot = surf(x_1_grid, x_2_grid, reshape(V_real, length(x_2), length(x_1)),...
                 "EdgeColor", "none", "FaceColor", "interp");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Ground Truth Function to be Optimized")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
grid on

nexttile(2)
Pred_plot = surf(x_1_grid, x_2_grid, reshape(zeros(height(Omega),1), length(x_2), length(x_1)),...
                "EdgeColor", "none", "FaceColor", "interp");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("GPR Prediction")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
grid on

nexttile(3)
pcolor(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "EdgeColor", "none", "FaceColor", "interp")
hold on
Traj_sp_plot = plot(X_e_sp(:,1), X_e_sp(:,2),'g','LineWidth', 2,...
                    'Marker','*', 'MarkerSize',6);
Traj_plot = plot(X_e(:,1), X_e(:,2),'k','LineWidth', 2,...
                'Marker','*', 'MarkerSize',8);
Traj_act_plot = plot(X_e(:,1), X_e(:,2),'red','LineWidth', 2,...
                'Marker','o', "MarkerFaceColor", "red", 'MarkerSize',8);
plot(X_e(1,1), X_e(1,2), 'rsq', 'LineWidth',4, 'MarkerSize',12)
hold off
axis equal
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Evolution in Omega")

nexttile(4)
sd_plot = surf(x_1_grid, x_2_grid, reshape(zeros(height(Omega),1), length(x_2), length(x_1)),...
        "EdgeColor", "none", "FaceColor", "interp");
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title('Uncertainty')

nexttile(5)
EI_plot = surf(x_1_grid, x_2_grid, reshape(zeros(height(Omega),1), length(x_2), length(x_1)),...
                "EdgeColor", "none", "FaceColor", "interp");
hold on
xPosEI_plot = plot3([1, 1],xrange(2,:),0*[1 1],'--k','LineWidth',1.5);
yPosEI_plot = plot3(xrange(1,:),[1, 1],0*[1 1],'--k','LineWidth',1.5); 
hold off; 
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title('Expected Improvement'); 
grid on;

nexttile(6)
Est_plot = pcolor(x_1_grid, x_2_grid, reshape(zeros(height(Omega),1), length(x_2), length(x_1)),...
                 "EdgeColor", "none", "FaceColor", "interp");
axis equal
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Estimation")

nexttile(7, [1, 2])
loglikelihood = 0;
loglikelihood_plot = plot(t(1), loglikelihood, "LineWidth", 2);
title("Log-likelihood")
xlabel('Time [s]')
grid on

% nexttile(8)
% d_likelihood = 1;
% d_likelihood_plot = plot(t(1), d_likelihood, "LineWidth", 2);
% title("change in Log-likelihood")
% xlabel('Time [s]')
% grid on
% 
% nexttile(9)
% Loss = 1;
% loss_plot = plot(t(1), Loss, "LineWidth", 2);
% title("Loss")
% xlabel('Time [s]')
% grid on

set(findall(fig1h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig1h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig1h, "-property", "FontSize"), "FontSize", 20)

%% Loop

for i = 2:length(t)-1

    % avgV = mean(V);
    % stdV = std(V);
    % 
    % if stdV <= 1e-8 || isnan(stdV)
    %     stdV = 1;
    % end
    % 
    % V_normalized = (V - avgV)/stdV;
    % mdl = fitrgp(X_e_sp, V_normalized);

    mdl = fitrgp(X_e_sp, V,...
                'KernelFunction', 'squaredexponential',...
                'Standardize',true,...
                'SigmaLowerBound',0.1,...
                'BasisFunction','constant',...
                'PredictMethod','exact');
    [V_pred, sd] = predict(mdl, Omega);
    % Loss = resubLoss(mdl);

    %% Expected Improvement
    % This EI is from http://krasserm.github.io/2018/03/21/bayesian-optimization/
    
    xi = 0;  % Exploration-exploitation parameter (greek letter, xi)
             % High xi = more exploration
             % Low xi = more exploitation (can be < 0)

    d = V_pred - max(V) - xi; % (y - f*) if maximization

    EI = (sd ~= 0).*(d.*normcdf(d./sd) + sd.*normpdf(d./sd));

    % [eimax,posEI] = max(EI); 
    % xEI = Omega(posEI,:);
    
    % extract next feasible points (constraint)
    idx_feasible = sum((Omega - X_e(end,:)).^2, 2) <= dmax^2;

    feasible_Points = Omega(idx_feasible, :);
    feasible_EI = EI(idx_feasible);

    % sort EI (first element is the maximum)
    [eisorted, idx_ei] = sort(feasible_EI, 1, "descend");

    % re-arrange Omega
    feasiblePtsSortedEI = feasible_Points(idx_ei,:);

    % Algorithm for selecting the next point X
    for j = 1:height(eisorted)
        % save the maximum EI (first elements)
        eimax = eisorted(j);
        % check for equal EI values
        idxEqEI = eimax == eisorted;
        % Extract the set of posible x's in which EI is maximum
        posibleX = feasiblePtsSortedEI(idxEqEI, :);
        % Select a random X from the set (if there are many)
        idx_rand = randi(height(posibleX));
        xEI = posibleX(idx_rand, :);
        % check for already visited point, if already visited, iterate
        % again
        flag_visitedP = ismember(xEI, X_e, "rows");
        if ~flag_visitedP
            break;
        end
    end
    
    % % sort EI (first element is the maximum)
    % [eisorted, idx_ei] = sort(feasible_EI, 1, "descend");
    % 
    % % Avoid already visited points
    % for j = 1:height(eisorted)
    %     eimax = eisorted(j);
    %     xEI = feasible_Points(idx_ei(j), :);
    %     % check for already visited point
    %     XYcheck = X_e == xEI;
    %     check = XYcheck(:,1) & XYcheck(:,2);
    %     if ~any(check)
    %         break;
    %     end
    % end

    X_e(end + 1,:) = xEI;

    % new interpolation
    t_interp = ((i-1)*T_s:T_s_interp:i*T_s-T_s_interp)';
    x_e_sp = spline(t(i:i+1), X_e(i:i+1,1), t_interp);
    y_e_sp = spline(t(i:i+1), X_e(i:i+1,2), t_interp);
    X_e_sp_new = [x_e_sp, y_e_sp];
    X_e_sp(end+1:end+height(X_e_sp_new), :) = X_e_sp_new;
    % Update the measurement vector
    V_new = a + b*pdf(gm_dist, X_e_sp_new) + c*randn(height(X_e_sp_new), 1);
    V = [V; 
        V_new];
    
    i

    %% update plots

    Pred_plot.ZData = reshape(V_pred, length(x_2), length(x_1));
    sd_plot.ZData = reshape(sd, length(x_2), length(x_1));
    EI_plot.ZData = reshape(EI, length(x_2), length(x_1));
    Est_plot.CData = reshape(V_pred, length(x_2), length(x_1));
    xPosEI_plot.XData = xEI([1, 1]);
    xPosEI_plot.ZData = eimax*[1 1];
    yPosEI_plot.YData = xEI([2, 2]);
    yPosEI_plot.ZData = eimax*[1 1];
    
    Traj_sp_plot.XData(end + 1 : end + height(x_e_sp)) = x_e_sp';
    Traj_sp_plot.YData(end + 1 : end + height(y_e_sp)) = y_e_sp';
    Traj_plot.XData(end + 1) = X_e(end,1);
    Traj_plot.YData(end + 1) = X_e(end,2);
    Traj_act_plot.XData = X_e(end,1);
    Traj_act_plot.YData = X_e(end,2);

    % Update the log-likelihood plot

    loglikelihood = [loglikelihood;...
                    mdl.LogLikelihood];
    % changeLogL = abs(loglikelihood(end) - loglikelihood(end-1));
    % d_likelihood = [d_likelihood;...
    %                 changeLogL];

    % d_likelihood_plot.XData(end + 1) = t(i);
    % d_likelihood_plot.YData(end + 1) = d_likelihood(end);
    loglikelihood_plot.XData(end + 1) = t(i);
    loglikelihood_plot.YData(end + 1) = loglikelihood(end);
    % loss_plot.XData(end + 1) = t(i);
    % loss_plot.YData(end + 1) = loss(mdl, X_e_sp, V); %Loss

    % Exit the loop when the LogLikelihood is low enough
    if mdl.LogLikelihood <= -1000
        disp("LogLikelihood reached the condition")
        break;
    end

    pause(0.1)

end

%% More plots xd

t_plot = linspace(0, T_s*i, height(X_e_sp))';

fig2h = figure(2);
tiledlayout(3,1)

nexttile(1)
% plot(t(1:length(X_e)), X_e, "LineWidth", 3)
plot(t_plot, X_e_sp, "LineWidth", 3)
xlabel('Time [s]')
ylabel('Position [m]')
title("Position")
legend("$x_{e_1}$", "$x_{e_2}$")

v_normplot = sqrt(sum((diff(X_e_sp)/T_s_interp).^2, 2));
v_normplot(end + 1) = v_normplot(end);
nexttile(2)
% plot(t(1:length(v_normplot)), v_normplot, "LineWidth", 3)
plot(t_plot, v_normplot, "LineWidth", 3)
xlabel('Time [s]')
ylabel('Velocity [m/s]')
title('Velocity norm')
grid on;

set(findall(fig2h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig2h,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fig2h, "-property", "FontSize"), "FontSize", 20)