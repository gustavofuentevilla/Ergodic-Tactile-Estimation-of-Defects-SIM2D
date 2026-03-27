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
t_f = 10;           %Tiempo final por iteración
T_s = 0.1;                  % Tiempo de muestreo
t = (0:T_s:t_f)';   %Vector de tiempo por iteración

% Initial state
z_0 = [X0(1,1); 0; X0(1,2); 0]; 

% t for interpolation
freq_interp = 200;
t_interp = (0:1/freq_interp:T_s)';

%% Measurement parameters
a = 5; %Reference force
b = 0.2;
c = 0.05; %Amplitud del ruido en la medición

V_real = a + b*Phi_x;

figh = figure;
tiledlayout(2,2)
nexttile(1)
surf(x_1_grid, x_2_grid, reshape(V_real, length(x_2), length(x_1)),...
    "EdgeColor", "none", "FaceColor", "interp")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Ground Truth Function to be Optimized")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
grid on

% Random initial positions
P_start = 2;

% X_e = [rand(P_start,1)*diff(xrange(1,:))+xrange(1,1),...
%      rand(P_start,1)*diff(xrange(2,:))+xrange(2,1)]; 

X_e = zeros(2, 2);
X_e(1,:) = Mu(1,:);
X_e(2,:) = Mu(1,:) + [0.005, 0.005];

x_e_sp = spline(t(1:2), X_e(:,1), t_interp);
y_e_sp = spline(t(1:2), X_e(:,2), t_interp);

X_e_sp = [x_e_sp, y_e_sp];

% V = zeros(size(X_e,1),1);

% Take the measurement at X_e
V = a + b*pdf(gm_dist, X_e_sp) + c*randn(height(X_e_sp), 1);


%% Loop

for i = 2:length(t)-1

    mdl = fitrgp(X_e_sp, V, 'KernelFunction', 'ardmatern52');
    [V_pred, sd] = predict(mdl, Omega);

    nexttile(2)
    surf(x_1_grid, x_2_grid, reshape(V_pred, length(x_2), length(x_1)),...
        "EdgeColor", "none", "FaceColor", "interp")
    xlabel('$x_1$ [m]')
    ylabel('$x_2$ [m]')
    zlabel('Predicted Value')
    grid on
    hold on
    scatter3(X_e_sp(:,1), X_e_sp(:,2), V, 20,...
        'k','filled','MarkerEdgeColor','k')
    title("Observed points")
    hold off

    nexttile(3)
    surf(x_1_grid, x_2_grid, reshape(sd, length(x_2), length(x_1)),...
        "EdgeColor", "none", "FaceColor", "interp")
    title('Uncertainty');

    %% Expected Improvement
    % This EI is from http://krasserm.github.io/2018/03/21/bayesian-optimization/
    
    xi = 3;  % Exploration-exploitation parameter (greek letter, xi)
             % High xi = more exploration
             % Low xi = more exploitation (can be < 0)

    d = V_pred - max(V) - xi; % (y - f*) if maximization

    EI = (sd ~= 0).*(d.*normcdf(d./sd) + sd.*normpdf(d./sd));

    [eimax,posEI] = max(EI); 
    xEI = Omega(posEI,:);
    X_e(end+1,:) = xEI;

    % new interpolation
    t_interp = ((i-1)*T_s:1/freq_interp:i*T_s)';
    x_e_sp = spline(t(i:i+1), X_e(i:i+1,1), t_interp);
    y_e_sp = spline(t(i:i+1), X_e(i:i+1,2), t_interp);
    X_e_sp(end+1:end+height(x_e_sp), :) = [x_e_sp, y_e_sp];
    % Update the measurement vector
    V = a + b*pdf(gm_dist, X_e_sp) + c*randn(height(X_e_sp), 1);

    nexttile(4)
    surf(x_1_grid, x_2_grid, reshape(EI, length(x_2), length(x_1)),...
        "EdgeColor", "none", "FaceColor", "interp")
    hold on
    plot3(xEI([1 1]),xrange(2,:),eimax*[1 1],'--k','LineWidth',1.5);
    plot3(xrange(1,:),xEI([2 2]),eimax*[1 1],'--k','LineWidth',1.5);
    title('Expected Improvement'); 
    grid on; 
    hold off; 
    
    i

    % if mean(sd) < 2
    %     break;
    % end

    pause(0.1)

end
%%
figure
pcolor(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "EdgeColor", "none", "FaceColor", "interp")
hold on
plot(X_e_sp(:,1), X_e_sp(:,2),'g','LineWidth', 2, 'Marker','*', 'MarkerSize',6)
plot(X_e(:,1), X_e(:,2),'k','LineWidth', 2, 'Marker','*', 'MarkerSize',8)
plot(X_e(1,1), X_e(1,2), 'rsq', 'LineWidth',4, 'MarkerSize',12)
hold off
axis equal
xlim(xrange(1,:))
ylim(xrange(2,:))



