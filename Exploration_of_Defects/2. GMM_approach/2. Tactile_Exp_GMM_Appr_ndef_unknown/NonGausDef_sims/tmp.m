close all
clear
clc

n = 2; 

L_1_l = 0.0;
L_1_u = 0.28;
dx_1 = (L_1_u - L_1_l)/50;

L_2_l = 0.0;
L_2_u = 0.20;
dx_2 = (L_2_u - L_2_l)/50;

% Dimensiones \mathbf{x} = [x_1 x_2]^T
x_1 = (L_1_l:dx_1:L_1_u)';
x_2 = (L_2_l:dx_2:L_2_u)';

% vector de límites inferior y superiores de las dimensiones
L_i_l = [L_1_l, L_2_l];
L_i_u = [L_1_u, L_2_u];

% Longitudes
L_1 = L_1_u - L_1_l;
L_2 = L_2_u - L_2_l;

[x_1_grid, x_2_grid] = meshgrid(x_1, x_2);

% Espacio de búsqueda discretizado
Omega = [reshape(x_1_grid,[],1), reshape(x_2_grid,[],1)];

% centroids
C_r_1 = [0.23, 0.165];   
% rotation angles
theta_r_1 = [0; 0];   
% Lenghts (base, altura)
rect_lens_1 = [0.003, 0.02];     

def_1 = rectangular_def(Omega, C_r_1, theta_r_1, rect_lens_1);

C_r_2 = [0.01, 0.165];   
theta_r_2 = 0;   
rect_lens_2 = [0.01, 0.01];    
def_2 = rectangular_def(Omega, C_r_2, theta_r_2, rect_lens_2);

C_r_3 = [0.12, 0.04;
         0.12, 0.036]; 
theta_r_3 = [0, 0];   
rect_lens_3 = [0.01, 0.003;
               0.003, 0.005];    
def_3 = rectangular_def(Omega, C_r_3, theta_r_3, rect_lens_3);

def_centers = [def_1.Center; def_2.Center; def_3.Center];

pi_1 = 1;
pi_2 = 1;
pi_3 = 1;

Phi_x = pi_1*def_1.pdf + pi_2*def_2.pdf + pi_3*def_3.pdf;

fig1h = figure(1);
tiledlayout(1,2)

nexttile
surf(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "FaceColor","interp","EdgeColor","none")
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Real PDF")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
grid on

nexttile
pcolor(x_1_grid, x_2_grid, reshape(Phi_x, length(x_2), length(x_1)),...
    "FaceColor","interp","EdgeColor","none")
hold on

plot(def_centers(:,1),def_centers(:,2), "r+", "MarkerSize",20, "LineWidth",5)
hold off
xlim([L_1_l, L_1_u])
ylim([L_2_l, L_2_u])
title("Real PDF")
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
axis equal tight
grid on
colorbar

set(findall(fig1h,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig1h,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(fig1h, "-property", "FontSize"), "FontSize", 20)