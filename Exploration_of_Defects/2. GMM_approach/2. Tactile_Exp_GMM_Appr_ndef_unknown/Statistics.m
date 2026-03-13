close all
clearvars -except Erg_traj_ipopt
clc

%%
% Case (1: full algoritm, 2: No repeating points, 3: No removing gaussians)
Case_i = 2;
% Number of Defectss
def = 7;
% Number of initial positions
totalX0 = 20;
% Number of tests per initial position
totalOut = 5;

% Buffers for:
nbIter = zeros(totalX0,totalOut);   % Number of iterations
nbMu = zeros(totalX0,totalOut);     % Number of found defects
Tsim = zeros(totalX0,totalOut);     % Simulation times
T_traj = cell(totalX0,totalOut);    % Time spend computing ergodic trajectory
T_PDF = cell(totalX0,totalOut);     % Time spend computing the PDF
T_iter = cell(totalX0,totalOut);    % Times per iteration
Mu_buffer = cell(totalX0,totalOut);
Mu_found_buffer = cell(totalX0,totalOut);

checkidx = false(totalX0,totalOut);

for X0_i = 1:totalX0
    for out_i = 1:totalOut
        % Loading file and extract variables
        filename = "Test2/Case" + Case_i +...
                   "/" + def + "Def/X0_" + X0_i + ...
                    "/output_" + out_i + ".mat";
        load(filename, "n_iter", "Mu_found", "T_sim_toc", "Estim_sol",...
             "T_ErgC_i", "T_PDF_i", "Mu")
        
        % Save data into matrix (row = initial position; column = test)
        nbIter(X0_i, out_i) = n_iter;
        nbMu(X0_i, out_i) = height(Mu_found);

        % Checking
        checkidx(X0_i, out_i) = height(Estim_sol(end).Mu_found) < def;

        % Times
        Tsim(X0_i, out_i) = T_sim_toc;
        T_traj{X0_i, out_i} = T_ErgC_i;
        T_PDF{X0_i, out_i} = T_PDF_i;
        T_iter{X0_i, out_i} = T_ErgC_i + T_PDF_i;

        % Defect locations
        Mu_buffer{X0_i, out_i} = Mu;
        Mu_found_buffer{X0_i, out_i} = Mu_found;
    end
end

% Success rate (Number of tests in which all the defects are found)
successRate = sum(nbMu == def, 'all') / numel(nbMu) * 100;

% Success rate (porcentage of overall found defects)
TotalFoundDef = 0;
for X0_i = 1:totalX0
    for out_i = 1:totalOut
        TotalFoundDef = TotalFoundDef + height(Mu_found_buffer{X0_i, out_i});
    end
end
SR_def = TotalFoundDef/(def*totalX0*totalOut) * 100;

% Calculate the average values and standard deviations
avgIter = mean(nbIter, "all");
avgTsim = mean(Tsim, "all");

stdIter = std(nbIter, 0, "all");
stdTsim = std(Tsim, 0, "all");

% Simulation time per iteration
TimePerIter = [];

for X0_i = 1:totalX0
    for out_i = 1:totalOut
        TimePerIter = [TimePerIter; 
                       T_iter{X0_i, out_i}(1:nbIter(X0_i, out_i))];
    end
end

avgTimePerIter = mean(TimePerIter);
stdTimePerIter = std(TimePerIter);

% Mean error distance between estimated and real defect locations
e_centers = [];
for X0_i = 1:totalX0
    for out_i = 1:totalOut
        % If found defects buffer its not empty
        if ~isempty(Mu_found_buffer{X0_i, out_i})
            k_idx = dsearchn(Mu_buffer{X0_i, out_i},...
                             Mu_found_buffer{X0_i, out_i});
            Mu_tmp = Mu_buffer{X0_i, out_i}(k_idx, :);
            e_centers = [e_centers;...
                         sum((Mu_tmp - Mu_found_buffer{X0_i, out_i}).^2, 2).^0.5];
        end
    end
end

avgLocationError = mean(e_centers);
stdLocationError = std(e_centers);