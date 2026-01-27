close all
clearvars -except Erg_traj_ipopt
clc

%%
% Case (1: full algoritm, 2: No repeating points, 3: No removing gaussians)
Case_i = 2;
% Number of Defectss
def = 3;
% Number of initial positions
totalX0 = 20;
% Number of tests per initial position
totalOut = 5;

% Buffers for:
nbIter = zeros(totalX0,totalOut);   % Number of iterations
nbMu = zeros(totalX0,totalOut);     % Number of found defects
Tsim = zeros(totalX0,totalOut);     % Simulation times

for X0_i = 1:totalX0
    for out_i = 1:totalOut
        filename = "Test2/Case" + Case_i +...
                   "/" + def + "Def/X0_" + X0_i + ...
                    "/output_" + out_i + ".mat";
        load(filename, "n_iter", "Mu_found", "T_sim_toc")

        nbIter(X0_i, out_i) = n_iter;
        nbMu(X0_i, out_i) = height(Mu_found);
        Tsim(X0_i, out_i) = T_sim_toc;
    end
end

% Success rate (Number of tests in which all the defects are found)
successRate = sum(nbMu == def, 'all') / numel(nbMu) * 100;

% Calculate the average values and standard deviations
avgIter = mean(nbIter, "all");
avgTsim = mean(Tsim, "all");

stdIter = std(nbIter, 0, "all");
stdTsim = std(Tsim, 0, "all");