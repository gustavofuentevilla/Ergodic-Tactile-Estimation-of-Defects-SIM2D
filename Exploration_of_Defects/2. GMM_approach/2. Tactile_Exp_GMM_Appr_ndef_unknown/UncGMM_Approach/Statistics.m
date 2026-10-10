close all
clearvars -except Erg_traj_ipopt
clc

%%
% Test (2: Same severity levels, 3: Different severity levels)}
Test_i = 3;
% Case (1: TEDD full algoritm, 2: TEDD-R No repeating points, 3: TEDD-T No removing gaussians)
% Case_i = 4;
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
Sigma_buffer = cell(totalX0,totalOut);
Sigma_found_buffer = cell(totalX0,totalOut);

checkidx = false(totalX0,totalOut);

for X0_i = 1:totalX0
    for out_i = 1:totalOut
        % Loading file and extract variables
        filename = "UncGMM_Test/Test" + Test_i + ...
                   "/" + def + "Def/X0_" + X0_i + ...
                    "/output_" + out_i + ".mat";
        load(filename, "n_iter", "Mu_found", "T_sim_toc", "Estim_sol",...
             "T_ErgC_i", "T_PDF_i", "Mu", "Sigma", "Sigma_found")
        
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
        % Defect Covariances
        Sigma_buffer{X0_i, out_i} = Sigma;
        Sigma_found_buffer{X0_i, out_i} = Sigma_found;
    end
end

%% Success rate (Number of tests in which all the defects are found)
successRate = sum(nbMu == def, 'all') / numel(nbMu) * 100;

% Success rate (porcentage of overall found defects)
% TotalFoundDef = 0;
% for X0_i = 1:totalX0
%     for out_i = 1:totalOut
%         TotalFoundDef = TotalFoundDef + height(Mu_found_buffer{X0_i, out_i});
%     end
% end
% 
% TotalDef = def*totalX0*totalOut;
%
% SR_def = TotalFoundDef/(TotalDef) * 100;

% SR Considering false positives
idx_fp = nbMu > def;
nbFalsePositives = zeros(size(nbMu));
nbFalsePositives(idx_fp) = nbMu(idx_fp) - def;
nbTotalFoundDef = nbMu;
nbTotalFoundDef(idx_fp) = def; 

TotalDef = def*totalX0*totalOut;
TotalFoundDef = sum(nbTotalFoundDef, "all");
FalsePositives = sum(nbFalsePositives, "all");

SR_def = TotalFoundDef / (TotalDef + FalsePositives) * 100;

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


%% Shape metrics

EllipseIoU = [];
OrientationError = [];

for X0_i = 1:totalX0
    for out_i = 1:totalOut
        % If found defects buffer its not empty
        if ~isempty(Mu_found_buffer{X0_i, out_i})
            k_idx = dsearchn(Mu_buffer{X0_i, out_i},...
                             Mu_found_buffer{X0_i, out_i});
            Mu_tmp = Mu_buffer{X0_i, out_i}(k_idx, :);
            Sigma_tmp = Sigma_buffer{X0_i, out_i}(:,:,k_idx);
            % Mu_tmp matches with current Mu_found_buffer{X0_i, out_i}
            % Sigma_tmp matches with current Sigma_found_buffer{X0_i, out_i}

            % Compute vector of Ellipse IoU metric on found defects
            EllipseIoU_current = getEllipseIoU(...
                                        Mu_tmp,...
                                        Sigma_tmp,...
                                        Mu_found_buffer{X0_i, out_i},...
                                        Sigma_found_buffer{X0_i, out_i},...
                                        100);
            EllipseIoU = [EllipseIoU;...
                          EllipseIoU_current];

            % Compute vector of Major axis orientation errors on found
            % defects
            OrientationErrorCurrent = getOrientationError(...
                                        Sigma_tmp,...
                                        Sigma_found_buffer{X0_i, out_i});
            OrientationError = [OrientationError; OrientationErrorCurrent];

        end
    end
end

% Mean and standard deviations
avgEIoU = mean(EllipseIoU);
stdEIoU = std(EllipseIoU);

avgOrientationError = mean(OrientationError);
stdOrientationError = std(OrientationError);


