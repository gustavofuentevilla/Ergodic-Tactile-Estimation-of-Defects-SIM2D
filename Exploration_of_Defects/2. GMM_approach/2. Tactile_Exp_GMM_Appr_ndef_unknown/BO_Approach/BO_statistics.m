close all
clear
clc

%%
% Test (2: Same severity levels, 3: Different severity levels)}
Test_i = 2;
% Case (1: full algoritm, 2: No repeating points, 3: No removing gaussians)
Case_i = 1;
% Number of Defects
def = 3;
% Number of initial positions
totalX0 = 20;
% Number of tests per initial position
totalOut = 5;

% Buffers for:
nbSamples = zeros(totalX0,totalOut);   % Number of iterations
nbMu = zeros(totalX0,totalOut);     % Number of found defects
T_PDF_buffer = zeros(totalX0,totalOut);
T_iter = cell(totalX0,totalOut);    % Times per iteration
Mu_buffer = cell(totalX0,totalOut);
Mu_found_buffer = cell(totalX0,totalOut);
Sigma_buffer = cell(totalX0,totalOut);
Sigma_found_buffer = cell(totalX0,totalOut);

for X0_i = 1:totalX0
    for out_i = 1:totalOut
        % Loading file and extract variables
        filename = "BO_unconsTest/Test" + Test_i + ...
                   "/Case1/" + def + "Def/X0_" + X0_i + ...
                    "/output_" + out_i + ".mat";
        load(filename, "samples", "Mu_found", "timer_iter", "Estim_sol",...
             "Mu", "Sigma", "Sigma_found", "timer_PDF", "BO_samples")
        
        % Save data into matrix (row = initial position; column = test)
        % nbSamples(X0_i, out_i) = samples;  % For BO_stdTD 
        nbSamples(X0_i, out_i) = BO_samples; % For the rest of BO
        nbMu(X0_i, out_i) = height(Mu_found);

        % Iteration times
        T_iter{X0_i, out_i} = timer_iter;
        T_PDF_buffer(X0_i, out_i) = timer_PDF; 

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

% Calculate the average values and standard deviations of total samples
avgSamples = mean(nbSamples, "all");
stdSamples = std(nbSamples, 0, "all");

% Simulation time per iteration
TimePerIter = [];
for X0_i = 1:totalX0
    for out_i = 1:totalOut
        TimePerIter = [TimePerIter; 
                       T_iter{X0_i, out_i}(1:nbSamples(X0_i, out_i)) + ...
                       T_PDF_buffer(X0_i, out_i)/nbSamples(X0_i, out_i)];
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





