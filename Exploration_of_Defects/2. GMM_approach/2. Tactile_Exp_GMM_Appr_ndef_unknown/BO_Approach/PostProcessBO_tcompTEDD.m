close all
clear
clc

%% 

% Test (2: Same severity levels, 3: Different severity levels)}
Test_i = 3;
% Case (1: full algoritm (TEDD), 2: No repeating points, 3: No removing gaussians)
Case_i = 1;
% Number of Defectss
def = 7;
% Number of initial positions
totalX0 = 20;
% Number of tests per initial position
totalOut = 5;

for X0_i = 1:totalX0
    for out_i = 1:totalOut
        
        % Load TEDD iteration time
        filename = "/home/gustavo-fuentevilla/MATLAB/" + ...
                   "Ergodic_Tactile_Estimation_of_Defects_SIM2D/" + ...
                   "Exploration_of_Defects/2. GMM_approach/" + ...
                   "2. Tactile_Exp_GMM_Appr_ndef_unknown/" + ...
                   "Test" + Test_i + "/Case" + Case_i +...
                   "/" + def + "Def/X0_" + X0_i + ...
                   "/output_" + out_i + ".mat";
        load(filename, "T_ErgC_i", "T_PDF_i", "n_iter", "X_e_spline_reg")

        Traveled_distance = 0;
        for j = 1:n_iter
            P_diff = diff(X_e_spline_reg(:,:,j));
            Traveled_distance = Traveled_distance + ...
                                sum( sqrt( P_diff(:,1).^2 + P_diff(:,2).^2 ) );
        
        end
        TEDD_time = sum(T_ErgC_i + T_PDF_i);
        Tact_meas = n_iter * size(X_e_spline_reg, 1);

        % Load BO data
        filename = "BO_unconsTest/Test" + Test_i + "/Case" + Case_i +...
                   "/" + def + "Def/X0_" + X0_i + ...
                    "/output_" + out_i + ".mat";
        load(filename)
        
        % Compute the number of BO samples achieved when BO simulation time
        %  is larger than TEDD simulation time
        % for BO_samples = 1:samples
        %     if sum(timer_iter(1:BO_samples)) > TEDD_time
        %         break;
        %     end
        % end

        % Compute the number of BO samples achieved when BO traveled
        % distance is larger than TEDD's
        % BO_trav_dist = 0;
        % for BO_samples = 1:samples
        %     BO_P_diff = diff(X_e(BO_samples:BO_samples+1, :));
        %     BO_trav_dist = BO_trav_dist + ...
        %                    sum(sqrt(BO_P_diff(:,1).^2 + BO_P_diff(:,2).^2));
        %     if BO_trav_dist > Traveled_distance
        %         break;
        %     end
        % end

        % Compute the number of BO samples when the number of tactile
        % measurements are larger than TEDD's

        for BO_samples = 1:samples
            [~, idx_init] = ismember(X_e(BO_samples,:), X_e_sp, "rows");
            [~, idx_final] = ismember(X_e(BO_samples+1,:), X_e_sp, "rows");
            X_extract = X_e_sp(idx_init:idx_final, :);
            BOtactileMeasurements = size(X_extract, 1);
            if BOtactileMeasurements >= Tact_meas
                break;
            end
        end
        
        % Remove posterior data to that sample if BO takes more time than
        % TEDD
        if BO_samples < samples
            timer_iter(BO_samples+1:end) = 0;
            X_e(BO_samples+2:end,:) = [];
            t_exec(BO_samples+2:end) = [];
            loglikelihood(BO_samples+1:end) = [];
    
            [~, idx_Xesp] = ismember(X_e(end,:), X_e_sp, "rows");
            X_e_sp(idx_Xesp:end, :) = [];
            V(idx_Xesp:end, :) = [];
        end

        % Post-processing

        timer_PDF_init = tic;
        
        Par_PDF.thres_meas = a + max(delta);
        Par_PDF.D_max = 10;
        Par_PDF.Thres_Variation = max(sum(r_elips_Phi));
        Par_PDF.OneClustDistLimit = 2*max(r_elips_Phi,[],"all") + 0.07;
        
        Estim_sol = BO_Postprocessing(X_e_sp, V, Par_PDF);
        
        Priors_found = Estim_sol.Priors_found;
        Mu_found = Estim_sol.Mu_found;
        Sigma_found = Estim_sol.Sigma_found_a;
        
        timer_PDF = toc(timer_PDF_init);


        % Save new simulations
        filename = "BO_unconsTest/Test" + Test_i + "_tm/Case" + Case_i +...
                   "/" + def + "Def/X0_" + X0_i + ...
                    "/output_" + out_i + ".mat";

        save(filename, "-regexp",...
            "^(?!(Test_i|Case_i|def|totalX0|totalOut|X0_i|out_i)$).");


    end
end




