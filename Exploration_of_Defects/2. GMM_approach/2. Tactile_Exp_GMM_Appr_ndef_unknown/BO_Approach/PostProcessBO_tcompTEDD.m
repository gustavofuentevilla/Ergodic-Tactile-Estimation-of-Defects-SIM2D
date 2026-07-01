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

for X0_i = 18:18
    for out_i = 5:5
        
        % Load TEDD iteration time
        filename = "/home/gustavo-fuentevilla/MATLAB/" + ...
                   "Ergodic_Tactile_Estimation_of_Defects_SIM2D/" + ...
                   "Exploration_of_Defects/2. GMM_approach/" + ...
                   "2. Tactile_Exp_GMM_Appr_ndef_unknown/" + ...
                   "Test" + Test_i + "/Case" + Case_i +...
                   "/" + def + "Def/X0_" + X0_i + ...
                   "/output_" + out_i + ".mat";
        load(filename, "T_ErgC_i", "T_PDF_i")

        TEDD_time = sum(T_ErgC_i + T_PDF_i);

        % Load BO data
        filename = "BO_unconsTest/Test" + Test_i + "/Case" + Case_i +...
                   "/" + def + "Def/X0_" + X0_i + ...
                    "/output_" + out_i + ".mat";
        load(filename)
        
        % Compute the number of BO samples achieved when BO simulation time
        %  is larger than TEDD simulation time
        for BO_samples = 1:samples
            if sum(timer_iter(1:BO_samples)) > TEDD_time
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
        % filename = "BO_unconsTest/Test" + Test_i + "_new/Case" + Case_i +...
        %            "/" + def + "Def/X0_" + X0_i + ...
        %             "/output_" + out_i + ".mat";
        % 
        % save(filename, "-regexp",...
        %     "^(?!(Test_i|Case_i|def|totalX0|totalOut|X0_i|out_i)$).");


    end
end




