function [Priors, Mu, Sigma, Pix] = GMM_inc_EM(Data, K, Prev_Priors, Prev_Mu, Prev_Sigma, Prev_Pix, iteration, Opts)
%
% Incremental Expectation-Maximization estimation of GMM parameters.
% This source code is the implementation of the algorithms described in 
% Section 2.6.1, p.47 of the book "Robot Programming by Demonstration: A 
% Probabilistic Approach".
%
% Author:	Sylvain Calinon, 2009
%			http://programming-by-demonstration.org
%
% This function learns the parameters of a Gaussian Mixture Model 
% (GMM) using a recursive Expectation-Maximization (EM) algorithm, starting 
% from an initial kmeans estimation of the parameters.
%
%
% Inputs -----------------------------------------------------------------
%   o Data:     N x D array representing N datapoints of D dimensions.
%   o K:        Number of GMM components.
%   o Prev_Priors: 1 x K Prior probabilities from last iteration (feedback)
%   o Prev_Mu:  K x D Means from last iteration for initializing the next
%               (feedback)
%   o Prev_Sigma:  D x D x K Covariance matrix from last iteration for
%                   initializing the next (feedback)
%   o Prev_Pix: N x K Posterior probability from last iteration (feedback)
%   o iteration: 1 x 1 Iteration index
%   o Opts:     1 x 1 Struct of options
% Outputs ----------------------------------------------------------------
%   o Priors:   1 x K array representing the prior probabilities of the
%               K GMM components.
%   o Mu:       K x D array representing the centers of the K GMM components.
%   o Sigma:    D x D x K array representing the covariance matrices of the 
%               K GMM components.
%   o Pix:      N x K Posterior probability of having the i-th gaussian 
%               given the data p(i|x)
% Comments ---------------------------------------------------------------
%   o This function uses standard EM algorithm in the first iteration and
%   with that solution as a feedback, computes the next iteration solutions
%   with the incremental EM algorithm. No historical Data is used, only the
%   current Data and the previous solution of the algorithm.
%
% This source code is given for free! However, I would be grateful if you refer 
% to the book (or corresponding article) in any academic publication that uses 
% this code or part of it. Here are the corresponding BibTex references: 
%
% @book{Calinon09book,
%   author="S. Calinon",
%   title="Robot Programming by Demonstration: A Probabilistic Approach",
%   publisher="EPFL/CRC Press",
%   year="2009",
%   note="EPFL Press ISBN 978-2-940222-31-5, CRC Press ISBN 978-1-4398-0867-2"
% }
%
% @article{Calinon07,
%   title="On Learning, Representing and Generalizing a Task in a Humanoid Robot",
%   author="S. Calinon and F. Guenter and A. Billard",
%   journal="IEEE Transactions on Systems, Man and Cybernetics, Part B",
%   year="2007",
%   volume="37",
%   number="2",
%   pages="286--298",
% }

% for n dimensional Data --> Data (:,n)
arguments
    Data (:,2) {mustBeNumeric}
    K    (1,1) {mustBeInteger,mustBePositive}
    Prev_Priors
    Prev_Mu
    Prev_Sigma
    Prev_Pix
    iteration (1,1) {mustBeInteger,mustBePositive}
    Opts.Max_iter (1,1) {mustBeInteger,mustBePositive} = 100
    Opts.loglik_threshold (1,1) {mustBeNumeric,mustBePositive} = 1e-10
    Opts.Min_var (1,1) {mustBeNumeric,mustBePositive} = 1e-5
end

%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%
loglik_old = -realmax;
nbStep = 0;
[nbData, nbDim] = size(Data);

%%%%%%%%%%%%%%%%%%%% Primera iteración %%%%%%%%%%%%%%%%%%%%%%%
if iteration == 1 
    
    %Standard GMM with k-means initialization
    [Priors, Mu, Sigma, Pix] = GMM_EM(Data, K, "Max_iter", Opts.Max_iter,...
        "loglik_threshold", Opts.loglik_threshold, "Min_var", Opts.Min_var); 

%%%%%%%%%%%%%%%%%%%% Resto de iteraciones %%%%%%%%%%%%%%%%%%%%%%%
else   
    % Initialization of Incremental EM: direct method
    Priors0 = Prev_Priors;
    Mu0 = Prev_Mu;
    Sigma0 = Prev_Sigma;
    Pix0 = Prev_Pix;
    nbData0 = size(Pix0, 1);
    Mu = Mu0;
    Sigma = Sigma0;
    Priors = Priors0;
    % Compute cumulated posterior probability of previous iteration
    E0 = sum(Pix0);
    
    % Incremental EM update algorithm
    while 1

      % E-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Pxi = zeros(nbData, K);
      for i = 1:K
        %Compute probability p(x|i)
        Pxi(:,i) = gaussPDF(Data, Mu(i,:), Sigma(:,:,i));
      end
      Pxi = Pxi + realmin; %To avoid zeros (and NaN's in next calculations)
      
      Pix_tmp = repmat(Priors, [nbData 1]).*Pxi; 
      Pix = Pix_tmp ./ repmat(sum(Pix_tmp,2), [1 K]);

      E = sum(Pix, 1);

      % M-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for i = 1:K
        %Update the priors
        Priors(i) = (E0(i) + E(i)) / (nbData0 + nbData);
        %Update the centers
        Mu(i,:) = (Mu0(i,:).*E0(i) + (Pix(:,i)')*Data) / (E0(i) + E(i));
        %Update the covariance matrices
        covtmp = zeros(nbDim, nbDim);
        for j = 1:nbData
            covtmp = covtmp + (Data(j,:) - Mu(i,:))'*(Data(j,:) - Mu(i,:)).*Pix(j,i);
        end
        Sigma(:,:,i) = ((Sigma0(:,:,i) + (Mu0(i,:) - Mu(i,:))*...
                (Mu0(i,:) - Mu(i,:))').*E0(i) + covtmp) / (E0(i) + E(i));
      end
      % Add a tiny variance to avoid numerical instability
      Sigma = Sigma + eye(nbDim)*Opts.Min_var;

      % Stopping criterion %%%%%%%%%%%%%%%%%%%%
      for i = 1:K
        %Compute the new probability p(x|i)
        Pxi(:,i) = gaussPDF(Data, Mu(i,:), Sigma(:,:,i));
      end

      % Compute the log likelihood
      F = Pxi*Priors';
      F(F < realmin) = realmin;
      loglik = mean(log(F));

      % Stop the process depending on the increase of the log likelihood
      nbStep = nbStep + 1;
      if abs((loglik/loglik_old) - 1) < Opts.loglik_threshold
          disp("Iterations in E-M algorithm: "+ nbStep)
          break;
      end
      loglik_old = loglik;
      if nbStep >= Opts.Max_iter
          warning("Maximum number of iterations in E-M algorithm reached")
          break;
      end
      
    end
    
    % Update posterior probability
    Pix = [Pix0; Pix];

end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function prob = gaussPDF(Data, Mu, Sigma)
%
% This function computes the Probability Density Function (PDF) of a
% multivariate Gaussian represented by means and covariance matrix.
%
% Author:	Sylvain Calinon, 2009
%			http://programming-by-demonstration.org
%
% Inputs -----------------------------------------------------------------
%   o Data:  N x D array representing N datapoints of D dimensions.
%   o Mu:    K x D array representing the centers of the K GMM components.
%   o Sigma: D x D x K array representing the covariance matrices of the 
%            K GMM components.
% Outputs ----------------------------------------------------------------
%   o prob:  N x 1 array representing the probabilities for the 
%            N datapoints.     

[nbData, nbDim] = size(Data);

Data_Mu = Data - repmat(Mu, nbData, 1); % (Data - Mu)
prob = sum( (Data_Mu/(Sigma)).*Data_Mu, 2 );
prob = exp(-0.5*prob) / sqrt( (2*pi)^nbDim * (abs(det(Sigma)) + realmin) );

end

function [Priors, Mu, Sigma, Pix] = GMM_EM(Data, K, Opts)

% Inputs -----------------------------------------------------------------
%   o Data:     N x D array representing N datapoints of D dimensions.
%   o K:        Number of GMM components.
% Outputs ----------------------------------------------------------------
%   o Priors:   1 x K array representing the prior probabilities of the
%               K GMM components.
%   o Mu:       K x D array representing the centers of the K GMM components.
%   o Sigma:    D x D x K array representing the covariance matrices of the 
%               K GMM components.
%   o Pix:      N x K Posterior probability of having the i-th gaussian 
%               given the data p(i|x)
% Comments ---------------------------------------------------------------
%   o This function uses the 'kmeans' function from the MATLAB Statistics 
%     toolbox. If you are using a version of the 'netlab' toolbox that also
%     uses a function named 'kmeans', please rename the netlab function to
%     'kmeans_netlab.m' to avoid conflicts. 

% for n dimensional Data --> Data (:,n)
arguments
    Data (:,2) {mustBeNumeric}
    K    (1,1) {mustBeInteger,mustBePositive}
    Opts.Max_iter (1,1) {mustBeInteger,mustBePositive} = 100
    Opts.loglik_threshold (1,1) {mustBeNumeric,mustBePositive} = 1e-10
    Opts.Min_var (1,1) {mustBeNumeric,mustBePositive} = 1e-5
end

%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%
loglik_old = -realmax;
nbStep = 0;
[nbData, nbDim] = size(Data);

%Use of the 'kmeans' function from the MATLAB Statistics toolbox
[Data_id, Centers] = kmeans(Data, K, "Replicates", 5); 
Mu0 = Centers;
Priors0 = zeros(1, K);
Sigma0 = zeros(nbDim, nbDim, K);
for i = 1:K
  idtmp = find(Data_id == i);
  Priors0(i) = length(idtmp); %Número de datos que contiene el iésimo grupo Priors = [datos_grp1, datos_grp2, ...]
  Sigma0(:,:,i) = cov([Data(idtmp,:); Data(idtmp,:)]); %Se concatena 2 veces para obtener menor varianza
end
% Add a tiny variance to avoid numerical instability
Sigma0 = Sigma0 + eye(nbDim)*Opts.Min_var;
Priors0 = Priors0 / sum(Priors0); % Porcentaje de datos en cada grupo, Pi

% First guesses
Mu = Mu0; %K x D array
Sigma = Sigma0; %D x D x K array
Priors = Priors0; %1 x K array

% Standard EM algorithm
while 1

  % E-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Pxi = zeros(nbData, K);
  for i = 1:K
    %Compute probability p(x|i)
    Pxi(:,i) = gaussPDF(Data, Mu(i,:), Sigma(:,:,i)); %Probabilidad de que se presenten los puntos x por el gaussiano i
  end
  Pxi = Pxi + realmin; %To avoid zeros (and NaN's in posterior calculations)
  % Compute posterior probability p(i|x): bayes rule
  Pix_tmp = repmat(Priors, [nbData 1]).*Pxi; 
  Pix = Pix_tmp ./ repmat(sum(Pix_tmp,2), [1 K]);  %(nbData, K) Probabilidad de que un x (dado) pertenezca al gaussiano i
  % Compute cumulated posterior probability
  E = sum(Pix); %(1,K) Probabilidad total sobre cada gaussiano

  % M-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for i = 1:K
    % Update the priors
    Priors(i) = E(i) / nbData; %Probabilidad de que un punto x (sin observar) pertenezca al iésimo gaussiano
    % Update the centers
    Mu(i,:) = (Pix(:,i)')*Data / E(i);
    % Update the covariance matrices
    Data_Mu = Data - repmat(Mu(i,:), [nbData 1]);
    Sigma(:,:,i) = (repmat(Pix(:,i)', [nbDim, 1]) .* (Data_Mu')*Data_Mu) / E(i);
  end
  % Add a tiny variance to avoid numerical instability
  Sigma = Sigma + eye(nbDim)*Opts.Min_var;

  % Stopping criterion %%%%%%%%%%%%%%%%%%%%
  for i = 1:K
    % Compute the new probability p(x|i)
    Pxi(:,i) = gaussPDF(Data, Mu(i,:), Sigma(:,:,i));
  end
  % Compute the log likelihood
  F = Pxi*Priors';
  F(F < realmin) = realmin;
  loglik = mean(log(F));
  % Stop the process depending on the increase of the log likelihood 
  nbStep = nbStep + 1;
  if abs((loglik/loglik_old) - 1) < Opts.loglik_threshold
    disp("Iterations in E-M algorithm: "+ nbStep)
    break;
  end
  loglik_old = loglik;
  
  if nbStep >= Opts.Max_iter
      warning("Maximum number of iterations in E-M algorithm reached")
      break;
  end

end

end

