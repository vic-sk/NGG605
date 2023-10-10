% LATER Fitting Exercise
%
% Copyright 2023 by Joshua I. Gold, University of Pennsylvania

%edits by Victoria Subritzky Katz, 10/10/2023

% The basic idea in fitting a model to data is to find the parameters of
% the model that provide in some sense the best match of the model to the
% data. This match is provided by the "objective function." 
% This exercise is intended to demystify this process by getting you
% to define the initial conditions and objective function for 
% fitting the LATER model to RT data. For a much more thorough, but still 
% very accessible, overview of model fitting (to behavioral data), here
% is a great place to start:
%
% https://elifesciences.org/articles/49547
%
% For this exercise, recall that the point of the LATER model is that 1/RT is
% distributed as a Gaussian, where we can define the parameters
% of the Gaussian (mu and sigma) with respect to the standard parameters
% of the LATER model (muR and deltaS):
%       mu = muR/deltaS
%       sigma = 1/deltaS
%
% So fitting LATER to behavioral data involves finding parameters
% muR and deltaS that provide the best match to the data, according to 
% the appropriate objective function.
%
% Follow along the steps below, some of which will require you to complete
% the code (and therefore hopefully think about how to relate the high-
% level concepts discussed above with the nitty-gritty part of getting 
% everything to actually work.

%%  1. Get the data
%   
%   Use this code to get a data set (array of RTs from a single condition) 
%   to fit, already preprocessed to include correct trials only and remove
%   outliers (including express saccades). See later_getData for details
data = later_getData([], [], 0.2);
RTs = data{1};
clear data
RTs_rec = 1./RTs; 
RTs_rec_neg = -RTs_rec;

%%  2. Define the objective function
%
% The objective function typically defines the error that you want to 
% minimize between your data and the model predictions. A common objective 
% function is the negative of the sum of the log-likelihoods of the data, 
% given the model parameters. To unpack that for the LATER model:
%
%   1. For each data point (RT from a single trial, in this case) and given
%       set of model parameters, compute the probability of the data, given
%       the model (i.e., the likelihood)
%   2. Take the logarithm
%   3. Sum all these log-likelihoods from all the data points
%   4. Take the negative, because we want to find the minimum (thus
%        corresponding to the maximum likelihood)
%
%   You can define the function simply using an "anonymous function"
%   (https://www.mathworks.com/help/matlab/matlab_prog/anonymous-functions.html), 
%   using this template that assumes that "fits" is a 2x1 vector of
%   [muR, deltaS]:
 
% EXERCISE:
%laterErrFcn = @(params, RTs_rec_neg) neg_log_liklihood(params, RTs_rec_neg);
%laterErrFcn = @(params, RTs) neg_log_liklihood(params, RTs_rec_neg);
laterErrFcn = @(params, RTs) neg_log_liklihood(params, RTs_rec);
%laterErrFcn =  @(fits) -sum(log(normpdf(RTs_rec_neg, fits(1)./fits(2), 1./fits(2))));
%laterErrFcn =  @(fits) -sum(log(normpdf(RTs_rec, fits(1)./fits(2), 1./fits(2))));


%%  3. Define initial conditions
%   
%   For the actual fitting, we will use fmincon
%   (https://www.mathworks.com/help/optim/ug/fmincon.html), which is 
%   "function minimization with constraints." This function allows for 
%   constraints that include upper and lower bounds on the parameters.
%   So here we define those bounds, along with the initial values.
%   We'll use fairly arbitrary values for the lower and upper
%   bounds, but we should pick the initial values more judiciously. HINT: 
%   Recall that the muR and deltaS should be strongly related to 
%   empirical summary statistics of `the (reciprocal) RT distribution.
lowerBounds = [0.001 0.001];
upperBounds = [1000 1000]; 

% EXERCISE:
% mu_init = mean(RTs_rec_neg);
% sigma_init = std(RTs_rec_neg); 
mu_init = mean(RTs_rec);
sigma_init = std(RTs_rec); 
initialValues = [mu_init,sigma_init];

%%  4. Run the fits
% 
%   We will be using GlobalSearch . The general advantage of this approach 
%   is to avoid local minima; for details, see:
%   https://www.mathworks.com/help/gads/how-globalsearch-and-multistart-work.html
%  
%   These options seem to work well, but I don't have a stronger
%   rationale for using them. See the Matlab documentation if you really
%   want to dive in and understand them, and let me know if you find
%   better settings!
opts = optimoptions(@fmincon,    ... % "function minimization with constraints"
   'Algorithm',   'active-set',  ...
   'MaxIter',     3000,          ...
   'MaxFunEvals', 3000);

% Definine the "optimization problem" using variables defined above
problem = createOptimProblem('fmincon',    ...
    'objective',   laterErrFcn,     ... % Use the objective function
    'x0',          initialValues,   ... % Initial conditions
    'lb',          lowerBounds,     ... % Parameter lower bounds
    'ub',          upperBounds,     ... % Parameter upper bounds
    'options',     opts);                % Options defined above

% Create a GlobalSearch object
gs = GlobalSearch;
   
% Run it, returning the best-fitting parameter values and the negative-
% log-likelihood returned by the objective function
%[fits(ii,:), nllk] = run(gs,problem);
[fits(:), nllk] = run(gs,problem);
%%  5. Evaluate the fits
% How do you know if you got a reasonable answer?

%visual the fit of the model on the data
laterModelFit_plot(RTs, fits)
deltaS = 1/fits(2); 
muR = fits(1)*deltaS; 

%plotting the model on the data, both in the form of a histogram of
%1/(reaction times) allows for eye check of whether the fit is reasonable,
%judging by how well it matches the data. For a better eveluation of
%the quality of the fit you could calculate error bars and other statistics


function ngll = neg_log_liklihood(params, data)
    % params(1) is the mean and params(2) is the standard deviation
    mu = params(1);
    sigma = params(2);
    
    % Compute the negative log likelihood of each data point under the normal distribution
    ngll = -sum(log(normpdf(data, mu, sigma)));
end