function [c, d, beta, alpha, delta, nD, rep_step4, gamma, color, ...
          MaxIterations, MaxFunctionEvaluations, OptimalityTolerance] = parameter_dummy()
global n_1 n_2
%% +++++ dimensions of the variables +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
n_1 = 0;        % dimension of x
n_2 = 0;        % dimension of y

%% +++++ parameter for the algorithm +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% parameters for the box containing G_tilde (see explanation after Algorithm 1):
c = [];                 % n_2 dimensional vector, lower bound of G_tilde
d = [];                 % n_2 dimensional vector, upper bound of G_tilde
% discretization parameter for G_tilde:
steps = 10;             % with how many steps G_tilde should be discretized (can be n_2 dimensional 
                        % or scalar)
beta = (d-c)./steps;    % n_2 dimensional vector, stepsize in every direction

% parameters for step 2:
% distance for nearly equidistant approximation of the set of feasible points of the 
% upper level problem (see (5.1)):
alpha = 0.1; 
% small number for slight mutation of the initial a for numerical reasons (see step 2a and 
% explanation below Algorithm 1):
delta = 0.1;

% parameters for step 4:
rep_step4 = 2;              % how often step 4 should be repeated (A^i refined)
nD = 2;                     % number of desired new discretization points (see (5.5))
% vector of length rep_step4 defining how large the variation of a_hat in step 4 should be (see
% before (5.5)):
gamma = [0.3, 0.15];  
% vector of length rep_step4 + 1 containing colors for plotting A^i after step 4 
% resp. M(A^i) after step 5, if empty no plot will be produced:
color = ['b', 'r', 'g'];   

% parameters for fmincon:
MaxIterations = 1000;
MaxFunctionEvaluations = 3000;
OptimalityTolerance = 1e-06;
end