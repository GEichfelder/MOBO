function [c, d, beta, alpha, delta, nD, rep_step4, gamma, color, ...
          MaxIterations, MaxFunctionEvaluations, OptimalityTolerance] = parameter_paper_example_2()
global n_1 n_2

global V 
load ('V.mat','V')

%% +++++ dimensions of the variables +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
n_1 = 14;        % dimension of x
n_2 = 1;        % dimension of y

%% +++++ parameter for the algorithm +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% parameters for the box containing G_tilde (see explanation after Algorithm 1):
c = [0.0001];                 % n_2 dimensional vector, lower bound of G_tilde
d = [pi];                 % n_2 dimensional vector, upper bound of G_tilde
% discretization parameter for G_tilde:
steps = 0;             % with how many steps G_tilde should be discretized (can be n_2 dimensional 
                        % or scalar)
beta = pi/8;    % n_2 dimensional vector, stepsize in every direction

% parameters for step 2:
% distance for nearly equidistant approximation of the set of feasible points of the 
% upper level problem (see (5.1)):
alpha = 0.2; 
% small number for slight mutation of the initial a for numerical reasons (see step 2a and 
% explanation below Algorithm 1):delta = 0.1;

% parameters for step 4:
rep_step4 = 3;              % how often step 4 should be repeated (A^i refined)
nD = 3;                     % number of desired new discretization points (see (5.5))
% vector of length rep_step4 defining how large the variation of a_hat in step 4 should be (see
% before (5.5)):
gamma = [1/35, 1/50, 1/70];
% vector of length rep_step4 + 1 containing colors for plotting A^i after step 4 
% resp. M(A^i) after step 5, if empty no plot will be produced:
color = ['b', 'r' ,'g', 'y'];

% parameters for fmincon:
MaxIterations = 5000;
MaxFunctionEvaluations = 10000;
OptimalityTolerance = 1e-06;
end