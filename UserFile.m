%% +++++ Cleanup +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clc;
clear;
close all;

%% +++++ Please enter your functions and parameters ++++++++++++++++++++++++++++++++++++++++++++++++
global Model n_1 n_2 m_1 m_2 p

% name of file holding functions:
Model = 'model_example';
% is it a model from BOLIBv2 ?
BOLIB = false;
% name of file holding parameters for the algorithm:
parameter = 'parameter_example';

%% +++++ No changes below here +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
global ModelBOLIB

if BOLIB
    ModelBOLIB = Model;
    Model = 'transformBOLIB';
end

% call solver
min_set = Solver(parameter);

% print results
fprintf('\n')
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('   F values        |    x values       |       y values     | exitflag\n')
for i = 1:size(min_set, 1)
    fprintf([repmat('%0.4f ', 1, m_2), '|'] , min_set(i, 1:m_2));
    fprintf([repmat('%0.4f ', 1, n_1), '|'], min_set(i, m_2+1:m_2+n_1));
    fprintf([repmat('%0.4f ', 1, n_2), '|'], min_set(i, m_2+n_1+1:m_2+n_1+n_2));
    fprintf('%d \n', min_set(i, m_2+n_1+n_2+m_1+m_1+p+2));
end
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
