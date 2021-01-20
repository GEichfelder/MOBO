function min_A_i = Solver(parameter)
%SOLVER Solves Mulitcriterial Bilevel Optimization Problems with the algorithm introduced by 
%G. Eichfelder in 2007.
% Inputs:
%   parameter:  string, name of a .m file where parameters for the algorithm are defined
% Outputs:
%   min_A_i:    Matrix of size ( . , m_2+n_1+n_2+m_1+m_1+p+2), structure of the columns:
%               F(x,y) | x | y | a | mu | nu | t | exitflag.

    %% +++++ supress MATLAB warnings, enable custom warnings +++++++++++++++++++++++++++++++++++++++
    warning('off', 'all')
    warning('backtrace', 'off')
    warning('on', 'Custom:Model:G_tilde')
    warning('on', 'Custom:dimensions')
    warning('on', 'Custom:step2a:minf_j:exitflag0')
    warning('on', 'Custom:step2a:minf_j:exitflag2')
    warning('on', 'Custom:step2a:SP:exitflag0')
    warning('on', 'Custom:step2b:exitflag0')
    warning('on', 'Custom:step2b:exitflag2')
    warning('on', 'Custom:step2b:Matrixsingular')
    warning('on', 'Custom:step2:A_0empty')
    warning('on', 'Custom:step4:exitflag0')
    warning('on', 'Custom:step4:exitflag2')
    warning('on', 'Custom:step4:Matrixsingular')

    %% +++++ make subfolders available +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    addpath(genpath('parameters'));
    addpath(genpath('models'));
    
    %% +++++ get the parameters ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    [c, d, beta, alpha, delta, nD, rep_step4, gamma, color, ...
        MaxIterations, MaxFunctionEvaluations, OptimalityTolerance] ...
        = feval(parameter);
    want_plot = ~isempty(color);
    global n_1 n_2 m_1 m_2 p p_tilde Model
    if n_2 < 1
        error('y must be at least one dimensional.')
    end

    %% +++++ get dimensions of the functions +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    set_function_dimensions();
    
    %% +++++ test if Model is valid ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if p_tilde < 1
        error('G_tilde must be at least one dimensional.')
    end
    % test, if dimensions are given correctly:
    check_dimensions();
    % test, if G_tilde is independent of y:
    check_G_tilde();
  
    %% +++++ some basic definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    E_m_1 = eye(m_1);
    E_m_1_n_2 = eye(m_1+n_2);
    Legend = {};  % legend for plot built up iteratively
    r = E_m_1(:, m_1);
    x0 = ones(n_1,1);
    options = optimoptions('fmincon','Display','none', ...
                           'MaxIterations', MaxIterations, ...
                           'MaxFunctionEvaluations', MaxFunctionEvaluations, ...
                           'OptimalityTolerance', OptimalityTolerance);

    %% +++++ step 1 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    k = zeros(n_2, 1);      % vector for values of k which the program is using to calculate the
                            % discretisation points of [c, d] (grid) 
                            % (y_k = c + k_1*beta_1*e_1 + ... + k_n_2*beta_n_2*e_n_2)
    k_max = ((d-c)./beta);  % vector holding maximal values for each entry of k
    A_0 = [];
    fprintf('step 1 finished \n')

    %% +++++ step 2 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    while k(n_2) <= k_max(n_2)
        y_k = c + k.*beta;
        % test, if y_k is feasible for G_tilde, necessary if G_tilde is not the whole box [c, d]
        if smaller_eq(feval(Model, zeros(n_1, 1), y_k, 'G_tilde'), zeros(p_tilde, 1))
            % min {f_{m_1} | (x, y_k) \in G}
            [x_m_1_bar, ~] = ...
                    fmincon(@(x)f(x, y_k, m_1), x0, [], [], [], [], [], [], ...
                            @(x)G(x, y_k), options);
            a_E = feval(Model, x_m_1_bar, y_k, 'f');
            j = 1;
            % for every direction in H repeat Step 2a and 2b
            while j <= m_1-1 || (j == 1 && m_1 == 1) 
                %% +++++ step 2a +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                [x_j_bar, ~, exitflag] = ...
                    fmincon(@(x)f(x, y_k, j), x0, [], [], [], [], [], [], ...
                            @(x)G(x, y_k), options);
                if exitflag > 0
                    f_x_j_bar = feval(Model, x_j_bar, y_k, 'f');
                    a = [f_x_j_bar(1:m_1-1); 0];
                    % determine mu, nu Lagrange multipliers of (SP(a, r, a_tilde)) with
                    % minimal solution f(x_bar_1, y_k, m_1) for step 4 of the algorithm
                    t0x0 = [f_x_j_bar(m_1); x_j_bar];
                    [tx, t, exitflag, ~, lambda] = ...
                        fmincon(@(tx)SPfun(tx), t0x0, [], [], [], [], [], [], ...
                                @(tx)SPcon(tx, a, r, y_k), options);
                    % exitflag == -2 can't happen because initial point is feasible
                    if exitflag == 0
                        warning('Custom:step2a:SP:exitflag0', ...
                                ['fmincon: MaxIterations or MaxFunctionEvaluations ', ...
                                 'exceeded in step 2a at solving SP(a, r, y^k) to get ' ...
                                 'the Lagrange Multipliers for x_%d_bar ', ...
                                 'minimal solution of min(f_%d)\n', ...
                                 'nevertheless using result as ', ...
                                 'first point for y^k = %s in A_0\n'], j, j, mat2str(y_k, 5))
                    end
                    % if (SP(a, r, a_tilde)) returns another x than x_j_bar:
                    x_j_bar = tx(2:1+n_1);
                    f_x_j_bar = feval(Model, x_j_bar, y_k, 'f');
                    a = [f_x_j_bar(1:m_1-1); 0];
                    mu = lambda.ineqnonlin(1:m_1);
                    nu = lambda.ineqnonlin(m_1+1:m_1+p);
                    A_0 =  [A_0; ...
                            feval(Model, x_j_bar, y_k, 'F').', x_j_bar.', y_k.', ...
                            a.', mu.', nu.', t, exitflag];   
                    % update a:
                    v_a = E_m_1(:, j); % Direction of variation of a depends on j
                    a = a + delta*v_a;
                    v_tilde = [v_a; zeros(m_1, 1)];
                    %% +++++ step 2b +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    while a(j) <= a_E(j)
                        % if m_1==1 this will not happen
                        % solve (SP(a, r, y_k)):
                        [tx, t, exitflag, ~, lambda] = ...
                            fmincon(@(tx)SPfun(tx), ones(1+n_1, 1), [], [], [], [], [], [], ...
                                    @(tx)SPcon(tx, a, r, y_k), options);
                        if exitflag > 0
                            x_l = tx(2:n_1+1);
                            mu_l = lambda.ineqnonlin(1:m_1);
                            nu_l = lambda.ineqnonlin(m_1+1:m_1+p);
                            A_0 = [A_0; ...
                                   feval(Model, x_l, y_k, 'F').', x_l.', y_k.', ...
                                   a.', mu_l.', nu_l.', t, exitflag];
                            % compute M, N as in Theorem 4.3:
                            M = [M_1(t, x_l, a, r, y_k, mu_l, nu_l), M_2(x_l, y_k)];
                            N = N_mat(t, mu_l);
                            % compute lambda_bar as in (5.2):
                            lastwarn('', '');
                            invMNv_tilde = M\(N*v_tilde);
                            [~, warnId] = lastwarn();
                            if(isempty(warnId))
                                lambda_bar = alpha / norm(invMNv_tilde(2:n_1+1));
                                % update a:
                                a = a + lambda_bar * v_a;
                            else
                                warning('Custom:step2b:Matrixsingular', ...
                                        ['(M)^{-1}*N is close to singular or badly scaled \n', ...
                                         'skipped update of a for y^k = %s\n'], mat2str(y_k, 5))
                                % break while loop
                                a(j) = a_E(j) + 1;
                            end
                        else
                            if exitflag == 0
                                warning('Custom:step2b:exitflag0', ...
                                        ['fmincon: MaxIterations or MaxFunctionEvaluations ', ...
                                         'exceeded in step 2b for y^k = %s ', ...
                                         'in direction %d, a_%d = %0.4f \n ', ...
                                         'left a out \n'], mat2str(y_k, 5), j, j, a(j))
                            elseif exitflag == -2
                                warning('Custom:step2b:exitflag2', ...
                                        ['fmincon: no feasible point was found ', ...
                                         'in step 2b for y^k = %s ', ...
                                         'in direction %d, a_%d = %0.4f \n ', ...
                                         'left a out \n'],  mat2str(y_k, 5), j, j, a(j))
                            end
                            % break while loop
                            a(j) = a_E(j) + 1;
                        end
                    end
                else
                    if exitflag == 0
                        warning('Custom:step2a:minf_j:exitflag0', ...
                                ['fmincon: MaxIterations or MaxFunctionEvaluations exceeded ', ...
                                 'in step 2a at solving min(f_%d) for y^k = %s \n', ...
                                 'left y^k out \n'], j, mat2str(y_k, 5))
                    elseif exitflag == -2
                        warning('Custom:step2a:minf_j:exitflag2', ...
                                ['fmincon: no feasible point was found ', ...
                                 'in step 2a at solving min(f_%d) for y^k = %s \n', ...
                                 'left y^k out \n'], j, mat2str(y_k, 5))
                    end
                end
                j = j + 1;
            end
        end
        % compute new k for grid of [c, d]:
        k(1) = k(1) + 1;
        for nextindex = 1:n_2
            if k(nextindex) > k_max(nextindex)
                k(nextindex) = 0;
                if nextindex < n_2
                    k(nextindex+1) = k(nextindex+1)+1;
                else    % k(n_2) reached maximum
                    k(n_2) = k_max(n_2)+1;      % break while-loop
                end
            end
        end
    end
    
    if isempty(A_0)
        error('Custom:step2:A_0empty', ['there was no point found for A^0 in step 2, ', ...
                                        'check your parameters\n'])
    end
    
    if want_plot
        Legend{1} = 'A^0';
        plot_f_y(A_0, color(1), 'x', Legend)
        plot_F(A_0, color(1), 'x', Legend)
    end

    fprintf('step 2 finished \n')

    %% +++++ step 3 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % i = 0:
    min_A_i = JGY(A_0);

    if want_plot
        Legend{2} = 'M(A^0)';
        plot_f_y(min_A_i, color(1), 'o', Legend)
        plot_F(min_A_i, color(1), 'o', Legend)
    end

    fprintf('step 3 finished \n')

    for i = 0:rep_step4-1
        size_min_A_i_old = size(min_A_i, 1);
        j = 1;
        % vary v_1 in H_hat -> variation of a_hat in the components wrt f
        while j <= m_1-1 || (j == 1 && m_1 == 1)   
            if m_1 > 1
                v_1 = E_m_1_n_2(:, j);      % direction of variation of a depends on j
            else 
                v_1 = zeros(m_1+n_2, 1);    % if m_1 == 1 do not adjust a_hat in the first component
            end
            for k = 1:n_2   % vary v_2 in H_hat -> variation of a_hat in the components wrt y
                v_2 = E_m_1_n_2(:, m_1+k);  % direction of variation of a depends on k
    %% +++++ step 4 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                for l = 1:size_min_A_i_old
                    % pick variables out of min_A_i:
                    t_0 = min_A_i(l, m_2+n_1+n_2+m_1+n_2+p+1);
                    x_0 = min_A_i(l, m_2+1:m_2+n_1).';
                    y_0 = min_A_i(l, m_2+n_1+1:m_2+n_1+n_2).';
                    a_0 = min_A_i(l, m_2+n_1+n_2+1:m_2+n_1+n_2+m_1).';
                    % Lagrange Multipliers to (SP(a_hat, r_hat)):
                    mu_0 = min_A_i(l, m_2+n_1+n_2+m_1+1:m_2+n_1+n_2+m_1+m_1).';
                    nu_0 = min_A_i(l, m_2+n_1+n_2+m_1+m_1+1:m_2+n_1+n_2+m_1+m_1+p).';
                    % compute mu_0_hat as in Thm 4.2:
                    mu_0_hat = mu_hat(mu_0, nu_0, x_0, y_0);
                    % compute M_hat and N_hat as in Thm 5.2 to compute s_i as in before (5.5) to 
                    % compute a_hat as in (5.5):
                    M_hat = M_hat_mat(t_0, x_0, y_0, a_0, r, mu_0_hat, nu_0);
                    N_hat = N_hat_mat(mu_0_hat);
                    % test, if matrix is close to singular or badly scaled
                    lastwarn('', '');
                    invM_hatN_hat = M_hat\N_hat;
                    [~, warnId] = lastwarn();
                    if(isempty(warnId))
                        if m_1 > 1
                            invM_hatN_hatv_1 = invM_hatN_hat*v_1;
                            s_1 = gamma(i+1) / norm(invM_hatN_hatv_1(2:1+n_1+n_2));
                        else
                            s_1 = 0;
                        end
                        invM_hatN_hatv_2 = invM_hatN_hat*v_2;
                        s_2 = gamma(i+1) / norm(invM_hatN_hatv_2(2:1+n_1+n_2));
                        for l_1 = -nD:nD
                            for l_2 = -nD:nD
                                if l_1 ~= 0 || l_2 ~= 0
                                    % compute a_hat = [a_hat_1; a_tilde] as in (5.5)
                                    a_hat_1 = a_0 + l_1*s_1*v_1(1:m_1);
                                    a_tilde  = y_0 + l_2*s_2*v_2(m_1+1:m_1+n_2);
                                    % test if a_tilde \in G_tilde, if not take the next point
                                    if smaller_eq(feval(Model, zeros(n_1, 1), a_tilde, 'G_tilde'), ...
                                                  zeros(p_tilde, 1))
                                        % solve (SP(a, r, a_tilde))
                                        tx0 = ones(1+n_1, 1);
                                        [tx, t, exitflag, ~, lambda] = ...
                                            fmincon(@(tx)SPfun(tx), tx0, [], [], [], [], [], [], ...
                                                    @(tx)SPcon(tx, a_hat_1, r, a_tilde), options);
                                        if exitflag > 0
                                            x = tx(2:1+n_1);          
                                            mu = lambda.ineqnonlin(1:m_1);
                                            nu = lambda.ineqnonlin(m_1+1:m_1+p);
                                            % A_i+1 as in eq after (5.6), already safe it as
                                            % min_A_i, because M(A_i) subset of A_i+1 -> no new
                                            % variable
                                            min_A_i = [min_A_i; ...
                                                feval(Model, x, a_tilde, 'F').', x.', a_tilde.', ...
                                                a_hat_1.', mu.', nu.', t, exitflag];
                                        elseif exitflag == 0
                                            warning('Custom:step4:exitflag0', ...
                                                    ['fmincon: MaxIterations or ', ...
                                                     'MaxFunctionEvaluations exceeded ', ...
                                                     'in step 4 for x = %s, ', ...
                                                     'y = %s, ', ...
                                                     'a_hat = %s \n', ...
                                                     'left a_hat out \n'], ...
                                                     mat2str(x_0.', 5), mat2str(y_0.', 5), ...
                                                     mat2str([a_hat_1; a_tilde].', 5))
                                        elseif exitflag == -2
                                            warning('Custom:step4:exitflag2', ...
                                                    ['fmincon: no feasible point was found ', ...
                                                     'in step 4 ', ...
                                                     'for x = %s, ', ...
                                                     'y = %s, ', ...
                                                     'a_hat = %s \n', ...
                                                     'left a_hat out \n'], ...
                                                     mat2str(x_0.', 5), mat2str(y_0.', 5), ...
                                                     mat2str([a_hat_1; a_tilde].', 5))
                                        end
                                    end
                                end
                            end
                        end
                    else 
                        warning('Custom:step4:Matrixsingular', ...
                                ['(M_hat)^{-1}*N_hat is close to singular or badly scaled \n', ...
                                 'left resulting a_hat out\n'])
                    end
                end
            end
            j = j + 1;
        end
        
        if want_plot
            Legend{2*i+3} = strcat('A^', num2str(i+1));
            plot_f_y(min_A_i, color(i+2), 'x', Legend)
            plot_F(min_A_i, color(i+2), 'x', Legend)
        end

        fprintf('step 4 finished for i=%d \n', i)
        
        if size_min_A_i_old == size(min_A_i, 1)
            fprintf(['no improvement in step 4 for i = %d, program stoppped\n', ...
                     'if only got Warning regarding (M_hat)^{-1}*N_hat then step 4 will never ', ...
                     'improve A^%d \n', ...
                     'if got Warning regarding fmincon in step 4, then try a smaller gamma\n'], ...
                     i, i);
            break;
        end
    %% +++++ step 5 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % min_A_i+1 as in (5.7), already safe it as min_A_i
        min_A_i = JGY(min_A_i);
        
        if want_plot
            Legend{2*i+4} = strcat('M(A^', num2str(i+1), ')');
            plot_f_y(min_A_i, color(i+2), 'o', Legend)
            plot_F(min_A_i, color(i+2), 'o', Legend)
        end
        
        fprintf('step 5 finished for i=%d \n', i)
    end
    
    fprintf('program finished \n')

end

%% +++++++++++++++++++ subfunctions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function f_val = f(x, y, i)
%F compute value of each f function (lower level)
% Inputs:
%   x       :   n_1 dimensional vector            
%   y       :   n_2 dimensional vector
%   i       :   int optional index which component of f should be returned
% Outputs:
%   f_val   :   m_1 dimensional vector.
    global Model
    f_val = feval(Model, x, y, 'f');
    f_val = f_val(i);
end

function [c,ceq] = G(x,y)
%G compute nonlinear constraint for fmincon (x, y) \in G
% Inputs:
%   x       :   n_1 dimensional vector            
%   y       :   n_2 dimensional vector
% Outputs:
%   c       :   p dimensional vector, holds values of g
%   ceq     :   empty, return parameter for fmincon.
    global Model
    c = feval(Model, x, y, 'G');
    ceq = [];
end

function t = SPfun(tx)
%SPFUN extract first entry of tx.
    t = tx(1);
end

function [c, ceq] = SPcon(tx, a, r, a_tilde)
%SPCON Compute nonlinear constraints for fmincon, see (SP(a, r, a_tilde))
% Inputs:
%   tx      :   1+n_1 dimensional vector, holding (t, x) of current point
%   a       :   m_1 dimensional vector
%   r       :   m_1 dimensional vector
%   a_tilde :   n_2 dimensional vector, element of G_tilde
% Outputs:
%   c       :   m_1+p dimensional vector of nonlinear inequality 
%               constraints at (t, x)
%   ceq     :   0 dimensional vector of nonlinear equality constraints 
%               at (t, x).
    global n_1 m_1 p Model
    t = tx(1);
    x = tx(2:1+n_1);
    c(1:m_1) = -(a + t * r - feval(Model, x, a_tilde, 'f'));
    c(m_1+1:m_1+p) = feval(Model, x, a_tilde, 'G');
    ceq = [];
end

function M_1_val = M_1(t, x, a, r, a_tilde, mu, nu)
%M_1 Compute M_1 as in Theorem 4.3
% Inputs:
%   t       :   scalar
%   x       :   n_1 dimensional vector
%   a       :   m_1 dimensional vector
%   r       :   m_1 dimensional vector
%   a_tilde :   n_2 dimensional vector
%   mu      :   m_1 Lagrange multipliers to the constraint 
%               a+tr-f(x, a_tilde)
%   nu      :   p Lagrange multipliers to the constraint (x, a_tilde) \in G
% Outputs:
%   M_1_val :   (1+n_1+m_1+p)x(1+n_1+m_1) dimensional Matrix, see Thm 4.3.
    global n_1 m_1 p Model
    k0 = a + t*r - feval(Model, x, a_tilde, 'f');
    M_1_val = zeros(1+n_1+m_1+p, 1+n_1+m_1);
    M_1_val(1, 1+n_1+1:1+n_1+m_1) = -r;
    M_1_val(1+n_1+1:1+n_1+m_1, 1) = mu .* r;      
    M_1_val(2:1+n_1, 2:1+n_1) = d2dxL(x, a_tilde, mu, nu);   
    dxf_val(:, 1:n_1) = feval(Model, x, a_tilde, 'f', 'x');
    for i = 1:m_1
        M_1_val(2:1+n_1, 1+n_1+i) = dxf_val(i, :).';
        M_1_val(1+n_1+i, 2:1+n_1) = -mu(i)*dxf_val(i, :);
        M_1_val(1+n_1+i, 1+n_1+i) = k0(i);
    end
    dxG_val = feval(Model, x, a_tilde, 'G', 'x');
    for i = 1:p
        M_1_val(1+n_1+m_1+i, 2:1+n_1) = nu(i)*dxG_val(i, :);
    end
end

function d2dxL_val = d2dxL(x, a_tilde, mu, nu)
%D2DXL Compute Hessian of Lagrange function of (SP(a, r, a_tilde)) wrt x
% Inputs:
%   x           :   n_1 dimensional vector
%   a_tilde     :   n_2 dimensional vector
%   mu          :   m_1 Lagrange multipliers to the constraint 
%                   a+tr-f(x, a_tilde)
%   nu          :   p Lagrange multipliers to the constraint 
%                   (x, a_tilde) \in G
% Outputs:
%   d2dxL_val   :   n_1 x n_1 dimensional matrix.
    global n_1 m_1 p Model
    d2dxL_val = zeros(n_1, n_1);
    d2dxf_val = feval(Model, x, a_tilde, 'f', 'xx');
    for i = 1:m_1
        d2dxL_val = d2dxL_val + ...
            mu(i)*d2dxf_val((i-1)*n_1+1:(i-1)*n_1+n_1, :);
    end
    d2dxG_val = feval(Model, x, a_tilde, 'G', 'xx');
    for i = 1:p
        d2dxL_val = d2dxL_val + ...
            nu(i)*d2dxG_val((i-1)*n_1+1:(i-1)*n_1+n_1, :);
    end
end

function d2L_val = d2L(x, y, mu, nu)
%D2L Compute Hessian of Lagrange function of (SP(a_hat, r_hat))
% Inputs:
%   x           :   n_1 dimensional vector
%   y           :   n_2 dimensional vector
%   mu          :   m_1+n_2 Lagrange multipliers to the constraints
%                   a__hat+tr_hat-f_hat(x, y)
%   nu          :   p Lagrange multipliers to the constraint 
%                   (x, y) \in G
% Outputs:
%   d2dxL_val   :   (1+n_1+n_2) x (1+n_1+n_2) dimensional matrix.
    global n_1 n_2 m_1 p Model
    d2L_val = zeros(1+n_1+n_2);
    d2dxf = feval(Model, x, y, 'f', 'xx');
    d2dxyf = feval(Model, x, y, 'f', 'xy');
    d2dyf = feval(Model, x, y, 'f', 'yy');
    for i = 1:m_1
        d2L_val(1+1:1+n_1, 1+1:1+n_1) = ...
            d2L_val(1+1:1+n_1, 1+1:1+n_1) + mu(i)*d2dxf((i-1)*n_1+1:(i-1)*n_1+n_1, :);
        d2L_val(1+1:1+n_1, 1+n_1+1:1+n_1+n_2) = ...
            d2L_val(1+1:1+n_1, 1+n_1+1:1+n_1+n_2) + mu(i)*d2dxyf((i-1)*n_2+1:(i-1)*n_2+n_2, :).';
        d2L_val(1+n_1+1:1+n_1+n_2, 1+1:1+n_1) = ...
            d2L_val(1+n_1+1:1+n_1+n_2, 1+1:1+n_1) + mu(i)*d2dxyf((i-1)*n_2+1:(i-1)*n_2+n_2, :);
        d2L_val(1+n_1+1:1+n_1+n_2, 1+n_1+1:1+n_1+n_2) = ...
            d2L_val(1+n_1+1:1+n_1+n_2, 1+n_1+1:1+n_1+n_2) ...
            + mu(i)*d2dyf((i-1)*n_2+1:(i-1)*n_2+n_2, :);
    end
    d2dxG = feval(Model, x, y, 'G', 'xx');
    d2dxyG = feval(Model, x, y, 'G', 'xy');
    d2dyG = feval(Model, x, y, 'G', 'yy');
    for i = 1:p
        d2L_val(1+1:1+n_1, 1+1:1+n_1) = ...
            d2L_val(1+1:1+n_1, 1+1:1+n_1) + nu(i)*d2dxG((i-1)*n_1+1:i*n_1, :);
        d2L_val(1+1:1+n_1, 1+n_1+1:1+n_1+n_2) = ...
            d2L_val(1+1:1+n_1, 1+n_1+1:1+n_1+n_2) + nu(i)*d2dxyG((i-1)*n_2+1:i*n_2, :).';
        d2L_val(1+n_1+1:1+n_1+n_2, 1+1:1+n_1) = ...
            d2L_val(1+n_1+1:1+n_1+n_2, 1+1:1+n_1) + nu(i)*d2dxyG((i-1)*n_2+1:i*n_2, :);
        d2L_val(1+n_1+1:1+n_1+n_2, 1+n_1+1:1+n_1+n_2) = ...
            d2L_val(1+n_1+1:1+n_1+n_2, 1+n_1+1:1+n_1+n_2) ...
            + nu(i)* d2dyG((i-1)*n_2+1:i*n_2, :);
    end
end

function M_2_val = M_2(x, a_tilde)
%M_2 Compute M_2 as in Theorem 4.3
% Inputs:
%   x       :   n_1 dimensional vector
%   a_tilde :   n_2 dimensional vector
% Outputs:
%   M_2_val :   (1+n_1+m_1+p) x p dimensional Matrix, see Thm 4.3.
    global n_1 m_1 p Model
    M_2_val = zeros(1+n_1+m_1+p,p);
    dxG_val = feval(Model, x, a_tilde, 'G', 'x');
    for i = 1:p
        M_2_val(2:1+n_1, i) = - dxG_val(i, :);
    end
    G_val = feval(Model, x, a_tilde, 'G');
    for i = 1:p
        M_2_val(1+n_1+m_1+i, i) = G_val(i);
    end
end

function N_val = N_mat(t, mu)
%N_mat Compute N as in Theorem 4.3
% Inputs:
%   t       :   scalar
%   mu      :   m_1 Lagrange multipliers to the constraint 
%               a+tr-f(x, a_tilde)
% Outputs:
%   N_val   :   (1+n_1+m_1+p) x (2*m_1) dimensional Matrix, see Thm 4.3.
    global n_1 m_1 p
    N_val = zeros(1+n_1+m_1+p, 2*m_1);
    E = eye(m_1);
    for i = 1:m_1
        N_val(1+n_1+i, :) = -mu(i)*[E(i, :), t*E(i, :)];
    end
end

function mu_hat_val = mu_hat(mu, nu, x, a_tilde)
%MU_HAT Compute Lagrange Multiplier to (SP(a_hat, r_hat)) as in Thm 4.2
% Inputs:
%   mu          :   m_1 Lagrange multipliers to the constraint 
%                   a+tr-f(x, a_tilde) in (SP(a, r, a_tilde))
%   nu          :   p Lagrange multipliers to the constraint
%                   (x, a_tilde) \in G in (SP(a, r, a_tilde))
%   x           :   n_1 diemnsional vector
%   a_tilde     :   n_2 dimensional vector
% Outputs:
%   mu_hat_val  :   m_1+n_2 Lagrange multipliers to the constraint
%                   a_hat+tr_hat-f_hat(x, y) in (SP(a_hat, r_hat).
global n_2 m_1 p Model
    mu_tilde = zeros(n_2, 1);
    dyf_val = feval(Model, x, a_tilde, 'f', 'y');
    for i = 1:m_1
        mu_tilde = mu_tilde - mu(i)*dyf_val(i, :).';
    end
    dyG_val = feval(Model, x, a_tilde, 'G', 'y');
    for i = 1:p
        mu_tilde = mu_tilde + nu(i)*dyG_val(i,:).';
    end
    mu_hat_val = [mu; mu_tilde];
end

function M_hat_val = M_hat_mat(t, x, y, a, r, mu_hat, nu)
%M_HAT_MAT Compute M_hat as in Theorem 5.2
% Inputs:
%   t           :   scalar
%   x           :   n_1 dimensional vector
%   y           :   n_2 dimensional vector
%   a           :   m_1 dimensional vector, first components of a_hat
%   mu_hat      :   m_1+n_2 Lagrange multipliers  to the constraints
%                   a_hat+t*r_hat-f_hat(x, y)
%   nu          :   p Lagrange multipliers to the constraint 
%                   (x, y) \in G
% Outputs:
%   M_hat_val   :   (1+n_1+n_2+m_1+n_2+p) x (1+n_1+n_2+m_1+n_2+p) dimensional 
%                   Matrix, see Thm 5.2.
    global n_1 n_2 m_1 p Model
    M_hat_val = zeros(1+n_1+n_2+m_1+n_2+p);                                                     
    M_hat_val(1:1+n_1+n_2, 1:1+n_1+n_2) = d2L(x, y, mu_hat, nu);
    M_hat_val(1, 1+n_1+n_2+m_1) = -1;
    M_hat_val(1+n_1+n_2+m_1, 1) = mu_hat(m_1);
    df_val = zeros(m_1, n_1+n_2);    
    df_val(:, 1:n_1) = feval(Model, x, y, 'f', 'x');
    df_val(:, n_1+1:n_1+n_2) = feval(Model, x, y, 'f', 'y');
    k = a + t*r - feval(Model, x, y, 'f');
    for i = 1:m_1
        M_hat_val(1+n_1+n_2+i, 1+1:1+n_1+n_2) = -mu_hat(i)*df_val(i, :);
        M_hat_val(1+1:1+n_1+n_2, 1+n_1+n_2+i) = df_val(i, :).';
        M_hat_val(1+n_1+n_2+i, 1+n_1+n_2+i) = k(i);
    end
    for i = 1:n_2
        M_hat_val(1+n_1+n_2+m_1+i, 1+n_1+i) = -mu_hat(m_1+i);
        M_hat_val(1+n_1+i, 1+n_1+n_2+m_1+i) = 1;       
    end
    dG_val = zeros(p, n_1+n_2);   
    dG_val(:, 1:n_1) = feval(Model, x, y, 'G', 'x');
    dG_val(:, n_1+1:n_1+n_2) = feval(Model, x, y, 'G', 'y');
    G_val = feval(Model, x, y, 'G');
    for i = 1:p
        M_hat_val(1+n_1+n_2+m_1+n_2+i, 1+1:1+n_1+n_2) = nu(i)*dG_val(i, :);
        M_hat_val(1+1:1+n_1+n_2, 1+n_1+n_2+m_1+n_2+i) = -dG_val(i, :).';
        M_hat_val(1+n_1+n_2+m_1+n_2+i, 1+n_1+n_2+m_1+n_2+i) = G_val(i);
    end
end

function N_hat_val = N_hat_mat(mu_hat)
%N_HAT_MAT Compute M_hat as in Theorem 5.2
% Inputs:
%   mu_hat      :   m_1+n_1 Lagrange multipliers  to the constraints
%                   a_hat+t*r_hat-f_hat(x, y)
% Outputs:
%   N_hat_val   :   (1+n_1+n_2+m_1+n_2+p) x (m_1+n_2) dimensional Matrix.
    global n_1 n_2 m_1 p
    N_hat_val = zeros(1+n_1+n_2+m_1+n_2+p, m_1+n_2);
    for i = 1:m_1+n_2
        N_hat_val(1+n_1+n_2+i, i) = -mu_hat(i);
    end
end

function B = JGY(A)
%JGY Jahn-Graf-Younes modified Algorithm, see "New algorithms for discrete
%vector optimization based on Greaf-Younes method and cone-monotone sorting
%functions" (by C. Günther & N. Popvici), Algorithm 4 with 
%phi(a_i) = (1, ..., 1)^T * a_i.
    global m_2
    k = size(A, 1);
    % Phase 1:
    phi_A = sum(A(:, 1:m_2), 2);
    [~, index_sorted_phi_A] = sort(phi_A);
    % Phase 2:
    B = A(index_sorted_phi_A(1), :);
    for i = 2:k
        attach = true;
        for j = 1:size(B, 1)
            if bigger_neq(A(index_sorted_phi_A(i), 1:m_2), B(j, 1:m_2))
                attach = false;
                break
            end
        end
        if attach
            B = [B; ...
                 A(index_sorted_phi_A(i), :)];
        end
    end
end

function test_bigger_than_ref = bigger_neq(test_val, ref_val)
%BIGGER_THAN test if test_val \in ref_val + K\setminus{0}
% Inputs:
%   test_val                :   n dimensional vector
%   ref_val                 :   n dimensional vector
% Outputs:
%   test_bigger_than_ref    :   bool, true if test_val ~= ref_val and test_val >= ref_val wrt K
%                               (see fct smallereq_than).
    test_bigger_than_ref = smaller_eq(ref_val, test_val);
    if test_bigger_than_ref && isequal(test_val, ref_val)
        test_bigger_than_ref = false;
    end
end

function test_smaller_than_ref = smaller_eq(test_val, ref_val)
%BIGGER_THAN test if test_val \in ref_val - R^n_+
% Inputs:
%   test_val                :   n dimensional vector
%   ref_val                 :   n dimensional vector
% Outputs:
%   test_smaller_than_ref   :   bool, true if test_val_i <= ref_val_i for 
%                               all i\in{1,...,n}.
% change function if you use another cone
    test_smaller_than_ref = true;
    for i = 1:length(test_val)
        if test_val(i) > ref_val(i)
            test_smaller_than_ref = false;
            break
        end
    end
end

function plot_f_y (set, color, marker, Legend)
%PLOT_F_Y if possible plot f(x,y) or f(x,y) and y
%Inputs: 
%   set     :   Matrix holding in rows x, y at column position m_2+1:m_2+n_1 resp.
%               m_2+n_1+1:m_2+n_1+n_2
%   color   :   color for the plot
%   marker  :   marker for the plot
%   Legend  :   legend for the plot.
    global n_1 n_2 m_1 m_2 Model
    figure(1);
    hold on
    if m_1 == 2 && n_2 == 1
        f = zeros(size(set, 1), m_1);
        for l = 1:size(set, 1)
            f(l, :) = feval(Model, set(l, m_2+1:m_2+n_1).', set(l, m_2+n_1+1), 'f');
        end
        plot3(f(:, 1),...
              f(:, 2), ...
              set(:, m_2+n_1+1), 'color', color, 'marker', marker, 'LineStyle', 'none');
        xlabel('f_1')
        ylabel('f_2')
        zlabel('y')
    elseif m_1 == 1 && n_2 == 2
        f = zeros(size(set, 1), m_1);
        for l = 1:size(set, 1)
            f(l, :) = feval(Model, set(l, m_2+1:m_2+n_1).', set(l, m_2+n_1+1:m_2+n_1+2).', 'f');
        end
        plot3(f(:, 1),...
              set(:, m_2+n_1+1), ...
              set(:, m_2+n_1+2), 'color', color, 'marker', marker, 'LineStyle', 'none');
        xlabel('f_1')
        ylabel('y_1')
        zlabel('y_2')
    elseif m_1 == 1 && n_2 == 1
        f = zeros(size(set, 1), 1);
        for l = 1:size(set, 1)
            f(l, :) = feval(Model, set(l, m_2+1:m_2+n_1).', set(l, m_2+n_1+1), 'f');
        end
        plot(f, set(:, m_2+n_1+1), 'color', color, 'marker', marker, 'LineStyle', 'none');
        xlabel('f')
        ylabel('y')
    elseif m_1 == 1
        f = zeros(size(set, 1), 1);
        for l = 1:size(set, 1)
            f(l, :) = feval(Model, set(l, m_2+1:m_2+n_1).', set(l, m_2+n_1+1:m_2+n_1+n_2), 'f');
        end
        plot(f(:, 1),  zeros(size(set,1), 1), 'color', color, 'marker', marker, 'LineStyle', 'none');
        xlabel('f_1')
    elseif m_1 == 2
        f = zeros(size(set, 1), m_1);
        for l = 1:size(set, 1)
            f(l, :) = feval(Model, set(l, m_2+1:m_2+n_1).', set(l, m_2+n_1+1:m_2+n_1+n_2), 'f');
        end
        plot(f(:, 1), f(:, 2), 'color', color, 'marker', marker, 'LineStyle', 'none');
        xlabel('f_1')
        ylabel('f_2')
    elseif m_1 == 3
        f = zeros(size(set, 1), m_1);
        for l = 1:size(set, 1)
            f(l, :) = feval(Model, set(l, m_2+1:m_2+n_1).', set(l, m_2+n_1+1:m_2+n_1+n_2), 'f');
        end
        plot3(f(:, 1), f(:, 2), f(:, 3), 'color', color, 'marker', marker, 'LineStyle', 'none');
        xlabel('f_1')
        ylabel('f_2')
        zlabel('f_3')
    else
        fprintf('not able to plot f or f and y, because m_1 > 3 or m_1+n_2 > 3\n');
    end
    axis equal
    legend(Legend)
    hold off
end

function plot_F (set, color, marker, Legend)
%PLOT_F if possible plot F(x, y)
%Inputs: 
%   set     :   Matrix holding in rows F(x,y) at column position 1:m_2
%   color   :   color for the plot
%   marker  :   marker for the plot
%   Legend  :   legend for the plot.
    global m_2
    figure(2);
    hold on
    if m_2 == 1
        plot(set(:, 1), zeros(size(set,1), 1), 'color', color, 'marker', marker, 'LineStyle', 'none');
        xlabel('F')
    elseif m_2 == 2
        plot(set(:, 1), set(:, 2), 'color', color, 'marker', marker, 'LineStyle', 'none');
        xlabel('F_1')
        ylabel('F_2')
    elseif m_2 == 3
        plot3(set(:, 1), set(:, 2), set(:,3), 'color', color, 'marker', marker, 'LineStyle', 'none');
        xlabel('F_1')
        ylabel('F_2')
        zlabel('F_3')
    else
        fprintf('not able to plot F, because m_2 > 3\n');
    end
    axis equal
    legend(Legend)
    legend boxon 
    hold off
end

function check_G_tilde()
%CHECK_G_TILDE check if function G_tilde is independent of x.
    global Model n_1 n_2
    G_test_1 = feval(Model, zeros(n_1, 1), ones(n_2, 1), 'G_tilde'); 
    G_test_2 = feval(Model, ones(n_1, 1), ones(n_2, 1), 'G_tilde');
    if ~isequal(G_test_1, G_test_2)
        error('Custom:Model:G_tilde', 'Model is invalid because G_tilde depends on x.\n')
    end
end

function set_function_dimensions()
    global Model n_1 n_2 m_1 m_2 p p_tilde
    x = ones(n_1, 1);
    y = ones(n_2, 1);
    m_1 = size(feval(Model, x, y, 'f'), 1);
    m_2 = size(feval(Model, x, y, 'F'), 1);
    p = size(feval(Model, x, y, 'G'), 1);
    p_tilde = size(feval(Model, x, y, 'G_tilde'), 1);
    if p_tilde < 1
        error('G_tilde(y) must be at least one dimensional')
    end
end

function check_dimensions()
%CHECK_DIMENSIONS check if dimensions of the derivatives of f and G are correct.
    global Model n_1 n_2 m_1 p
    x = ones(n_1, 1);
    y = ones(n_2, 1);
    lastwarn('','');
    keyf =    {'f', 'G'};
    dimkeyf = [ m_1, p];
    keyxy =    {'x', 'y', 'xx', 'xy', 'yy'};
    dimkeyxy = [ 1,   1,   n_1,  n_2,  n_2; ...
                 n_1, n_2, n_1,  n_1,  n_2];
    for i = 1:2       
        for j = 1:5
            size_derivative = size(feval(Model, x, y, keyf{i}, keyxy{j}));
            if ~isequal(size_derivative, [dimkeyf(i)*dimkeyxy(1, j) , dimkeyxy(2, j)])
                warning('Custom:dimensions', ...
                        'Dimension of d%s/d%s is inaccurate, should be %d x %d but is %d x %d', ...
                        keyf{i}, keyxy{j}, dimkeyf(i)*dimkeyxy(1, j), dimkeyxy(2, j), ...
                        size_derivative(1), size_derivative(2))
            end
        end
    end
    [~, warnId] = lastwarn();
    if ~isempty(warnId)
        error('Custom:dimensions', 'Dimensions doesn''t fit, see warning(s) above')
    end
end
