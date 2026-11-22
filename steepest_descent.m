function [x, history] = steepest_descent(x0, params)
%STEEPEST_DESCENT  Steepest Descent with Armijo backtracking.
%
%   [x, history] = STEEPEST_DESCENT(x0, params) minimizes the Rosenbrock
%   function starting from x0, using the search direction d = -grad f(x).
%
%   params is a struct with the fields
%       epsilon   gradient-norm tolerance
%       max_iter  maximum number of iterations
%       alpha0    initial step length of the line search
%       beta      backtracking contraction factor
%       c         Armijo parameter
%
%   history contains the iterates (x, one column per iterate), the
%   objective values (f), gradient norms (grad_norm), accepted step
%   lengths (alpha), line-search evaluations per iteration (ls_evals),
%   the evaluation counters (fevals, gevals, hevals), and the
%   termination status ('converged', 'max_iter' or 'line_search_failed').

    f = @rosenbrock;

    x   = x0(:);
    f_x = f(x);
    g   = grad(x);

    history.x         = x;
    history.f         = f_x;
    history.grad_norm = norm(g);
    history.alpha     = [];
    history.ls_evals  = [];
    history.gevals    = 1;      % gradient at x0
    history.hevals    = 0;
    status = 'max_iter';

    for k = 1:params.max_iter

        if norm(g) <= params.epsilon
            status = 'converged';
            break;
        end

        % Steepest Descent direction
        d = -g;

        % Backtracking line search
        [alpha, x_new, f_new, accepted, n_ls] = backtracking( ...
            f, x, f_x, g, d, params.alpha0, params.beta, params.c);
        history.ls_evals(end+1) = n_ls;

        if ~accepted
            status = 'line_search_failed';
            break;
        end

        % Update the iterate
        x   = x_new;
        f_x = f_new;
        g   = grad(x);
        history.gevals = history.gevals + 1;

        % Store results
        history.alpha(end+1)     = alpha;
        history.x(:, end+1)      = x;
        history.f(end+1)         = f_x;
        history.grad_norm(end+1) = norm(g);
    end

    if strcmp(status, 'max_iter') && norm(g) <= params.epsilon
        status = 'converged';
    end

    history.status = status;
    history.fevals = 1 + sum(history.ls_evals);   % f(x0) + line-search trials
end
