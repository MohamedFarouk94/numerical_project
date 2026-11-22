function [x, history] = damped_newton(x0, params)
%DAMPED_NEWTON  Newton's method with Armijo backtracking (damped Newton).
%
%   [x, history] = DAMPED_NEWTON(x0, params) minimizes the Rosenbrock
%   function starting from x0. The search direction solves
%
%       H(x) d = -grad f(x).
%
%   If d is not a descent direction (or the linear system is singular),
%   the method falls back to d = -grad f(x) and counts the event in
%   history.fallbacks.
%
%   params and history are as described in STEEPEST_DESCENT.

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
    history.fallbacks = 0;
    status = 'max_iter';

    for k = 1:params.max_iter

        if norm(g) <= params.epsilon
            status = 'converged';
            break;
        end

        % Newton direction: solve H*d = -g (no explicit inverse)
        H = hessian(x);
        history.hevals = history.hevals + 1;
        d = H \ (-g);

        % Fall back to steepest descent if d is not a descent direction
        if any(~isfinite(d)) || g' * d >= 0
            d = -g;
            history.fallbacks = history.fallbacks + 1;
        end

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
