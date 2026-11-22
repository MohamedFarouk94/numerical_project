function [x, history] = bfgs(x0, params)
%BFGS  BFGS quasi-Newton method with Armijo backtracking.
%
%   [x, history] = BFGS(x0, params) minimizes the Rosenbrock function
%   starting from x0. B approximates the INVERSE Hessian, starts at the
%   identity, and is updated with
%
%       B <- (I - rho*s*y') * B * (I - rho*y*s') + rho*s*s',
%       rho = 1 / (y'*s),   s = x_new - x,   y = grad(x_new) - grad(x).
%
%   The update is skipped when the curvature condition y'*s > 0 fails
%   (Armijo backtracking does not guarantee it), and B is reset to the
%   identity if -B*g is not a descent direction. Both events are counted
%   in history.skipped_updates and history.resets.
%
%   params and history are as described in STEEPEST_DESCENT.

    curvature_tol = 1e-10;      % relative threshold for y'*s > 0

    f = @rosenbrock;

    x   = x0(:);
    n   = length(x);
    I   = eye(n);
    B   = I;
    f_x = f(x);
    g   = grad(x);

    history.x               = x;
    history.f               = f_x;
    history.grad_norm       = norm(g);
    history.alpha           = [];
    history.ls_evals        = [];
    history.gevals          = 1;    % gradient at x0
    history.hevals          = 0;
    history.skipped_updates = 0;
    history.resets          = 0;
    status = 'max_iter';

    for k = 1:params.max_iter

        if norm(g) <= params.epsilon
            status = 'converged';
            break;
        end

        % BFGS search direction
        d = -B * g;

        % Reset to steepest descent if d is not a descent direction
        if g' * d >= 0
            B = I;
            d = -g;
            history.resets = history.resets + 1;
        end

        % Backtracking line search
        [alpha, x_new, f_new, accepted, n_ls] = backtracking( ...
            f, x, f_x, g, d, params.alpha0, params.beta, params.c);
        history.ls_evals(end+1) = n_ls;

        if ~accepted
            status = 'line_search_failed';
            break;
        end

        % Changes in the iterate and in the gradient
        s     = x_new - x;
        g_new = grad(x_new);
        history.gevals = history.gevals + 1;
        y     = g_new - g;
        ys    = y' * s;

        % BFGS inverse-Hessian update (only if y'*s > 0)
        if ys > curvature_tol * norm(y) * norm(s)
            rho = 1 / ys;
            B   = (I - rho * s * y') * B * (I - rho * y * s') ...
                  + rho * (s * s');
            B   = (B + B') / 2;     % remove numerical asymmetry
        else
            history.skipped_updates = history.skipped_updates + 1;
        end

        % Update the iterate
        x   = x_new;
        f_x = f_new;
        g   = g_new;

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
