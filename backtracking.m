function [alpha, x_new, f_new, accepted, n_evals] = ...
        backtracking(f, x, f0, g, d, alpha0, beta, c)
%BACKTRACKING  Backtracking line search with the Armijo condition.
%
%   Starting from alpha = alpha0, the step length is multiplied by beta
%   until the Armijo sufficient-decrease condition
%
%       f(x + alpha*d) <= f(x) + c * alpha * g'*d
%
%   holds, or until the maximum number of backtracks is reached.
%
%   Inputs
%       f       function handle of the objective
%       x       current iterate
%       f0      objective value f(x) (passed in to avoid re-evaluating it)
%       g       gradient at x
%       d       search direction
%       alpha0  initial step length
%       beta    contraction factor, 0 < beta < 1
%       c       Armijo parameter, 0 < c < 1
%
%   Outputs
%       alpha     accepted step length (0 if the search failed)
%       x_new     accepted point (x if the search failed)
%       f_new     objective value at x_new
%       accepted  true if an acceptable step was found
%       n_evals   number of objective evaluations used by the search

    max_backtracks = 50;

    alpha    = alpha0;
    x_new    = x;
    f_new    = f0;
    accepted = false;
    n_evals  = 0;

    % Backtracking is only meaningful for a descent direction.
    slope = g' * d;
    if slope >= 0
        alpha = 0;
        return;
    end

    for i = 1:max_backtracks
        x_trial = x + alpha * d;
        f_trial = f(x_trial);
        n_evals = n_evals + 1;

        % Armijo sufficient-decrease condition
        if f_trial <= f0 + c * alpha * slope
            x_new    = x_trial;
            f_new    = f_trial;
            accepted = true;
            return;
        end

        alpha = beta * alpha;
    end

    % No acceptable step was found.
    alpha = 0;
end
