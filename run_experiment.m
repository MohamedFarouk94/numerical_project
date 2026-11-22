function results = run_experiment(name, x0, params, n_runs, time_budget)
%RUN_EXPERIMENT  Run the three methods from x0, time them, and print a table.
%
%   results = RUN_EXPERIMENT(name, x0, params, n_runs, time_budget)
%
%   Each method is run once to collect its history and metrics, and is then
%   timed over repeated independent runs (a warm-up run is included in the
%   first call). Timing stops after n_runs runs or after time_budget
%   seconds, whichever comes first, but never before 20 runs. The timed
%   runs include the cost of recording the history.
%
%   results is a 1x3 struct array (Steepest Descent, Damped Newton, BFGS)
%   with the fields: method, x, history, iterations, f, dist, grad_norm,
%   fevals, gevals, hevals, status, n_runs, time_median, time_mean,
%   time_std and time_per_iter (all times in milliseconds).

    if nargin < 4, n_runs = 1000; end
    if nargin < 5, time_budget = 5; end

    min_runs = 20;
    x_star   = [1; 1];
    methods  = {'Steepest Descent', 'Damped Newton', 'BFGS'};
    solvers  = {@steepest_descent, @damped_newton, @bfgs};

    results = struct('method', {}, 'x', {}, 'history', {}, ...
        'iterations', {}, 'f', {}, 'dist', {}, 'grad_norm', {}, ...
        'fevals', {}, 'gevals', {}, 'hevals', {}, 'status', {}, ...
        'n_runs', {}, 'time_median', {}, 'time_mean', {}, ...
        'time_std', {}, 'time_per_iter', {});

    for m = 1:numel(methods)

        % --- Metrics from one reference run --------------------------
        [x, h] = solvers{m}(x0, params);

        % --- Timing over repeated independent runs -------------------
        times   = zeros(1, n_runs);
        n_done  = 0;
        elapsed = 0;                                  % seconds
        while n_done < n_runs && (elapsed < time_budget || n_done < min_runs)
            t_id = tic;
            solvers{m}(x0, params);
            t = toc(t_id);
            n_done = n_done + 1;
            times(n_done) = 1000 * t;                 % milliseconds
            elapsed = elapsed + t;
        end
        times = times(1:n_done);

        r.method        = methods{m};
        r.x             = x;
        r.history       = h;
        r.iterations    = numel(h.alpha);
        r.f             = h.f(end);
        r.dist          = norm(x - x_star);
        r.grad_norm     = h.grad_norm(end);
        r.fevals        = h.fevals;
        r.gevals        = h.gevals;
        r.hevals        = h.hevals;
        r.status        = h.status;
        r.n_runs        = n_done;
        r.time_median   = median(times);
        r.time_mean     = mean(times);
        r.time_std      = std(times);
        r.time_per_iter = median(times) / max(r.iterations, 1);

        results(m) = r;
    end

    print_report(name, x0, params, results);
end


%% ------------------------------------------------------------------------
function print_report(name, x0, params, res)
% Print a human-readable table followed by ready-to-paste LaTeX rows.

    line = repmat('-', 1, 78);

    fprintf('\n%s\n', repmat('=', 1, 78));
    fprintf('EXPERIMENT: %s   (x0 = [%g, %g], tol = %g)\n', ...
        name, x0(1), x0(2), params.epsilon);
    fprintf('%s\n', repmat('=', 1, 78));
    fprintf('%-26s %16s %16s %16s\n', 'Metric', res.method);
    fprintf('%s\n', line);

    print_row('Status',               '%16s',    {res.status});
    print_row('Iterations',           '%16d',    {res.iterations});
    print_row('Final objective',      '%16.6e',  {res.f});
    print_row('Distance to minimum',  '%16.6e',  {res.dist});
    print_row('Final gradient norm',   '%16.6e', {res.grad_norm});
    print_row('Function evaluations', '%16d',    {res.fevals});
    print_row('Gradient evaluations', '%16d',    {res.gevals});
    print_row('Hessian evaluations',  '%16d',    {res.hevals});
    print_row('Median time (ms)',     '%16.4f',  {res.time_median});
    print_row('Std of time (ms)',     '%16.4f',  {res.time_std});
    print_row('Median time/iter (ms)', '%16.4f', {res.time_per_iter});
    print_row('Timed runs',           '%16d',    {res.n_runs});
    fprintf('%s\n', line);

    h_b = res(3).history;
    h_n = res(2).history;
    fprintf('BFGS: %d skipped update(s), %d reset(s).  ', ...
        h_b.skipped_updates, h_b.resets);
    fprintf('Damped Newton: %d steepest-descent fallback(s).\n', ...
        h_n.fallbacks);

    fprintf('\nFinal points:\n');
    for m = 1:numel(res)
        fprintf('  %-17s (%.8f, %.8f)\n', [res(m).method ':'], ...
            res(m).x(1), res(m).x(2));
    end

    fprintf('\nLaTeX rows (Metric & SD & Damped Newton & BFGS):\n');
    print_tex('Iterations',                  {res.iterations}, 'int');
    print_tex('Final objective value',       {res.f},          'sci');
    print_tex('Distance to minimum',         {res.dist},       'sci');
    print_tex('Final gradient norm',         {res.grad_norm},  'sci');
    print_tex('Function evaluations',        {res.fevals},     'int');
    print_tex('Gradient evaluations',        {res.gevals},     'int');
    print_tex('Hessian evaluations',         {res.hevals},     'int');
    print_tex('Median time (ms)',            {res.time_median}, 'fix');
    print_tex('Median time per iteration (ms)', {res.time_per_iter}, 'fix');
    fprintf('\n');
end


function print_row(label, fmt, values)
    fprintf('%-26s', label);
    for k = 1:numel(values)
        fprintf([' ' fmt], values{k});
    end
    fprintf('\n');
end


function print_tex(label, values, kind)
    fprintf('%s', label);
    for k = 1:numel(values)
        fprintf(' & %s', tex_number(values{k}, kind));
    end
    fprintf(' \\\\\n');
end


function s = tex_number(v, kind)
% Format a number for LaTeX: integer, fixed point, or scientific notation.
    switch kind
        case 'int'
            s = sprintf('%d', v);
        case 'fix'
            s = sprintf('%.3f', v);
        otherwise
            if v == 0
                s = '0';
            else
                e = floor(log10(abs(v)));
                mant = v / 10^e;
                if round(mant * 100) / 100 >= 10      % e.g. 9.999 -> 10.00
                    mant = mant / 10;
                    e = e + 1;
                end
                s = sprintf('$%.2f\\times10^{%d}$', mant, e);
            end
    end
end
