function plot_results(h_sd, h_newton, h_bfgs, x0, params, files)
%PLOT_RESULTS  Convergence plots and trajectories for the three methods.
%
%   PLOT_RESULTS(h_sd, h_newton, h_bfgs, x0, params, files) creates
%       - the objective value versus iteration (semilogarithmic),
%       - the gradient norm versus iteration, with the tolerance line,
%       - the trajectories over the contour map of the Rosenbrock function.
%
%   files is an optional struct with the fields objective, gradient and
%   trajectory. Each non-empty field is the name of a PNG file in which
%   the corresponding figure is saved (300 dpi).
%
%   If any history is long (more than 100 points), the iteration axis of
%   the convergence plots is logarithmic and shows k + 1.

    if nargin < 6 || isempty(files), files = struct(); end

    setup_graphics();

    output_dpi  = 300;
    value_floor = 1e-16;      % zeros cannot be drawn on a logarithmic axis

    %% Styles (SD is drawn last, thickest and with a hollow marker, so
    %% that BFGS stays visible when the two curves coincide)
    st_sd     = make_style([0 0 1],   ':', 'o', 2,   10);
    st_newton = make_style([0 0.6 0], '-', 's', 1.5, 7);
    st_bfgs   = make_style([1 0 0],   '-', '^', 1.5, 5);

    max_len  = max([numel(h_sd.f), numel(h_newton.f), numel(h_bfgs.f)]);
    long_run = max_len > 100;
    if long_run
        % Markers and dots are unreadable with hundreds of points.
        st_sd.marker = '';  st_newton.marker = '';  st_bfgs.marker = '';
        st_sd.line = '-';
    end

    ttl = sprintf('x_0 = (%g, %g)', x0(1), x0(2));

    %% Figure: objective function convergence
    fig = figure('Name', 'Objective function convergence');
    hold on;
    p_bfgs   = curve(h_bfgs.f,   value_floor, long_run, st_bfgs);
    p_newton = curve(h_newton.f, value_floor, long_run, st_newton);
    p_sd     = curve(h_sd.f,     value_floor, long_run, st_sd);
    finish_axes(long_run);
    ylabel('Objective Function Value');
    title(['Objective Function Convergence, ' ttl]);
    if long_run
        legend_location = 'southeast';
    else
        legend_location = 'northeast';
    end
    legend([p_sd, p_newton, p_bfgs], ...
        'Steepest Descent', 'Damped Newton', 'BFGS', ...
        'Location', legend_location);
    save_figure(fig, files, 'objective', output_dpi);

    %% Figure: gradient norm convergence
    fig = figure('Name', 'Gradient norm convergence');
    hold on;
    p_bfgs   = curve(h_bfgs.grad_norm,   value_floor, long_run, st_bfgs);
    p_newton = curve(h_newton.grad_norm, value_floor, long_run, st_newton);
    p_sd     = curve(h_sd.grad_norm,     value_floor, long_run, st_sd);

    % Stopping tolerance
    k_first = 0 + long_run;
    k_last  = max_len - 1 + long_run;
    p_tol = plot([k_first, k_last], params.epsilon * [1 1], '--', ...
        'Color', [0.3 0.3 0.3], 'LineWidth', 1);

    finish_axes(long_run);
    ylabel('Gradient Norm');
    title(['Gradient Norm Convergence, ' ttl]);
    % Choose an empty corner so the legend does not hide any curve:
    % bottom centre for short runs, bottom right for long runs.
    if long_run
        legend_location = 'southeast';
    else
        legend_location = 'south';
    end
    legend([p_sd, p_newton, p_bfgs, p_tol], ...
        'Steepest Descent', 'Damped Newton', 'BFGS', ...
        sprintf('Tolerance (%g)', params.epsilon), ...
        'Location', legend_location);
    save_figure(fig, files, 'gradient', output_dpi);

    %% Figure: trajectories over the contour map
    fig = figure('Name', 'Optimization trajectories');

    [xg, yg] = meshgrid(-2:0.02:2, -1:0.02:3);
    z = (1 - xg).^2 + 100 * (yg - xg.^2).^2;
    levels = logspace(-1, 3.5, 30);       % log-spaced level curves

    contour(xg, yg, z, levels);
    hold on;
    p_bfgs   = trajectory(h_bfgs,   st_bfgs);
    p_newton = trajectory(h_newton, st_newton);
    p_sd     = trajectory(h_sd,     st_sd);

    p_initial = plot(x0(1), x0(2), 'kp', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'y');
    p_minimum = plot(1, 1, 'r*', 'MarkerSize', 12, 'LineWidth', 2);

    axis([-2 2 -1 3]);
    xlabel('x');
    ylabel('y');
    title(['Optimization Trajectories on Rosenbrock Contours, ' ttl]);
    legend([p_sd, p_newton, p_bfgs, p_initial, p_minimum], ...
        'Steepest Descent', 'Damped Newton', 'BFGS', ...
        'Initial Point', 'Global Minimum', ...
        'Location', 'northwest');
    grid on;
    hold off;
    save_figure(fig, files, 'trajectory', output_dpi);
end


%% ------------------------------------------------------------------------
function setup_graphics()
% Use the gnuplot toolkit in Octave (when installed); MATLAB needs nothing.
    if exist('OCTAVE_VERSION', 'builtin') ~= 0
        if any(strcmp(available_graphics_toolkits(), 'gnuplot'))
            graphics_toolkit('gnuplot');
        end
    end
end


function st = make_style(color, line, marker, width, msize)
    st = struct('color', color, 'line', line, 'marker', marker, ...
                'width', width, 'size', msize);
end


function h = curve(y, value_floor, long_run, st)
% Draw one convergence curve on the current (logarithmic) axes.
    k = 0:numel(y) - 1;
    if long_run, k = k + 1; end
    h = plot(k, max(y, value_floor), [st.line st.marker], ...
        'Color', st.color, 'LineWidth', st.width, 'MarkerSize', st.size);
end


function finish_axes(long_run)
% Logarithmic y axis (and x axis for long runs), grid, labels.
    set(gca, 'YScale', 'log');
    if long_run
        set(gca, 'XScale', 'log');
        xlabel('Iteration + 1 (log scale)');
    else
        xlabel('Iteration');
    end
    grid on;
    hold off;
end


function h = trajectory(hist, st)
% Draw one optimization path on the current axes.
    h = plot(hist.x(1, :), hist.x(2, :), [st.line st.marker], ...
        'Color', st.color, 'LineWidth', st.width, 'MarkerSize', st.size);
end


function save_figure(fig, files, field, dpi)
% Save the figure if a file name was given for it.
    drawnow;
    if isfield(files, field) && ~isempty(files.(field))
        print(fig, files.(field), '-dpng', sprintf('-r%d', dpi));
    end
end
