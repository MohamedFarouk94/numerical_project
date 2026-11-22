%% MAIN  Steepest Descent, Damped Newton and BFGS on the Rosenbrock function.
%
%   Runs the three methods for each starting point listed below, prints the
%   comparison tables (with ready-to-paste LaTeX rows), and saves the
%   figures used in the report.
%
%   Required files: rosenbrock.m, grad.m, hessian.m, backtracking.m,
%   steepest_descent.m, damped_newton.m, bfgs.m, run_experiment.m,
%   plot_results.m (and rosenbrock3d.m for Figure 1).

clc;
close all;

% Set batch_mode = true before running to skip the final pause.
if ~exist('batch_mode', 'var'), batch_mode = false; end

%% Algorithm parameters (shared by all methods and experiments)
params.epsilon  = 1e-6;    % gradient-norm tolerance
params.max_iter = 1000;    % maximum number of iterations
params.alpha0   = 1;       % initial step length
params.beta     = 0.5;     % backtracking factor
params.c        = 1e-4;    % Armijo parameter

n_runs = 1000;             % maximum number of timed runs per method

%% Experiments
% Each experiment has a name, a starting point and (optionally) the names
% of the PNG files for its figures. An empty file struct means "no plots".
files_required = struct('objective', 'figure_2.png', ...
                        'gradient',  'figure_3.png', ...
                        'trajectory', 'figure_4.png');
files_supplementary = struct('objective', 'figure_5.png', ...
                             'gradient',  'figure_6.png', ...
                             'trajectory', 'figure_7.png');

experiments = struct( ...
    'name',  {'Required start',  'Slightly perturbed start', 'Standard start'}, ...
    'x0',    {[-1; 1],           [-1.01; 1.01],               [-1.2; 1]}, ...
    'files', {files_required,    files_supplementary,         []});

%% Run all experiments
for e = 1:numel(experiments)
    ex = experiments(e);

    results = run_experiment(ex.name, ex.x0, params, n_runs);

    if ~isempty(ex.files)
        plot_results(results(1).history, results(2).history, ...
                     results(3).history, ex.x0, params, ex.files);
    end
end

%% Surface plot (Figure 1)
rosenbrock3d('figure_1.png');

if ~batch_mode
    input('Press Enter to close the figures...');
end
close all;
