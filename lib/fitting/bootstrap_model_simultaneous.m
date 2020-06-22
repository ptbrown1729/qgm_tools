function [model_fit_params, model_std_errs, model_fit_params_distribution,...
    curve_fit_params, curve_std_errs, curve_fit_params_distribution, chi_sqr, num_failed_fits] ...
    = bootstrap_model_simultaneous(independent_vars, dependent_vars, errs,...
    model_fn,...
    init_shared_params, fixed_shared_params, shared_lbs, shared_ubs,...
    init_curve_params, fixed_curve_params, curve_lbs, curve_ubs, num_curve_params,...
    ntrials, mode)
% do boostrap fitting for fit_mode_simultaneous
% see the documentation for fit_mode_simultaneous
% 
% mode. If this parameter is 'data', the bootstrap is performed by, for
% each trial, selecting random values in the vicinity of each data point.
% These values are normally distributed about the datapoint with standard
% deviation given by the errorbar. This process if repeated for ntrials,
% and the fit parameters are recorded for each trial. Final fit parameters
% and their errorbars are taken as the mean and standard deviation of the
% resulting distributions. If mode is 'fit', then the model is first fit to
% the data, and then the model function is sampled at each datapoint.
% Then normally distributed values are chosen around these samples for each
% trial.

if ~iscell(independent_vars)
    independent_vars = {independent_vars};
end

if ~iscell(dependent_vars)
    dependent_vars = {dependent_vars};
end

if ~iscell(errs)
    errs = {errs};
end

if ~exist('ntrials', 'var') || isempty(ntrials)
    ntrials = 300;
end

if ~exist('mode', 'var') || isempty(mode)
    mode = 'data';
end

allowed_modes = {'data', 'fit'};
cmp_fun = @(c) strcmp(c, mode);
if ~any(cellfun(cmp_fun, allowed_modes))
    error('mode variable in bootstrap_model_simultaneous.m was not an allowed mode.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_sets = length(independent_vars);
num_shared_params = length(init_shared_params);

if strcmp(mode, 'fit')
    [model_fit_params, ~, curve_fit_params, ~, ~, exitflag] = ...
        fit_model_simultaneous(independent_vars, dependent_vars, errs,...
    model_fn,...
    init_shared_params, fixed_shared_params, shared_lbs, shared_ubs,...
    init_curve_params, fixed_curve_params, curve_lbs, curve_ubs, num_curve_params);

    for jj = 1:num_sets
        curve_fit_params_current = curve_fit_params((jj-1) * num_curve_params + 1 : jj * num_curve_params);
        dependent_vars{jj} = model_fn(model_fit_params, curve_fit_params_current, independent_vars{jj}); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('set', 'off');
model_fit_params_distribution = zeros(ntrials, num_shared_params);
curve_fit_params_distribution = zeros(ntrials, num_curve_params * num_sets);

max_trials = 2*ntrials;
n_successes = 0;
n_runs = 0;
% TODO: need to get rid of fits that fail to converge!
while (n_successes < ntrials) && (n_runs < max_trials)

    
    dependent_vars_trial = cell(1, num_sets);
    for jj = 1:num_sets
       dependent_vars_trial{jj} = normrnd(dependent_vars{jj}, errs{jj});
    end

    [model_fit_params, ~, curve_fit_params, ~, ~, exitflag] = ...
        fit_model_simultaneous(independent_vars, dependent_vars_trial, errs,...
    model_fn,...
    init_shared_params, fixed_shared_params, shared_lbs, shared_ubs,...
    init_curve_params, fixed_curve_params, curve_lbs, curve_ubs, num_curve_params);

    if exitflag ~= 0
        n_successes = n_successes + 1;
        model_fit_params_distribution(n_successes, :) = model_fit_params;
        curve_fit_params_distribution(n_successes, :) = curve_fit_params;
    end
  
    n_runs = n_runs + 1;
end
warning('set', 'on');

num_failed_fits = n_runs - n_successes;
if (n_successes < ntrials)
    model_fit_params_distribution = model_fit_params_distribution(1:n_successes, :);
    curve_fit_params_distribution = curve_fit_params_distribution(1:n_successes, :);
    warning('Maximum number of trials exceeded before ntrials fits succeeded in bootstrap_model_simultaneous.m. Returning successful fit rsults');
end

model_fit_params = mean(model_fit_params_distribution, 1);
model_std_errs = std(model_fit_params_distribution, [], 1);

curve_fit_params = mean(curve_fit_params_distribution, 1);
curve_std_errs = std(curve_fit_params_distribution, [], 1);

% TODO: set this!
chi_sqr = 0;

end