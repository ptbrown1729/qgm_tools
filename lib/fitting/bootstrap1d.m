function [fitp, err, fitp_distributions] = ...
    bootstrap1d(xvals, yvals, variance, functions,...
    initial_parameters, fixed_parameters, lower_limits, upper_limits,...
    num_params, fit_mode, ntrials, mode)
% Given data with errorbars, take normally sampled points from within the
% errorbar and fit them

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

if strcmp(mode, 'fit')
    [fps, ~, fit_fn_handle, ~, ~] = ...
        fit1D(xvals, yvals, variance, functions,...
        initial_parameters, fixed_parameters, lower_limits, upper_limits,...
        num_params, fit_mode);
    yvals = fit_fn_handle(xvals);
end

% ntrials = 1000;
yerr = 1 ./ sqrt(variance);

warning('set', 'off');
fitp_distributions = zeros(ntrials, length(initial_parameters));
for ii = 1:ntrials
    
    yvals_trial = normrnd(yvals, yerr);

    [fps, ~, ~, stderrs, ~] = ...
        fit1D(xvals, yvals_trial, variance, functions,...
        initial_parameters, fixed_parameters, lower_limits, upper_limits,...
        num_params, fit_mode);
    fitp_distributions(ii, :) = fps;
end
warning('set', 'on');

fitp = mean(fitp_distributions, 1);
err = std(fitp_distributions, [], 1);

end