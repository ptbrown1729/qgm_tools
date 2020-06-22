function [model_fit_params, model_std_errs, curve_fit_params, curve_std_errs, chi_sqr, exitflag] ...
    = fit_model_simultaneous(independent_vars, dependent_vars, errs,...
    model_fn,...
    init_shared_params, fixed_shared_params, shared_lbs, shared_ubs,...
    init_curve_params, fixed_curve_params, curve_lbs, curve_ubs, num_curve_params)
% This function is useful for fitting a family of curves to a shared model.
% We assume that the model has two kinds of parameters, (1) those that all
% of the curves share (e.g. temperature) and (2) those that vary from curve
% to curve (e.g. amplitude, offset, etc.).
%
% independent_vars, a cell array of 1D vectors giving independent variable
% values for each curve. (e.g. time points for different traces). Different
% traces are not required to have the same number of time points.
%
% dependent_vars, a cell array of 1D vectors giving dependent values for
% each curve. (e.g. amplitude at each time point). Different curves are not
% required to have the same number of dependent variable points. However,
% dependent_vars{ii} and independent_vars{ii} must be the same size.
%
% errs, a cell array of 1D vectors giving the uncertainty in the
% dependenent variable for each curve. variances{ii} must be the same size 
% as indepenent_vars{ii}
%
% model_fn(shared_parameters, curve_parameters, independent_var) defines
% the model to fit. It takes three arguments. The first argument, shared
% parameters, accepts an array of parameters which will be common to all
% traces. The second argument, curve_parameters, accepts an array of
% parameters which will differ between curves. The final argument,
% independent_var, is the independent variable. An example of such a
% function might be sine([Phase], [Amplitude, Offset], Time), if you are
% fitting a family of sine curves which are all required to have the same
% phase, but might differ in amplitude.
%
% init_shared_params, 
%
% fixed_shared_params, 
%
% shared_lbs, 
%
% shared_ubs,
%
% init_curve_params,
%
% fixed_curve_params,
%
% curve_lbs, 
% 
% curve_ubs, 
% 
% num_curve_params

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Argument checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(independent_vars)
    independent_vars = {independent_vars};
end

num_sets = length(independent_vars);
num_shared_params = length(init_shared_params);

if ~iscell(dependent_vars)
    dependent_vars = {dependent_vars};
end

if ~iscell(errs)
    errs = {errs};
end


if ~exist('fixed_shared_params', 'var') || isempty(fixed_shared_params)
    fixed_shared_params = zeros(size(init_shared_params));
end
fixed_shared_params = reshape(fixed_shared_params, [1, length(fixed_shared_params)]);

if ~exist('fixed_curve_params', 'var') || isempty(fixed_curve_params)
    fixed_curve_params = zeros(size(init_curve_params));
end
fixed_curve_params = reshape(fixed_curve_params, [1, length(fixed_curve_params)]);

if ~exist('shared_lbs', 'var') || isempty(shared_lbs)
    shared_lbs = -inf * ones(size(init_shared_params));
end
shared_lbs = reshape(shared_lbs, [1, length(shared_lbs)]);

if ~exist('shared_ubs', 'var') || isempty(shared_ubs)
    shared_ubs = inf * ones(size(init_shared_params));
end
shared_ubs = reshape(shared_ubs, [1, length(shared_ubs)]);

% if init curve parameters 
init_curve_params = reshape(init_curve_params, [1, length(init_curve_params)]);
if length(init_curve_params) == num_sets * num_curve_params
elseif length(init_curve_params) == num_curve_params
    init_curve_params = repmat(init_curve_params, [1, num_sets]);
else
    error('init_curve_params was incorrect length in fit_model_simultaneous');
end

if ~exist('fixed_curve_params', 'var') || isempty(fixed_curve_params)
    fixed_curve_params = zeros(size(init_curve_params));
end

fixed_curve_params = reshape(fixed_curve_params, [1, length(fixed_curve_params)]);
if length(fixed_curve_params) == num_sets * num_curve_params
elseif length(fixed_curve_params) == num_curve_params
    fixed_curve_params = repmat(fixed_curve_params, [1, num_sets]);
else
    error('fixed_curve_params was incorrect length in fit_model_simultaneous');
end


if ~exist('curve_lbs', 'var') || isempty(curve_lbs)
    curve_lbs = -inf * ones(size(init_curve_params));
end
curve_lbs = reshape(curve_lbs, [1, length(curve_lbs)]);

if ~exist('curve_ubs', 'var') || isempty(curve_ubs)
    curve_ubs = inf * ones(size(init_curve_params));
end
curve_ubs = reshape(curve_ubs, [1, length(curve_ubs)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create model function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FitFnHandles = {};
for ii = 1:num_sets
    CurrFitFn = @(P, M) ( model_fn(P, M, independent_vars{ii}) - dependent_vars{ii} ) ./ errs{ii};
    FitFnHandles = horzcat(FitFnHandles, {CurrFitFn});
end

% define function handle for simultaneous fitting. This
% requires concatenating our cell array of function handles
% into a vector.
fit_fn = @(P, M) [];
%             num_curve_params = 2;
for ii = 1:num_sets
    fit_fn = @(P, M) [fit_fn(P, M), ...
        FitFnHandles{ii}(P,...
        M( num_curve_params * (ii-1) + 1 : num_curve_params * ii))];
end
fit_fn = @(Var) fit_fn(Var(1:num_shared_params), Var(num_shared_params + 1:end));
%[[ModelP1, ModelP2, ModelP3, ModelP4,...], [Curve11, Curve12,..., Curve1n], [Curve21, Curve22, ...], ... , [Curvem1,...]]

% initial parameters and limits
init_params = [init_shared_params, init_curve_params];
fixed_params = [fixed_shared_params, fixed_curve_params];
lbs = [shared_lbs, curve_lbs];
ubs = [shared_ubs, curve_ubs];

%fit function
[fit_params, std_errs, chi_sqr, exitflag] = ...
    lsq_fixedp(fit_fn, init_params, fixed_params, lbs, ubs);

% parameters that will be output arguments
model_fit_params = fit_params(1:num_shared_params);
curve_fit_params = fit_params(num_shared_params + 1 : end);

model_std_errs = std_errs(1:num_shared_params);
curve_std_errs = std_errs(num_shared_params + 1 : end);

end