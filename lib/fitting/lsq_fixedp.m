function [fit_parameters, std_errs, reduced_chi_sqr, exit_flag] = ...
    lsq_fixedp(function_handle, init_parameters, fixed_parameters, lower_bounds, upper_bounds)
% [fit_parameters, std_errs, reduced_chi_sqr, exit_flag] = ...
%     lsq_fixedp(function_handle, init_parameters, fixed_parameters, lower_bounds, upper_bounds)
% For easy use of lsqnonlin with some parameters fixed. Like lsqnonlin,
% this accepts a function, where the square fo this function is going to be
% minimized. So, e.g., if you want to fit points {y_i} to a model function
% f(P, x) for points {x_i}, you would need to supply function g(P) =
% [f(P,x_1) - y_1, f(P, x_2) - y_2, ...]
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~ischar(function_handle) && ~isa(function_handle, 'function_handle')
    errorStruct.message = 'Argument Functions was not a string of characters function handle.';
    errorStruct.identifier = 'lsq_fixedp:ArgType';
    error(errorStruct);
end

if ~exist('fixed_parameters', 'var') || isempty(fixed_parameters)
    fixed_parameters = zeros(1,length(init_parameters));
end

if length(fixed_parameters) ~= length(init_parameters)
    fixed_parameters = zeros(1, length(init_parameters));
    warning('FixedParameters was a different size than InitialParameters. Set to zeros.');
end

if any(fixed_parameters ~= 0 & fixed_parameters ~= 1)
    fixed_parameters(fixed_parameters ~= 0 & fixed_parameters ~= 1) = 0;
    warning('One or more fixed parameters was not a zero or one. Set offending parameters to zero');
end

if ~exist('upper_bounds', 'var') || isempty(upper_bounds)
    upper_bounds = inf(length(init_parameters),1);
end

if length(upper_bounds) ~= length(init_parameters)
    upper_bounds = inf(length(init_parameters),1);
end

if ~exist('lower_bounds', 'var') || isempty(lower_bounds)
    lower_bounds = -inf(length(init_parameters),1);
end

if length(lower_bounds) ~= length(init_parameters)
    lower_bounds = -inf(length(init_parameters), 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(function_handle)
    fn_handle = str2func(function_handle);
elseif isa(function_handle, 'function_handle')
    fn_handle = function_handle;
else
    errorStruct.message = 'Function was not a function handle or a string';
    errorStruct.identifier = 'lsq_fixedp:value';
    error(errorStruct);
end

%Now create a new anonymous function that doesn't include any of our fixed
%parameters. This is all we need for lsqnonlin
FullFunctionHandle = @(Pred) fn_handle(Reduced2FullParams(Pred,init_parameters,fixed_parameters));

%Do fitting and extract parameters to return.
ReducedInitialParameters = Full2ReducedParams(init_parameters,fixed_parameters); %define this here because I need to know how many there are.
redlowerbs = Full2ReducedParams(lower_bounds,fixed_parameters);
redupperbs = Full2ReducedParams(upper_bounds,fixed_parameters);
    %use lsqnonlin function
    [RedFitParameters, resnorm, residual, exit_flag, output, lambda, jacobian] = ...
        lsqnonlin(FullFunctionHandle, ReducedInitialParameters, redlowerbs, redupperbs);
    reduced_chi_sqr = resnorm / (length(FullFunctionHandle(RedFitParameters)) - length(RedFitParameters));
    %get standard error from confidence intervals
%     NinetyFivePer_ConfInterval = nlparci(ReducedInitialParameters, residual, 'jacobian', jacobian);
    NinetyFivePer_ConfInterval = nlparci(RedFitParameters, residual, 'jacobian', jacobian);
    std_errs = (NinetyFivePer_ConfInterval(:,2) - NinetyFivePer_ConfInterval(:,1)) / 3.92;
    %convert back to original params
    fit_parameters = Reduced2FullParams(RedFitParameters, init_parameters, fixed_parameters);
    std_errs = Reduced2FullParams(std_errs, init_parameters, fixed_parameters);

end

%Helper functions to deal with fixed parameters
function [FullParams] = Reduced2FullParams(RedP,InitParams,FixedParams)
FullParams = zeros(size(InitParams));
FullParams(FixedParams == 1) = InitParams(FixedParams == 1);
FullParams(FixedParams == 0) = RedP;
end

function [RedParams] = Full2ReducedParams(FullParams,FixedParams)
RedParams = FullParams(FixedParams == 0);
end