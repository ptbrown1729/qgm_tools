function [fit_parameters, param_names, fit_fn_handle, std_err] = ...
    fit2D(xvals, yvals, M, variance, functions_cell, ...
    initial_parameters, fixed_parameters, lower_bds, upper_bds,...
    num_params, fit_mode, ignore_nans)
%[FitParameters, PNames, FittedFunctionHandle, StandardErr] = ...
%    fit2D(xvals, yvals, M, Variance, Functions, ...
%    InitialParameters, FixedParameters, LowerLimits, UpperLimits,...
%    num_params, fit_mode)%General fit function. Allows linear combinations of simple defined
%
%functions, initial parameter choice, and fixing parameters.
%%%Arguments%%%
%%%Xvals is a 2d matrix the same size as M specifying the x coordinate of
%each point
%%%Yvals is a 2d matrix the same size as M specifying the y coordinate of
%each point
%%%M is the 2d image to be fit
%%%Variance_i = 1/sigma_i^2 = w_i for point Y_i. If give Variance = [],
%will assume all Variances are 1.
%%%Functions is a cell array of strings specifying which fitting functions
%to use. e.g. Functions = {'gaussian2D','gaussian2D','lorentzian2D'}
%If you specify multiple functions, fit2D uses the sum of all the functions
%you give it. This is one of the reasons fit2D.m is useful. All of the
%functions used must have the form Val = func(P,X,Y,PNamesBool) where P is a
%list of the function parameters to be used during fitting. See
%Fitting/Functions/2D for examples of these functions
%%%InitialParameters is a vector of initial guesses for parameters P of the function.
%%%FixedParameters is a list of ones and zeros. A one fixes that parameter.
%fit2D.m automatically removes these from the fit and uses the initial
%value you sepcified in InitialParameters for that parameter
%%%LowerLimits is a list of lower limit values for P to be used during
%fitting
%%%UpperLimits is a list of upper limit values for P to be used during
%%%fitting
%%%Outputs%%%
%%%FitParameters results of fitting values for P
%%%PNames names of parameters determined from function
%%%FittedFunctionHandle is the best fit function. e.g. you can evaluate
%FittedFunctionHandle(X,Y).
%%%StandardErr standard error values for the parameters P
%TODO: Check if weighted fitting is working

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismatrix(M)
    errorStruct.message = 'M must be a matrix.';
    errorStruct.identifier = 'fit2D:ArgType';
    error(errorStruct)
    %     error('M must be a matrix.');
end

nx = size(M, 2);
ny = size(M, 1);

if ~exist('xvals', 'var') || isempty(xvals)
    [xvals, ~] = meshgrid(1:nx, 1:ny);
end

if ~exist('yvals', 'var') || isempty(yvals)
    [~, yvals] = meshgrid(1:nx, 1:ny);
end

if ~isequal(size(xvals), size(M))
    errorStruct.message = 'Xvals and M were different sizes.';
    errorStruct.identifier = 'MATLAB:samelen';
    error(errorStruct)
end

if ~isequal(size(yvals), size(M))
    errorStruct.message = 'yvals and M were different sizes.';
    errorStruct.identifier = 'MATLAB:samelen';
    error(errorStruct)
end

if ~exist('ignore_nans', 'var')
    ignore_nans = 0;
end

if ignore_nans
    to_fit_indices = ~isnan(M);
    xvals = xvals( to_fit_indices );
    yvals = yvals( to_fit_indices );
    variance = variance( to_fit_indices );
    M = M(to_fit_indices);
end

if any( isnan(xvals(:)) ) || any( isinf(xvals(:)) )
    errorStruct.message = 'Xvals contained nans or infs.';
    errorStruct.identifier = 'fit2D:nanOrinf';
    error(errorStruct);
end

if any( isnan(yvals(:)) ) || any( isinf(yvals(:)) )
    errorStruct.message = 'Yvals contained nans or infs.';
    errorStruct.identifier = 'fit2D:nanOrinf';
    error(errorStruct);
end

if ~exist('variance', 'var') || isempty(variance)
    variance = ones(size(M));
    %disp('Set Variances to ones');
end

if ~isequal(size(variance), size(M))
    errorStruct.message = 'Variance and M were different sizes.';
    errorStruct.identifier = 'MATLAB:samelen';
    error(errorStruct)
end

if ischar(functions_cell) || isa(functions_cell, 'function_handle')
    functions_cell = {functions_cell};
end

if ~iscell(functions_cell)
    errorStruct.message = 'Argument Functions was not a string of characters or a cell array.';
    errorStruct.identifier = 'fit2D:ArgType';
    error(errorStruct);
end

if ~exist('fixed_parameters','var') || isempty(fixed_parameters)
    fixed_parameters = zeros(1,length(initial_parameters));
end

if length(fixed_parameters) ~= length(initial_parameters)
    errorStruct.message = 'FixedParameters was different length than InitialParameters.';
    errorStruct.identifier = 'fit2D:ArgType';
    error(errorStruct);
end

if any(fixed_parameters ~= 0 & fixed_parameters ~= 1)
    fixed_parameters(fixed_parameters ~= 0 & fixed_parameters ~= 1) = 0;
    warning('One or more fixed parameters was not a zero or one. Set offending parameters to zero');
end

if ~exist('upper_bds', 'var') || isempty(upper_bds)
    upper_bds = inf(1, length(initial_parameters));
end

if length(upper_bds) ~= length(initial_parameters)
    errorStruct.message = 'UpperLimits was different length than InitialParameters.';
    errorStruct.identifier = 'fit2D:ArgType';
    error(errorStruct);
end

if ~exist('lower_bds', 'var') || isempty(lower_bds)
    lower_bds = -inf(1, length(initial_parameters));
end

if length(lower_bds) ~= length(initial_parameters)
    errorStruct.message = 'LowerLimits was different length than InitialParameters.';
    errorStruct.identifier = 'fit2D:ArgType';
    error(errorStruct);
end

if (~exist('num_params', 'var') || isempty(num_params) ) && length(functions_cell) == 1
    num_params = length(initial_parameters);
end

if ~exist('fit_mode', 'var') || isempty(fit_mode)
    fit_mode = 'fit';
end

allowed_modes = {'fit', 'lsqnonlin', 'fminunc'};
isallowed = @(x) strcmp(x, fit_mode);
if ~any(cellfun(isallowed, allowed_modes))
    fit_mode = 'fit';
    disp('Warning fit_mode was not an allowed value. Set to "fit".');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin fitting function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NFunctions = length(functions_cell);
fn_handles = {};
ParameterLengths = zeros(1, NFunctions);
param_names = {};

for ii = 1 : NFunctions
    if ischar(functions_cell{ii})
        CurrentHandle = str2func(functions_cell{ii});
        [~, CurrentPNames, ~] = CurrentHandle([], [], []);
        ParameterLengths(ii) = length(CurrentPNames);
        %     FnStrings = cat(2,CurrentFnStr,FnStrings);
    elseif isa(functions_cell{ii}, 'function_handle')
        CurrentHandle = functions_cell{ii};
        ParameterLengths(ii) = num_params(ii);
        CurrentPNames = cell(ParameterLengths(ii), 1);
        [CurrentPNames{:}] = deal('');
    else
        errorStruct.message = 'Function was not a function handle or a string';
        errorStruct.identifier = 'fit2D:value';
        error(errorStruct);
    end
    fn_handles{ii} = CurrentHandle;
    param_names = cat(2,param_names,CurrentPNames);
end

%Create function handle from sum of functions.
%Final handle takes a single P argument. i.e.
%f(P,X,Y)
FullFunctionHandle = @(P,X,Y) 0;
for ii = 1 : NFunctions
    %     CurrentFunctionHandle = str2func(Functions{ii});
    CurrentFunctionHandle = fn_handles{ii};
    PStart = sum(ParameterLengths(1 : ii-1)) + 1;
    PEnd = PStart+ParameterLengths(ii) - 1;
    %CurrentFixed = FixedParameters(PStart:PEnd);
    %FixedValues = InitialParameters(PStart:PEnd).*CurrentFixed;
    FullFunctionHandle = @(P,X,Y) FullFunctionHandle(P, X, Y) + CurrentFunctionHandle(P(PStart : PEnd), X, Y);
    %FullFunctionHandle = @(P,X,Y) FullFunctionHandle(P,X,Y)+CurrentFunctionHandle(P(PStart:PEnd).*(~CurrentFixed)+FixedValues,X,Y,0);
end

FullFunctionHandle = @(Pred,X,Y) FullFunctionHandle(Reduced2FullParams(Pred, initial_parameters, fixed_parameters), X, Y);

%matlab fit function can accept an anonymous function instead of a
%string...however that function must be of the ferm
%f = @(A,B,C,x)...where our functions are of the form f = @([A,B,C],x);
%My ugly solution is to create a string that looks like 'f = @(A,B,C,x) g([A,B,C],x)'
%and then letting matlab evaluate it.
%create string 'ABCDEFG'
ReducedInitialParameters = Full2ReducedParams(initial_parameters, fixed_parameters); %define this here because I need to know how many there are.
Letters = char([1:length(ReducedInitialParameters)] + 'A' - 1);
%use regular expressions to substitute in commas 'A,B,C,D'
ArgString = regexprep(Letters, '[a-zA-Z](?!$)', '$&,');
%build string describing the expression I want to evalute...'f = @(A,B,C,x) g([A,B,C],x)'
StringToEval = strcat('FullFunctionHandle_ListedArguments = @(',ArgString,',x,y) FullFunctionHandle([',ArgString,'],x,y)');
evalc(StringToEval);
%Now convert this new anonymous function into something fit.m can handle.
FullFitFun = fittype(FullFunctionHandle_ListedArguments, 'independent', {'x','y'}, 'dependent', 'z');
RedLowerLims = Full2ReducedParams(lower_bds, fixed_parameters);
RedUpperLims = Full2ReducedParams(upper_bds, fixed_parameters);

if strcmp(fit_mode, 'fit')
    %[RedFitParameters,~,~,~,~,Hessian] = fminunc(FitFunctionHandle,ReducedInitialParameters,options);
    fitobj = fit([xvals(:), yvals(:)], M(:), FullFitFun,...
        'StartPoint', ReducedInitialParameters, 'Weights', variance(:),...
        'Lower', RedLowerLims, 'Upper', RedUpperLims);
    RedFitParameters = coeffvalues(fitobj);
    
    % %Get standard error. Taking advantage
    % of relationship between confidence intervals and standard error.
    NinetyFivePer_ConfInterval = confint(fitobj, 0.95);
elseif strcmp(fit_mode, 'lsqnonlin')
    %use lsqnonlin function
    fit_fn = @(P) (FullFunctionHandle(P, xvals(:), yvals(:)) - M(:)) .* sqrt(variance(:));
    [RedFitParameters, resnorm, residual, exitflag, output, lambda, jacobian] = ...
        lsqnonlin(fit_fn, ReducedInitialParameters, RedLowerLims, RedUpperLims);
    NinetyFivePer_ConfInterval = transpose(nlparci(RedFitParameters,...
                                          residual, 'jacobian', jacobian));

elseif strcmp(fit_mode, 'fminunc')
    warning('Selected fit_mode "fminunc" in fit2D. This does not support upper and lower bounds. Those parameters will be ignored.');
    %use fminunc function
    fit_fn = @(P) sum( (FullFunctionHandle(P, xvals(:), yvals(:)) - M(:)).^2 .* variance(:));
    [RedFitParameters, fval, exitflag , output, grad,hessian] = ...
        fminunc(fit_fn, ReducedInitialParameters);
    NinetyFivePer_ConfInterval = zeros([2, length(RedFitParameters)]);
    warning('Selected fit_mode "fminunc" in fit2D. Standard error calculation not implemented.');
    %StandardErr = transpose(sqrt(diag(Ssqrd*inv(Hessian))))*sqrt(2); %transpose so size is [1,N].

else
    error('fit2D:fit_mode parameter was not an allowed value');    
end

fit_parameters = Reduced2FullParams(RedFitParameters, initial_parameters, fixed_parameters);
fit_fn_handle = @(X, Y) FullFunctionHandle(RedFitParameters, X, Y);

std_err = (NinetyFivePer_ConfInterval(2,:)-NinetyFivePer_ConfInterval(1,:))/3.92;
std_err = Reduced2FullParams(std_err,zeros(size(initial_parameters)),fixed_parameters);
end

%Helper functions to deal with fixed parameters.
function [FullParams] = Reduced2FullParams(RedP, InitParams, FixedParams)
FullParams = zeros(size(InitParams));
FullParams(FixedParams == 1) = InitParams(FixedParams ==1);
FullParams(FixedParams == 0) = RedP;
end

function [RedParams] = Full2ReducedParams(FullParams, FixedParams)
RedParams = FullParams(FixedParams == 0);
end