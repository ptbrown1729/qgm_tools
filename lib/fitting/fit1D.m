function [FitParameters, PNames, FittedFunctionHandle, StandardErr, reduced_chi_sqr] = ...
    fit1D(xvals, yvals, variance, Functions,...
    InitialParameters, FixedParameters, LowerLimits, UpperLimits,...
    num_params, fit_mode)
%function [FitParameters, PNames, FittedFunctionHandle, StandardErr, RedChiSqr] = ...
%     fit1D(xvals, yvals, variance, Functions,...
%     InitialParameters, FixedParameters, LowerLimits, UpperLimits,...
%     num_params, fit_mode)
%
%General fit function. Allows linear combinations of simple define 
%functions, initial parameter choice, and fixing parameters.
%
%Functions is a (possibly mixed) list of function names and function handles. 
%Each function should have the 
%form f(P,X,Plot_Bool), where P is a vector of parameters
%InitialParameters is a vector with the same order as P.
%FixedParameters is a vector of ones and zeros. A one fixes that parameter.
%If give Variance = [], will assume all Variances are 1.
%Variance_i = 1/sigma_i^2 = w_i for point Y_i.
%LowerLimits is a list of LowerLimits to be used
%UpperLimits is a list of UpperLimits to be used
%num_params is a list of the same size as Functions indicating how many
%parameters each function uses. This is necessary in the case of function
%handles, because there is no way to tell how many parameters these take
%otherwise.
%fit_mode = "fit", "lsqnonlin", or "fminunc". Default is "fit".
%%%%
%FitParameters, a vector of best fit parameters with same order as P.
%PNames describes the fit parameters
%FittedFunctionHandle is the best fit function. e.g. you can evaluate 
%FittedFunctionHandle(X).
%StandardErr is a vector of standard errors with the same order as P.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isequal(size(xvals), size(yvals))
    if length(xvals) == length(yvals)
        yvals = transpose(yvals);
    else
        errorStruct.message = 'Xvals and Yvals were different sizes.';
        errorStruct.identifier = 'MATLAB:samelen';
        error(errorStruct)
    end
end

if sum(isnan(xvals)) ~= 0 || sum(isinf(xvals)) ~= 0
    errorStruct.message = 'Xvals contained nans or infs.';
    errorStruct.identifier = 'fit1D:nanOrinf';
    error(errorStruct);
end

if sum(isnan(yvals)) ~= 0 || sum(isinf(yvals))~=0
    errorStruct.message = 'Yvals contained nans or infs.';
    errorStruct.identifier = 'fit1D:nanOrinf';
    error(errorStruct);
end

if ~exist('variance', 'var') || isempty(variance)
    variance = ones(size(xvals));
    %disp('Set Variances to ones');
end

if length(xvals) ~= length(variance)
    variance = ones(size(xvals));
    %disp('Set Variances to ones');
end

if ischar(Functions) || isa(Functions, 'function_handle')
    Functions = {Functions};
end

if ~iscell(Functions)
    errorStruct.message = 'Argument Functions was not a string of characters, fn handle, or a cell array.';
    errorStruct.identifier = 'fit1D:ArgType';
    error(errorStruct);
end

if ~exist('InitialParameters', 'var') || isempty(InitialParameters)
    error('Initial parameters was empty');
end

if ~exist('FixedParameters', 'var') || isempty(FixedParameters)
    FixedParameters = zeros(1,length(InitialParameters));
end

if length(FixedParameters) ~= length(InitialParameters)
    FixedParameters = zeros(1, length(InitialParameters));
    warning('FixedParameters was a different size than InitialParameters. Set to zeros.');
end

if any(FixedParameters ~= 0 & FixedParameters ~= 1)
    FixedParameters(FixedParameters ~= 0 & FixedParameters ~= 1) = 0;
    warning('One or more fixed parameters was not a zero or one. Set offending parameters to zero');
end

if ~exist('UpperLimits', 'var') || isempty(UpperLimits)
    UpperLimits = inf(length(InitialParameters),1);
end

if length(UpperLimits) ~= length(InitialParameters)
    UpperLimits = inf(length(InitialParameters),1);
end

if ~exist('LowerLimits', 'var') || isempty(LowerLimits)
    LowerLimits = -inf(length(InitialParameters),1);
end

if length(LowerLimits) ~= length(InitialParameters)
    LowerLimits = -inf(length(InitialParameters), 1);
end

if ~exist('num_params', 'var') && length(Functions) == 1
    num_params = length(InitialParameters);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NFunctions = length(Functions);
fn_handles = {};
ParameterLengths = zeros(1,NFunctions);
PNames = {};
% FnStrings = {};

%Generate and store information about functions.
for ii=1:NFunctions
    if ischar(Functions{ii})
        CurrentHandle = str2func(Functions{ii});
        [~,CurrentPNames,~] = CurrentHandle([],[]);
        ParameterLengths(ii) = length(CurrentPNames);
        %     FnStrings = cat(2,CurrentFnStr,FnStrings);
    elseif isa(Functions{ii}, 'function_handle')
        CurrentHandle = Functions{ii};
        ParameterLengths(ii) = num_params(ii);
        CurrentPNames = cell(ParameterLengths(ii), 1);
        [CurrentPNames{:}] = deal('');
    else
        errorStruct.message = 'Function was not a function handle or a string';
        errorStruct.identifier = 'fit1D:value';
        error(errorStruct);
    end
    fn_handles{ii} = CurrentHandle;
    PNames = cat(2,PNames,CurrentPNames);
end

%Create anonymous function which is the sum of Functions.
FullFunctionHandle = @(P, x) 0;
for ii=1:NFunctions
%     CurrentFunctionHandle = str2func(Functions{ii});
    CurrentFunctionHandle = fn_handles{ii};
    PStart = sum(ParameterLengths(1:ii-1)) + 1;
    PEnd = PStart+ParameterLengths(ii) - 1;
    FullFunctionHandle = @(P, x) FullFunctionHandle(P, x) + CurrentFunctionHandle(P(PStart:PEnd), x);
    %FullFunctionHandle = @(P,x) FullFunctionHandle(P,x)+CurrentFunctionHandle(P(PStart:PEnd).*(~CurrentFixed)+FixedValues,x); 
end
%Now create a new anonymous function that doesn't include any of our fixed
%parameters. This is all we need for lsqnonlin
FullFunctionHandle = @(Pred,x) FullFunctionHandle(Reduced2FullParams(Pred,InitialParameters,FixedParameters),x);

%for the matlab "fit" function we have to do a bit more work
%that fn can accept an anonymous function instead of a
%string...however that function must be of the form
%f = @(A,B,C,x)...where our functions are of the form f = @([A,B,C],x);
%My ugly solution is to create a string that looks like 'f = @(A,B,C,x) g([A,B,C],x)'
%and then letting matlab evaluate it.
%create string 'ABCDEFG'
ReducedInitialParameters = Full2ReducedParams(InitialParameters,FixedParameters); %define this here because I need to know how many there are.
Letters = char([1:length(ReducedInitialParameters)]+'A'-1);
%use regular expressions to substitute in commas 'A,B,C,D'
ArgString = regexprep(Letters,'[a-zA-Z](?!$)','$&,');
%build string describing the expression I want to evalute...'f = @(A,B,C,x) g([A,B,C],x)'
StringToEval = strcat('FullFunctionHandle_ListedArguments = @(', ...
    ArgString, ',x) FullFunctionHandle([', ArgString, '],x)');
evalc(StringToEval);
%Now convert this new anonymous function into something fit.m can handle.

%Do fitting and extract parameters to return.
try
    
    RedLowerLims = Full2ReducedParams(LowerLimits, FixedParameters);
    RedUpperLims = Full2ReducedParams(UpperLimits, FixedParameters);
    if strcmp(fit_mode, 'fit')
        %use fit function
        FullFitFun = fittype(FullFunctionHandle_ListedArguments, 'independent', 'x', 'dependent', 'y');
        fitobj = fit(xvals(:), yvals(:), FullFitFun, ...
            'StartPoint', ReducedInitialParameters,...
            'Weights', variance,...
            'Lower', RedLowerLims, 'Upper', RedUpperLims);
        
        RedFitParameters = coeffvalues(fitobj);
        NinetyFivePer_ConfInterval = transpose(confint(fitobj,0.95));
        reduced_chi_sqr = 0;
        
    elseif strcmp(fit_mode, 'lsqnonlin')
        %use lsqnonlin function
        fit_fn = @(P) (FullFunctionHandle(P, xvals(:)) - yvals(:)) .* sqrt(variance(:));
        [RedFitParameters, resnorm, residual, exitflag, output, lambda, jacobian] = ...
            lsqnonlin(fit_fn, ReducedInitialParameters, RedLowerLims, RedUpperLims);
       reduced_chi_sqr = resnorm / (length(xvals) - length(RedFitParameters));

               % SHOULDNT THIS BE REDFITPARAMS???
%         NinetyFivePer_ConfInterval = nlparci(ReducedInitialParameters, residual, 'jacobian', jacobian);  
        NinetyFivePer_ConfInterval = nlparci(RedFitParameters, residual, 'jacobian', jacobian);
        
    elseif strcmp(fit_mode, 'fminunc')
        warning('Selected fit_mode "fminunc" in fit1D. This does not support upper and lower bounds. Those parameters will be ignored.');
        %use fminunc function
        fit_fn = @(P) sum( (FullFunctionHandle(P, xvals(:)) - yvals(:)).^2 .* variance(:));
        [RedFitParameters, fval, exitflag , output, grad,hessian] = ... 
            fminunc(fit_fn, ReducedInitialParameters);
        NinetyFivePer_ConfInterval = zeros([size(RedFitParameters), 2]);
        warning('Selected fit_mode "fminunc" in fit1D. Standard error calculation not implemented.');
        reduced_chi_sqr = 0;
        %StandardErr = transpose(sqrt(diag(Ssqrd*inv(Hessian))))*sqrt(2); %transpose so size is [1,N].
    else
        error('fit1D:fit_mode parameter was not an allowed value');
    end
    FitParameters = Reduced2FullParams(RedFitParameters, InitialParameters, FixedParameters);
    %TODO right now fit1D and fit2D return different shaped stderrs. Should
    %standardize these...
    StandardErr = (NinetyFivePer_ConfInterval(:,2)-NinetyFivePer_ConfInterval(:,1))/3.92;
    %Divide by two to get `one sided` confidence interval. i.e. +/-
    %Divide by 1.96 to get standard deviation, since integrating normal
    %distribution out to 1.96*sigma is what gives 95% confidence interval.
    StandardErr = Reduced2FullParams(StandardErr, InitialParameters, FixedParameters);
    FittedFunctionHandle = @(x) FullFunctionHandle(RedFitParameters,x);
catch Exception
    FitParameters = InitialParameters;
    StandardErr = zeros(1,length(InitialParameters));
    FittedFunctionHandle = @(x) 0;
    reduced_chi_sqr = 0;
    MsgString = getReport(Exception);
    disp(MsgString);
end
 
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