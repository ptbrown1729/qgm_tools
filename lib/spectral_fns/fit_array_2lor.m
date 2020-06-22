function [fit_params, std_errs, fit_ref] = fit_array_2lor(array, array_unc, omegas, xs, ys,...
                                                 init_params, fixed_params, lbs, ubs,...
                                                 num_lors, track_fitting, dist_fn, fit_order)
% Given a 3D array, for each (y, x) position specified by the first two
% dimensions, fit a lorentzian to the data along the third dimension.
%
% array: 3D array where the first dimension represents y-position, the
% second x-position, and the third some other variable (such as frequency).
% If array is a 2D array a singleton dimension will be appended to the
% first dimension.
%
% array_unc: Same size as array, representing uncertainty associated with
% each point.
%
% omegas: Either a 3D array of the same size as array or a vector with
% length equal to size(array, 3). z-axis positions for each pixel.
%
% xs: x positions, vector of length size(array, 2)
%
% ys: y positions, vector of length size(array, 1)
%
% init_params:
%
% fixed_params:
%
% lbs:
%
% ubs:
%
% num_lors: number of lorentzians to use when fitting each pixel
%
% track_fitting: if 1, will adjust initial parameters based on previous
% fit parameters.
%
% dist_fn: optionally define a function which determines the distance
% between array pixels. The initial parameters for a lorentzian fit of one
% pixel are determined by the fit parameters at a previously fit nearby
% pixel
%
% fit_order: an array the same size as the first two dimensions of array.
% The values define the order in which the pixels will be fit. Smaller
% values will be fit before larger values.

array_2d = ndims(array) == 2;

if array_2d
    % if array is 2D instead of 3D, create 3D array by adding singleton
    % y-dimension
    array = permute(array, [3, 1, 2]);
end

if exist('array_unc', 'var') && ~isempty(array_unc) && array_2d
    array_unc = permute(array_unc, [3, 1, 2]);
end

if ~exist('array_unc', 'var') || isempty(array_unc)
    array_unc = ones( size(array) );
end

if exist('omegas', 'var') && ~isempty(omegas) && array_2d
    omegas = permute(omegas, [3, 1, 2]);
end

if ~exist('omegas', 'var') || isempty(omegas)
    omegas = 1 : size(array, 3);
end

if ~exist('ys', 'var') || isempty(ys)
    ys = 1 : size(array, 1);
end

if ~exist('xs', 'var') || isempty(xs)
    xs = 1 : size(array, 2);
end

if ~exist('dist_fn', 'var') || isempty(dist_fn)
    dist_fn = @(x1, y1, x2, y2) sqrt( (x1 - x2).^2 + (y1 - y2).^2);
end

if ~exist('num_lors', 'var') || isempty(num_lors)
    num_lors = 1;
end

% for all lorentzian

if ~exist('init_params', 'var') || isempty(init_params)
    % for single lorentzian
    if ndims(omegas) == 3
        domega = omegas(:, :, 2:end) - omegas(:, :, 1:end-1);
    else
        domega = omegas(2:end) - omegas(1:end - 1);
    end
    init_params_1 = [ mean(omegas(:)), 1, max(array(:)), 0];
    fixed_params_1 = [0, 0, 0, 0];
    fixed_params_2 = [0, 0, 0, 1];
    lbs_1 = [ min(omegas(:)), min(domega(:)), -inf, -inf];
    ubs_1 = [ max(omegas(:)), inf, inf, inf];
    % all lorentzians
    init_params = repmat(init_params_1, [1, num_lors]);
    init_params(1:4:end) = linspace( min(omegas), max(omegas), num_lors);
end

if ~exist('fixed_params', 'var') || isempty(fixed_params)
    fixed_params = horzcat( fixed_params_1, repmat(fixed_params_2, [1, num_lors - 1]) );
end

if ~exist('ubs', 'var') || isempty(ubs)
    ubs = repmat(ubs_1, [1, num_lors]);
end

if ~exist('lbs', 'var') || isempty(lbs)
    lbs = repmat(lbs_1, [1, num_lors]);
end

if numel(init_params) ~= 4 * num_lors || ...
   numel(fixed_params) ~= 4 * num_lors || ...
   numel(ubs) ~= 4 * num_lors || ...
   numel(lbs) ~= 4 * num_lors

    error('number of lorentzians selected and initial parameters are not consistent');
end


fn_name = repmat({'lorentzian1D'}, [1, num_lors]);

% output variables
already_fit = zeros( size(array(:, :, 1)) );
fit_params = zeros( length(ys), length(xs), length(init_params) );
std_errs = zeros( length(ys), length(xs), length(init_params) );
fit_ref = zeros( length(ys), length(xs) );

[xx, yy] = meshgrid(xs, ys);

if ~exist('fit_order', 'var') || isempty(fit_order)
    % fit along each column
    [xinds, yinds] = meshgrid( 1 : size(array, 2), 1 : size(array, 1) );
else
    [~, order] = sort(fit_order(:));
    [yinds, xinds] = ind2sub( size(fit_order), order);
end

for kk = 1 : numel(yinds)
    jj = xinds(kk);
    ii = yinds(kk);
   
    if isequal( size(omegas), size(array) )
        omegas_curr = omegas(ii, jj, :);
    else
        omegas_curr = omegas;
    end
    
    if kk > 1 && track_fitting
        % find nearest previously fit point
        dists = dist_fn( xs(jj), ys(ii), xx, yy);
        [~, dist_ind] = sort( dists(:) );
        candidates = dist_ind .* already_fit( dist_ind(:) )';
        inds = candidates( candidates ~= 0);
        [ii_prev, jj_prev] = ind2sub( size(already_fit), inds(1));
        
        fit_ref(ii, jj) = inds(1);
        init_params = squeeze(fit_params(ii_prev, jj_prev, :))';
    end

    [fitp, ~, ffh, se, ~] = fit1D(omegas_curr, squeeze(array(ii, jj, :)), 1 ./ squeeze(array_unc(ii, jj, :)) .^2, fn_name, init_params, fixed_params, lbs, ubs);

    % sort fit parameters by center
    [centers, I] = sort( fitp(1:4:end) );
    hwhm = fitp(2:4:end);
    hwhm = hwhm(I);
    amps = fitp(3:4:end);
    amps = amps(I);
    fitp(1:4:end) = centers;
    fitp(2:4:end) = hwhm;
    fitp(3:4:end) = amps;

    % sort uncertainty by center
    center_se = se(1:4:end);
    center_se = center_se(I);
    hwhm_se = se(2:4:end);
    hwhm_se = hwhm_se(I);
    amps_se = se(3:4:end);
    amps_se = amps_se(I);
    se(1:4:end) = center_se;
    se(2:4:end) = hwhm_se;
    se(3:4:end) = amps_se;

    %
    fit_params(ii, jj, :) = fitp;
    std_errs(ii, jj, :) = se;
    already_fit(ii, jj) = 1;
end
           
if array_2d
    fit_params = permute(fit_params, [2, 3, 1]);
    std_errs = permute(std_errs, [2, 3, 1]);
    fit_ref = permute(fit_ref, [2, 1]);
end
    
end