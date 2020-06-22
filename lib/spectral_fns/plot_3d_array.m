function fig_handle = plot_3d_array(array, array_unc, xs, ys, zs, ref_vals, fig_handle, use_common_scale, str_spec)
% array: 3d array where the first dimension is represents the y-position,
% the second dimension represents the x-position, and the third dimension
% represents some other variable (e.g. energy, rf frequency, etc.)
%
% array_unc: 3d array the same size as array representing the uncertainty
% associated with the points in the array. If empty, ignored.
%
% xs: x-coordinates, 1d array with length equal to size(array, 2)
%
% ys: y-coordinates, 1d array with length equal to size(array, 1)
%
% zs: z-coordinates, either a 1d array with length equal to size(array, 3)
%     or a 3D array of the same size as array
%
% ref_vals: 2D array of the same size as the first two dimensions of array.
% Each value is a reference line which will be plotted 
%
% fig_handle: figure on which to plot the array. Supplying a figure is
% useful if plotting several arrays of the same size.
%
% use_common_scale: whether or not to force same scale for all plots
%
% str_spec: plotting specifier to be used with all arrays.

if ndims(array) == 2
    array = permute(array, [3, 1, 2]);
end

if exist('array_unc', 'var') && ~isempty(array_unc) && ndims(array_unc) == 2
    array_unc = permute(array_unc, [3, 1, 2]);
end

if ~exist('xs', 'var') || isempty(xs)
    xs = 1 : size(array, 2);
end

if ~exist('ys', 'var') || isempty(ys)
    ys = 1 : size(array, 1);
end

if ~exist('zs', 'var') || isempty(zs)
    zs = 1 : size(array, 3);
end

if ~exist('fig_handle', 'var') || isempty(fig_handle)
    fig_handle = figure;
end

if ~exist('use_common_scale', 'var') || isempty(use_common_scale)
    use_common_scale = 1;
end

if ~exist('str_spec', 'var') || isempty(str_spec)
    str_spec = 'b-';
end

if ~exist('ref_vals', 'var')
    ref_vals = '';
end

if ~isempty(ref_vals) && numel(ref_vals) == 1
    ref_vals = ref_vals * ones( size(array(:, :, 1)) );
end

is_2d = length(xs) > 1 && length(ys) > 1;
if is_2d
    nrows = numel(ys);
    ncols = numel(xs);
else
    nplots = max( numel(xs), numel(ys) );
    nrows = floor( sqrt(nplots / 1.5) );
    ncols = ceil( nplots / nrows);
end
    
min_val = min( array(:) );
max_val = max( array(:) );
y_limits = [min_val - 0.2 * abs(min_val), max_val + 0.2 * abs(max_val) ];

for ii = 1 : numel(ys)
    for jj = 1 : numel(xs)
        
        if isequal( size(zs), size(array))
            z_curr = squeeze( zs(ii, jj, :) );
        else
            z_curr = zs;
        end
        
        if is_2d
            plot_index = (ii-1) * ncols + jj;
        else
            plot_index = max(ii, jj);
        end
        
        ax = subplot(nrows, ncols, plot_index);
        hold on;
        
        if ~isempty(array_unc)
            errorbar(z_curr, squeeze(array(ii, jj, :)), squeeze(array_unc(ii, jj, :)), str_spec);
        else
            plot(z_curr, squeeze(array(ii, jj, :)), str_spec);
        end
        
        if ~isempty(ref_vals)
            plot([ref_vals(ii, jj), ref_vals(ii, jj)], [ax.YLim(1), ax.YLim(2)]);
        end
        
        if use_common_scale
            ax.YLim = y_limits;
        end
        
        ax.XTick = '';
        ax.YTick = '';
        
%         title( sprintf('(%0.2f, %0.2f)', xs(jj), ys(ii)));
    
    end
end


end