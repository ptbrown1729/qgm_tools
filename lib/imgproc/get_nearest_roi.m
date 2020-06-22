function [x_end_pts, x_size_binned, y_end_pts, y_size_binned] = ...
                get_nearest_roi(x_end_pts, y_end_pts, img_bin_size, mode)
%   Given an initial region of interest, compute the nearest larger region
%   of interest which will be evenly divisible into some number of bins.
%
%   arguments
%   -----------------
%   x_end_pts: [x_start, x_end] approximate start and end point of ROI
%
%   y_end_pts: [y_start, y_end] approximate start and end point of the roi
%
%   img_bin_size: desired size to bin the final image
%
%   mode: 'even', 'odd', or 'any'. If 'even', then forces the ROI to
%   contain an even number of bins of img_bin_size. If 'odd', forces an odd
%   number of bins. If 'any', then does not stipulate an even or odd number
%   of bins.
%
%   results
%   -----------------
%   x_end_pts:
%
%   x_size_binned:
%
%   y_end_pts:
%
%   y_size_binned:

if ~exist('mode', 'var') || isempty(mode)
    mode = 'even';
end

allowed_modes = { 'even', 'odd', 'any'};
fn = @(c) strcmp(c, mode);
if ~any(cellfun(fn, allowed_modes))
    error('mode direction argument was not one of the allowed keywords.');
end


x_start = round( x_end_pts(1) );
x_end = round( x_end_pts(2) );
y_start = round( y_end_pts(1) );
y_end = round( y_end_pts(2) );
% x_start = obj.cx - obj.x_bz_edge;
% x_end = obj.cx + obj.x_bz_edge;
% y_start = obj.cy - obj.x_bz_edge;
% y_end = obj.cy + obj.x_bz_edge;

% ensure image cropping coordinates give number divisible by bin size
% TODO: could also pick binning to have a well defined center coordinate,
% in that case need to enforce that the image size is an odd multiple of
% the bin size.

% add pixels to make the cropped image divisible by the desired bin size
x_size_init = x_end - x_start + 1;
x_to_add = img_bin_size - mod(x_size_init, img_bin_size);

% check if we have an even or odd multiple of our bin size
x_even_multiple = 1 - mod ( (x_size_init + x_to_add) / img_bin_size, 1 );


if (x_even_multiple && strcmp(mode, 'odd')) || (~x_even_multiple && strcmp(mode, 'even'))
% if our image is currently an odd multiple of the bin size, add an extra
% pixel
    x_to_add = x_to_add + img_bin_size;
end

% compute new coordinates by including extra pixels to either side of our
% image.
x_start = x_start - floor(x_to_add / 2);
x_end = x_end + ceil(x_to_add / 2);
x_size_final = x_end - x_start + 1;
x_size_binned = x_size_final / img_bin_size;

x_end_pts = [x_start, x_end];

% same for y direction
y_size_init = y_end - y_start + 1;
y_to_add = img_bin_size - mod( y_size_init, img_bin_size );

% check if we have an even or odd multiple of our bin size
y_even_multiple = 1 - mod ( (y_size_init + y_to_add) / img_bin_size, 1);

if (y_even_multiple && strcmp(mode, 'odd')) || (~y_even_multiple && strcmp(mode, 'even'))
    y_to_add = y_to_add + img_bin_size;
end

y_start = y_start - floor(y_to_add / 2);
y_end = y_end + ceil(y_to_add / 2);
y_size_final = y_end - y_start + 1;
y_size_binned = y_size_final / img_bin_size;

y_end_pts = [y_start, y_end];

end