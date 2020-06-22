function [array_quad, xs, ys] = get_quadrant(array, quadrant_direction, mode, xs_in, ys_in)
 % Get a quadrant of a 2D, 3D, or 4D array based on the first two dimensions
 % of the array
 %
 % quadrant_direction is 'northeast', 'northwest', 'southeast', or
 % 'southwest', and specifies which quadrant of a matrix to extract.
 %
 % mode is 'even only', 'include_center', 'exclude_center'. Note that
 % 'include center' and 'exclude center' are only meaninful for odd sized
 % matrices, where the matrix cannot be evenly divided into fourths.

% check direction argument
if ~exist('quadrant_direction', 'var')
    quadrant_direction = 'northeast';
end

quad_dirs = {'northeast', 'northwest', 'southeast', 'southwest'};
fn = @(c) strcmp(c, quadrant_direction);
if ~any(cellfun(fn, quad_dirs))
    error('quadrant direction argument was not one of the allowed keywords.');
end

% check mode argument

if ~exist('mode', 'var')
    mode = 'include_center';
end

allowed_modes = {'even only', 'include_center', 'exclude_center'};
fn = @(c) strcmp(c, mode);
if ~any(cellfun(fn, allowed_modes))
    error('mode argument was not one of the allowed keywords.');
end

% check size of array
xsize = size(array, 2);
ysize = size(array, 1);

if strcmp(mode, 'even only') && ( mod(ysize, 2) == 0 || mod(xsize, 2) == 0)
    error('ysize and xsize must be odd');
end

if ~exist('xs_in', 'var') || isempty(xs_in)
    xs_in = 1 : xsize;
end

if ~exist('ys_in', 'var') || isempty(ys_in)
    ys_in = 1 : ysize;
end

% deal with even and odd cases
y1 = 1;
y4 = ysize;
if mod(ysize, 2) == 0
    y2 = ysize/2;
    y3 = ysize/2 + 1;
else
    if strcmp(mode, 'include_center')
        y2 = (ysize - 1) / 2 + 1;
        y3 = y2;
    elseif strcmp(mode, 'exclude_center')
        y2 = (ysize - 1) / 2;
        y3 = y2 + 2;
    else
        error();
    end
end

x1 = 1;
x4 = xsize;
if mod(xsize, 2) == 0
    x2 = xsize/2;
    x3 = xsize/2 + 1;
else
   if strcmp(mode, 'include_center')
       x2 =  (xsize - 1) / 2 + 1;
       x3 = x2;
   elseif strcmp(mode, 'exclude_center')
        x2 = (xsize - 1) / 2;
        x3 = x2 + 2;
   else
       error();
   end
end

if strcmp(quadrant_direction, 'northeast')
    x_start = x3;
    x_end = x4;
    y_start = y1;
    y_end = y2;
elseif strcmp(quadrant_direction, 'northwest')
    x_start = x1;
    x_end = x2;
    y_start = y1; 
    y_end = y2;
elseif strcmp(quadrant_direction, 'southeast')
    x_start = x3;
    x_end = x4;
    y_start = y3; 
    y_end = y4;
elseif strcmp(quadrant_direction, 'southwest')
    x_start = x1;
    x_end = x2;
    y_start = y3; 
    y_end = y4;
else
    error();
end

if ndims(array) == 2
%     array_quad = array(center_y:end, center_x:end);
    array_quad = array(y_start:y_end, x_start:x_end);
elseif ndims(array) == 3
     array_quad = array(y_start:y_end, x_start:x_end, :);
elseif ndims(array) == 4
     array_quad = array(y_start:y_end, x_start:x_end, :, :);
else
    error('unsupported number of array dimensions');
end

xs = xs_in(x_start:x_end);
ys = ys_in(y_start:y_end);

end


