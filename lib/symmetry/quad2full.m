function [full_mat, xs, ys] = quad2full(quad, quadrant_direction, includes_center, xs_in, ys_in)
% assemble a full image from a single quadrant.
%
% quad: the array (or stack of arrays) used to construct the larger arrays
%
% quadrant_direction: indicates which quadrant 'quad' represents. Should be
% 'northeast', 'northwest', 'southeast', or 'southwest'
%
% includes_center: boolean, whether or not quad includes the center
% row/column of the matrix. If it does, this will not be duplicated. TODO:
% for a rectangular matrix it is possible that we include either of the x
% or y center but not the other. Could handle that case.

if ~exist('includes_center', 'var') || isempty(includes_center)
    includes_center = 0;
end

if ~exist('quadrant_direction', 'var') || isempty(quadrant_direction)
    quadrant_direction = 'southeast';
end

if ~exist('xs_in', 'var') || isempty(xs_in)
    xs_in = 1 : size(quad, 2);
    if includes_center
        xs_in = xs_in - 1;
    end
end

if ~exist('ys_in', 'var') || isempty(ys_in)
    ys_in = 1 : size(quad, 1);
    if includes_center
        ys_in = ys_in - 1;
    end
end

quad_dirs = {'northeast', 'northwest', 'southeast', 'southwest'};
fn = @(c) strcmp(c, quadrant_direction);
if ~any(cellfun(fn, quad_dirs))
    error('quadrant direction argument was not one of the allowed keywords.');
end

if size(quad, 3) > 1
    npics = size(quad, 3);

    % combine quadrants for first image to get required size of matrix
    [first_mat, xs, ys] = quad2full(quad(:, :, 1), quadrant_direction, includes_center, xs_in, ys_in);
    
    % array to hold results
    full_mat = zeros(size(first_mat, 1), size(first_mat, 2), npics);
    full_mat(:, :, 1) = first_mat;
    
%     xs = zeros( length(xs_first), npics);
%     ys = zeros( length(ys_first), npics);
%     xs(:, 1) = xs_first;
%     ys(:, 1) = ys_first;
    
    % do other matrices
    for ii = 2 : npics
        full_mat(:, :, ii) = quad2full(quad(:, :, ii), quadrant_direction, includes_center, xs_in, ys_in);
    end
    
else
    % write function for southeast quadrant, and if we have a different
    % quadrant first transform it to that
    xs = xs_in;
    ys = ys_in;
    if strcmp(quadrant_direction, 'northeast')
        quad = flip(quad, 1);
        ys = -flip(ys);
    elseif strcmp(quadrant_direction, 'northwest')
        quad = flip(flip(quad, 1), 2);
        xs = -flip(xs);
        ys = -flip(ys);
    elseif strcmp(quadrant_direction, 'southeast')
    elseif strcmp(quadrant_direction, 'southwest')
        quad = flip(quad, 2);
        xs = -flip(xs);
    end

    quad_flipx = flip(quad, 2);

    if includes_center
        full_mat = horzcat(quad_flipx, quad(:, 2:end)) ;
        full_mat = vertcat(flip(full_mat, 1), full_mat(2:end, :) );
        
        xs = horzcat( - flip(xs), xs(2:end) );
        ys = horzcat( - flip (ys), ys(2:end) );
    else
        full_mat = horzcat(quad_flipx, quad);
        full_mat = vertcat( flip(full_mat, 1), full_mat); 
        
        xs = horzcat( - flip(xs), xs);
        ys = horzcat( - flip(ys), ys);
    end
end

end