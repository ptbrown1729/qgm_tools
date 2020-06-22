function [ array_eighth ] = get_eighth( array )
%GET_EIGHTH Summary of this function goes here
%   Detailed explanation goes here

% get unique part of brillouin zone for d4 symmetry
% (kx, ky) in [0, pi] x [0, pi] and kx > ky

xsize = size(array, 2);
ysize = size(array, 1);

if mod(ysize, 2) == 0 || mod(xsize, 2) == 0
    error('ysize and xsize must be odd');
end

if ysize ~= xsize
    error('ysize and xsize must be equal');
end

center_x = (xsize - 1) / 2 + 1;
center_y = (ysize - 1) / 2 + 1;

x_indices = [];
y_indices = [];

for ii = center_x : size(array, 1)
    x_indices = horzcat(x_indices, ii : size(array, 1) );
    y_indices = horzcat(y_indices, ii * ones(1, length(ii : size(array, 1)) ));
end

indices_eight = sub2ind([size(array, 1), size(array, 2)],...
                        y_indices, x_indices);

if ndims(array) == 2
    array_eighth = array(indices_eighth);
elseif ndims(array) == 3
    array_eighth = permute(array, [3, 1, 2]);
    array_eighth = permute(array_eighth(:, indices_eight), [2, 1]);
elseif ndims(array) ==4
    array_eighth = permute(array, [3, 4, 1, 2]);
    array_eighth = permute(array_eighth(:, :, indices_eighth), [3, 1, 2]);
else
    error();
end
end

