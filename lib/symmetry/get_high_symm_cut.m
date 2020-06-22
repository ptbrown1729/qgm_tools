function [array_gxmg, array_indices_gxmg,...
          array_gymg, array_indices_gymg, linear_index] = ...
        get_high_symm_cut(array, mode, qxs, qys)
%[array_gxmg, array_indices_gxmg, array_gymg, array_indices_gymg, linear_index] = ...
%         get_high_symm_cut(array, mode)
%
% Given a square array, take the edges and diagonal. If instead you want
% the center and half the diagonal, use get_quadrant to first extract a
% quadrant of your matrix.
%
% Arguments:
% --------------
%
% array: A 2D, 3D, 4D, or 5D array where the first two dimensions represent
% equally spaced positions in a square. This function takes this 2D
% representation and produces a 1D representation by only taking the
% elements of the array along points of high symmetry.
%
% mode: can be 'lowerleft', which is default and takes the Gamma point as
% the lower left corner, the X point as the lower right corner, and the 
% M point as the upper right corner. Here 'lower' means smaller subscript
% index. If mode is 'upperleft', the Gamma point is the upper left, the X
% point is the upper right, and the M point is the lower right.
%
% qxs: 1d array of qx points. If not provided, assume these are 0,...pi
%
% qys: 1d array of qy points. If not provided, assume these are 0,...pi
%
% Returns:
% --------------
%
% array_gxmg: is the array along the Gamma-X-M-Gamma path, where is the
% lower-left corner, X is the lower-right corner, M is the upper right
% corner.
%
% array_indices_gxmg: give the indices in the first two dimensions of array
% for each element of array_gxmg. These are given as 1D indices which 
% can be transformed into subscripts with the matlab function ind2sub
%
% array_gymg: is the array along the Gamma-Y-M-Gamma path, where Y is the
% upper left corner.
%
% array_indices_gxmg: give the indices in the first two dimensions of array
% for each element of array_gxmg. These are given as 1D indices which 
% can be transformed into subscripts with the matlab function ind2sub
%
% linear_index: indicates where the given points are along each of these
% lines. Gamma-X maps to the interval [0,1], X-M to [1, 2], and M-Gamma to
% [2, 3].
%

% argument checking
if ~exist('mode', 'var')
    mode = 'lowerleft';
end

allowed_modes = {'lowerleft', 'upperleft'};
fn = @(c) strcmp(c, mode);
if ~any(cellfun(fn, allowed_modes))
    error('mode argument was %s, which is not one of the allowed keywords.', mode);
end

% array size checking
xsize = size(array, 2);
ysize = size(array, 1);

if ysize ~= xsize
    error('ysize and xsize must be equal');
end

% TODO: wouldn't it be better to simply flip array if the other amode is
% called?
% y and x indices at high symmetry points
if strcmp(mode, 'lowerleft')
    y_gamma = ysize;
    x_gamma = 1;

    y_xpt = y_gamma;
    x_xpt = xsize;

    y_ypt = 1;
    x_ypt = x_gamma;

    y_m = y_ypt;
    x_m = x_xpt;
elseif strcmp(mode, 'upperleft')
    y_gamma = 1;
    x_gamma = 1;

    y_xpt = y_gamma;
    x_xpt = xsize;

    y_ypt = ysize;
    x_ypt = x_gamma;

    y_m = y_ypt;
    x_m = x_xpt;
end

% if qxs and qys are not supplied, use default value
if ~exist('qxs', 'var') || isempty(qxs) || ~exist('qys', 'var') || isempty(qys)
    if strcmp(mode, 'lowerleft')
        qxs = linspace(0, pi, xsize);
        qys = linspace(pi, 0, ysize);
    elseif strcmp(mode, 'upperleft')
        qxs = linspace(0, pi, xsize);
        qys = linspace(0, pi, ysize);
    else
        error();
    end 
end

[qxqx, qyqy] = meshgrid(qxs, qys);

% get subscripts for array slicing
% xs_gx = x_gamma : x_xpt;
% Trick to get x_gamma : 1 : x_xpt or x_gamma : -1 : x_xpt depending on which is larger 
xs_gx = linspace(x_gamma, x_xpt, abs(x_xpt - x_gamma) + 1);
ys_gx = y_gamma * ones( size(xs_gx) );

% ys_xm = y_xpt - 1 : -1 : y_m;
ys_xm = linspace(y_xpt, y_m, abs(y_xpt - y_m) + 1);
ys_xm = ys_xm(2:end);
xs_xm = x_m * ones( size(ys_xm) );

% ys_gy = y_gamma : -1 : y_ypt;
ys_gy = linspace(y_gamma, y_ypt, abs(y_ypt - y_gamma) + 1);
xs_gy = x_ypt * ones( size(ys_gy) );

% xs_ym = x_ypt + 1 : x_m;
xs_ym = linspace(x_ypt, x_m, abs(x_ypt - x_m) + 1);
xs_ym = xs_ym(2:end);
ys_ym = y_m * ones( size(xs_ym) );

ys_mg = linspace(y_m, y_gamma, abs(y_m - y_gamma) + 1);
ys_mg = ys_mg(2:end);
xs_mg = linspace(x_m, x_gamma, abs(x_m - x_gamma) + 1);
xs_mg = xs_mg(2:end);

% get qx, qy points
qx_gx = qxqx(sub2ind([size(qxqx, 1), size(qxqx, 2)], ys_gx, xs_gx));
qx_xm = qxqx(sub2ind([size(qxqx, 1), size(qxqx, 2)], ys_xm, xs_xm));
qx_gy = qxqx(sub2ind([size(qxqx, 1), size(qxqx, 2)], ys_gy, xs_gy));
qx_ym = qxqx(sub2ind([size(qxqx, 1), size(qxqx, 2)], ys_ym, xs_ym));
qx_mg = qxqx(sub2ind([size(qxqx, 1), size(qxqx, 2)], ys_mg, xs_mg));

qy_gx = qyqy(sub2ind([size(qyqy, 1), size(qyqy, 2)], ys_gx, xs_gx));
qy_xm = qyqy(sub2ind([size(qyqy, 1), size(qyqy, 2)], ys_xm, xs_xm));
qy_gy = qyqy(sub2ind([size(qyqy, 1), size(qyqy, 2)], ys_gy, xs_gy));
qy_ym = qyqy(sub2ind([size(qyqy, 1), size(qyqy, 2)], ys_ym, xs_ym));
qy_mg = qyqy(sub2ind([size(qyqy, 1), size(qyqy, 2)], ys_mg, xs_mg));

% convert subscripts to indices for full array
indices_gx = sub2ind([size(array, 1), size(array, 2)], ys_gx, xs_gx);
indices_xm = sub2ind([size(array, 1), size(array, 2)], ys_xm, xs_xm);
indices_gy = sub2ind([size(array, 1), size(array, 2)], ys_gy, xs_gy);
indices_ym = sub2ind([size(array, 1), size(array, 2)], ys_ym, xs_ym);
indices_mg = sub2ind([size(array, 1), size(array, 2)], ys_mg, xs_mg);

% TODO: is there a way to write this for arbitrary dimensions?
if ndims(array) == 2
    array_gx = transpose(array(indices_gx));
    array_xm = transpose(array(indices_xm));
    array_gy = transpose(array(indices_gy));
    array_ym = transpose(array(indices_ym));
    array_mg = transpose(array(indices_mg));

elseif ndims(array) == 3
    array_perm = permute(array, [3, 1, 2]);
    array_gx = permute(array_perm(:, indices_gx), [2, 1]);
    array_xm = permute(array_perm(:, indices_xm), [2, 1]);
    array_gy = permute(array_perm(:, indices_gy), [2, 1]);
    array_ym = permute(array_perm(:, indices_ym), [2, 1]);
    array_mg = permute(array_perm(:, indices_mg), [2, 1]);

elseif ndims(array) == 4
    array_perm = permute(array, [3, 4, 1, 2]);
    array_gx = permute(array_perm(:, :, indices_gx), [3, 1, 2]);
    array_xm = permute(array_perm(:, :, indices_xm), [3, 1, 2]);
    array_gy = permute(array_perm(:, :, indices_gy), [3, 1, 2]);
    array_ym = permute(array_perm(:, :, indices_ym), [3, 1, 2]);
    array_mg = permute(array_perm(:, :, indices_mg), [3, 1, 2]);
    
elseif ndims(array) == 5
    array_perm = permute(array, [3, 4, 5, 1, 2]);
    array_gx = permute(array_perm(:, :, :, indices_gx), [4, 1, 2, 3]);
    array_xm = permute(array_perm(:, :, :, indices_xm), [4, 1, 2, 3]);
    array_gy = permute(array_perm(:, :, :, indices_gy), [4, 1, 2, 3]);
    array_ym = permute(array_perm(:, :, :, indices_ym), [4, 1, 2, 3]);
    array_mg = permute(array_perm(:, :, :, indices_mg), [4, 1, 2, 3]);
    
else
    error('get_high_symm_cut not implemented for arrays with 6 or more dimensions');
end        

% map the interval [a, b] to [c, d] using
% x -> (x-a)/(b-a) * (d - c) + c;
% map interval [1, size(array_gx)] -> [0, 1]
% linear_index_gx =  ( (1 : size(array_gx, 1)) - 1) / (size(array_gx, 1) - 1) + 0;
[linear_index_gx, ~] = kvect2index(qx_gx, qy_gx);
% map interval [1, size(array_xm)] -> [1 + di, 2]
% di = 1 / size(array_xm, 1);
% linear_index_xm = ( (1 : size(array_xm, 1)) - 1) / (size(array_xm, 1) - 1) * (1 - di) + (1 + di);
[linear_index_xm, ~] = kvect2index(qx_xm, qy_xm);
% map interval [1, size(array_mg)] -> [2 + di, 3]
% di = 1 / size(array_mg, 1);
% linear_index_mg = ( (1 : size(array_mg, 1)) - 1) / (size(array_mg, 1) - 1) * (1 - di) + (2 + di);
[linear_index_mg, ~] = kvect2index(qx_mg, qy_mg);
if linear_index_mg(end) == 0
    linear_index_mg(end) = 3;
end

array_gxmg = cat(1, array_gx, array_xm, array_mg);
array_indices_gxmg = horzcat(indices_gx, indices_xm, indices_mg);
array_gymg = cat(1, array_gy, array_ym, array_mg);
array_indices_gymg = horzcat(indices_gy, indices_ym, indices_mg);
linear_index = cat(2, linear_index_gx, linear_index_xm, linear_index_mg);

end