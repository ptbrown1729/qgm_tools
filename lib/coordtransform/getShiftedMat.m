function [shifted_mat, roi_start_x, roi_start_y] = getShiftedMat(mat, dx, dy, wrap, pad_value)
%getShiftedMat Translates an nD array along its first two
%dimensions by distance dy and dx respectively.
%
%   [shifted_mat] =
%   getShiftedMat(mat, dx, dy, wrap) 
%
%   Arguments:
%   -------------------
%   mat: the nD array to be shifted. We can think of this as having a fixed
%   field of view which selects some portion of an infinite matrix. 
%
%   dx: distance to shift the array along its x-dimension (2nd dimension).
%   A positive value for dx shifts the matrix in the positive x-direction.
%   If dx is a matrix, will return an array of shifted matrices for each
%   dx.
%
%   dy: distance to shift the array long its y-dimension (1st dimension).
%   A positive value shifts the matrix in the positive y-direction (i.e.
%   downward if you display the matrix using imagesc). dy must be the same
%   size as dx.
%
%   wrap: If = 0, the matrix is treated as if it is padded with pad_values.
%   If wrap = 1, a circular shift is performed.
%
%   pad_value: the matrix is treated as though all values beyond its limits
%   are this value
%
%
%   Returns:
%   ---------------
%   shifted_mat: The shifted matrix, defined by shifted_mat(i,j,...) =
%   mat(i - dy, j - dx, ...). If dx is an array, shifted_mat has size
%   [size(mat), size(dx)].
%
%see also circshift
%
%   TODO: allow to accept multiple Dxs and Dys (Dxs and Dys
%   matrices of the same shape). Return all requested shifts as
%   extra dimensions. So if Dx and Dy are KxL matrices and Mat is
%   NxM, would return a NxMxKxL matrix. ???

if ~exist('wrap', 'var') || isempty(wrap)
    wrap = 0;
end

if numel(dx) ~= numel(dy)
    error('dx and dy must be same size');
end

if ~exist('pad_value', 'var') || isempty(pad_value)
    pad_value = 0;
end

if numel(dx) > 1
    %if dx is a matrix, call this function on each element of dx.
    shifted_mat = zeros(size(mat), size(dx));
    roi_start_x = zeros(size(dx));
    roi_start_y = zeros(size(dy));
    
    for ii = 1:numel(dx)
        [shifted_mat(:, :, ii), roi_start_x(ii), roi_start_y(ii)] = ...
                        getShiftedMat(mat, dx(ii), dy(ii), wrap, pad_value);
    end
    
else
    if wrap
        shifted_mat = circshift(mat, [dy, dx]);
    else
%         shifted_mat = zeros(size(mat));
        shifted_mat = pad_value * ones(size(mat));
        
        % TODO: why is this here instead of at the start?
        if ~isnumeric(mat)
            disp('First argument not numeric')
        end
        
        if dx >= 0 && dy >= 0
            shifted_mat(dy + 1 : end, dx + 1 : end, :) = ...
                mat(1 : end - abs(dy), 1 : end - abs(dx), :);
        elseif dx >= 0 && dy < 0
            shifted_mat(1 : end - abs(dy), dx + 1 : end, :) = ...
                mat(abs(dy) + 1 : end, 1 : end - abs(dx) ,:);
        elseif dx < 0 && dy >= 0
            shifted_mat(dy + 1 : end, 1 : end - abs(dx), :) = ...
                mat(1 : end - abs(dy), abs(dx) + 1 : end, :);
        elseif dx < 0 && dy < 0
            shifted_mat(1 : end - abs(dy), 1 : end - abs(dx), :) = ...
                mat(abs(dy) + 1 : end, abs(dx) + 1 : end, :);
        end
        roi_start_x = 1 - dx;
        roi_start_y = 1 - dy;
    end
end
end