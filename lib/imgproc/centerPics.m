function [img_stack_centered, roi_start_x, roi_start_y] = centerPics(img_stack, cxs, cys)
% centerPics Centers a stack of pictures around their individual centers 
% cxs and cys. Centered pictures are the same size as the initial pictures,
% and may be padded with zeros.
%
%   TODO: deprecate in favor of get_centered_img ?
%
%   Arguments:
%   ----------
%
%   [img_stack_centered, roi_start_x, roi_start_y] = centerPics(img_stack, cxs, cys)
%
%   img_stack: an Ny x Nx x NImgs array of images. This function returns 
%   
%   cxs: an NImgs list of center positions. If none is provided, use the
%   center of mass of each picture as the center.
%
%   cys:
%
%   Returns:
%   ----------
%
%   img_stack_centered: an array of the same size which is centered on the
%   combined center of mass of the image stack. 
%
%   roi_start_x and roi_start_y: give the
%   coordinate of the upper left most point in CenteredPics in
%   terms of the initial array, Pics. This coordinate is not
%   necessarily in the initial array, which is treated as if it
%   is padded by zeros.

if ~exist('cxs', 'var') || isempty(cxs)
    %get array of centers
    [cxs, ~, ~] = get_moment(img_stack, 1, [], []);
end

if ~exist('cys', 'var') || isempty(cys)
    [~, cys, ~] = get_moment(img_stack, 1, [], []);
end

if ~isequal( length(cxs), length(cys) ) && length(cxs) == size(img_stack, 3)
    error('cxs and cys must be equal in length and the same as the number of images img_stack');
end

% center coordinates
% ix = (size(img_stack, 2) - 1) / 2; % 
% iy = (size(img_stack, 1) - 1) / 2;
% think of coordinate centered on each pixel.
% the goal is to minimize the difference between ix and cx
ix = (size(img_stack, 2) + 1) / 2;
iy = (size(img_stack, 1) + 1) / 2;
% xshifts = floor(ix - cxs);
% yshifts = floor(iy - cys);
xshifts = round(ix - cxs);
yshifts = round(iy - cys);



nimgs = size(img_stack,3);
img_stack_centered = zeros( size(img_stack) );
for ii = 1:nimgs
    try
        pixels_to_move_x = xshifts(ii);
        pixels_to_move_y = yshifts(ii);
        img_stack_centered(:, :, ii) = getShiftedMat(img_stack(:, :, ii), pixels_to_move_x, pixels_to_move_y, 0, 0);
    catch
        img_stack_centered(:, :, ii) = NaN(size(img_stack,1), size(img_stack,2));
    end
end

%remove failed images from stack
NanImages = squeeze( all(all(isnan(img_stack_centered), 1), 2) );
img_indices_to_remove = find(NanImages);
img_stack_centered(:, :, img_indices_to_remove) = [];

xshifts(img_indices_to_remove) = [];
yshifts(img_indices_to_remove) = [];

roi_start_x = 1 - xshifts;
roi_start_y = 1 - yshifts;
%TODO decide what to do with XShifts/YShifts...
% img_stack_centered = img_stack;
end