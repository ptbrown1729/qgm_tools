function [centered_img, roi_start_x, roi_start_y, cx_roi, cy_roi] = get_centered_img(img, cx, cy, stable_size)
%   get_centered_img centers an image or stack of images around a given
%   center as closely as possible
%
%   [centered_img, roi_start_x, roi_start_y, cx_roi, cy_roi] ...
%                               = get_centered_img(img, cx, cy, stable_size)
%
%   arguments:
% -----------------
%   img: A 2d matrix or a 3d stack of images.
%
%   cx: Center position in the x-dimension (second dimension). cx must be a
%   single number. 
%
%   cy: Center position in the y-dimension (first dimension)
%
%   stable_size: boolean value. If 1, centered_img will be the same size as
%   img. If 0, centered_img may be 1 pixel smaller than img in either or
%   both dimensions
%
%   outputs:
% -----------------
%   centered_img: img after being shifted so that cx and cy are as near the
%   geometric center as possible. Note that cx and cy must be the same for
%   all images because if stable_size=0 then the output images are not
%   guaranteed to be the same size for different choices of cx and cy.
%
%   roi_start_x: 
%   roi_start_y: centered_img(j, i) = img(j + roi_start_y - 1, i + roi_start_x - 1 )
%   These variables are named 'roi_start' because centered_img(1, 1) =
%   img(roi_start_y, roi_start_x). Note that if either of these coordinates
%   fall outside of img, we treat img as if it is padded with zeros.
%
%   cx_roi: cx in the coordinates of centered_img
%
%   cy_roi: cy in the coordinates of centered_img

if ~exist('cx', 'var') || isempty(cx)
    [cx, ~, ~] = get_moment(img, 1);
    cx = mean(cx(:));
end

if ~exist('cy', 'var') || isempty(cy)
    [~, cy, ~] = get_moment(img, 1);
    cy = mean(cy(:));
end

if ~exist('stable_size', 'var') || isempty(stable_size)
    stable_size = 0;
end

% first determine how much we need to shift matrix
% get geometric center
ix = (size(img, 2) + 1) / 2;
iy = (size(img, 1) + 1) / 2;
% distance of geometric center from actual center
xshift = round(ix - cx);
yshift = round(iy - cy);

centered_img = getShiftedMat(img, xshift, yshift, 0, 0);

% start of region of interest
roi_start_x = 1 - xshift;
roi_start_y = 1 - yshift;

if ~stable_size
    % determine if we need to crop matrix
    % if ix and cx are both close to the center of pixels or both near the edge
    % of pixels, no need to crop. If one is near the center of a pixel and the
    % other is near the edge of a pixel, we can improve our centering by
    % cropping.
    crop_parity_x = round( 2 * (ix - cx) );
    crop_parity_y = round( 2 * (iy - cy) );

    % better to crop pixel at the end?
    % think it doesn't matter because would crop from side that will be lost in
    % shift anyways.
    if mod(crop_parity_x, 2) == 1
        % if img is bettered centered after cropping one pixel, crop
        if crop_parity_x > 0
            centered_img = centered_img(:, 1 : end - 1, :);
        else
            centered_img = centered_img(:, 2 : end, :);
            roi_start_x = roi_start_x + 1;
        end
    end

    if mod(crop_parity_y, 2) == 1
        if crop_parity_y > 0
            centered_img = centered_img(1 : end - 1, :, :);
        else
            centered_img = centered_img(2 : end, :, :);
            roi_start_y = roi_start_y + 1;
        end
    end
end

cx_roi = cx - (roi_start_x - 1);
cy_roi = cy - (roi_start_y - 1);

end