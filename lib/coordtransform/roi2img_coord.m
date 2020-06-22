function [x_img, y_img, affine_mat] = roi2img_coord(transform_params, x_roi, y_roi)
% Get coordinates in the full-image space given coordinates in the
% region-of-interest space and transformation parameters which specify the
% choice of region of interest and choice of binning.
%
% amiguities of up to roughly nbin_x or nbin_y pixels can be created if one
% does not know the exact roi start point. For example, if you don't know
% if the bins were determined beginning from the upper-left or the
% lower-right, etc.
%
% We adopt the viewpoint that pixel M(i, j) occupies the coordinates 
% X x Y = [j-1, j] x [i-1, i], which means that the pixel is labelled by
% its lower-right hand coordinate and centered at (x,y) = (j-0.5,i-0.5). 
% This choice of coordinates simplifies thinking about binning and allows
% the upper-left hand corner of the image to be the origin. 
%
% Compare this convention with centering the pixels at (x,y) = (j,i),
% which causes the origin to be at (0.5, 0.5). A major disadvantage of this
% latter choice is that it leads to non-intuitive coordinate
% transformations when doing binning. For example, suppose you want to bin
% your image 2x2. Look at the first square of pixels associated with
% indices (1, 1); (1, 2); (2, 1); and (2, 2).
% Using pixels [j-1, j] x[i-1, i] these would transform to
% (0.5, 0.5); (0.5, 1); (1, 0.5); and (1, 1) as we might expect. And this
% makes it easy to determine which of the binned/region-of-interest pixels
% these full-space pixels are in. We can simply take the ceiling of the
% pixel coordinates. This is all sensible, because the ceil function is
% asking for the closest integer to the lower-right direction.
%
% If we use instead [j-0.5, j+0.5] x [j-0.5, j+0.5] these become
% (0.75, 0.75); (0.75, 1.25); (1.25, 0.75); (1.25, 1.25), which are much
% less intuitive. However, determining which binned pixel is still easy.
% Now instead of using the ceiling function, we would need to use the 
% round function.
%
% transform_params = [x_roi_start, y_roi_start, nbin_x, nbin_y]
%
% x_roi are the x-coordinates in the region-of-interest space.
%
% y_roi are the y-coordinates in the region-of-interest space.
%
% RETURNS
%
% x_img are the x-coordinates in the full-image space
%
% y_img are the y-coordinates in the full-image space
%
% See Also img2roi_coord


% GlobalCoordinates(ii) = ...
%     ( transform_params(ii) * NBin(ii) - floor(NBin(ii) / 2) ) + ROIs(2 * ii - 1) - 1;

x_roi_start = transform_params(1);
y_roi_start = transform_params(2);
nbin_x = transform_params(3);
nbin_y = transform_params(4);

% do scaling transformation (un-binning)
[x_intermediate, y_intermediate, affine_xform1] = ...
    affine_transform([nbin_x 0; 0 nbin_y],[0; 0], x_roi, y_roi);
% next translate our coordinates to be referenced to the original image. This
% transformation means that (1,1) -> (x_roi_start, y_roi_start) as we
% would expect.
[x_img, y_img, affine_xform2] = ...
    affine_transform([1 0; 0 1], [(x_roi_start - 1); (y_roi_start - 1)], x_intermediate, y_intermediate);

% kind of obnoxious that matlab uses affine transformation for row vectors
% instead of column vectors.
affine_mat = transpose(affine_xform1.T * affine_xform2.T);

end