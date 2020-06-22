function [x_roi, y_roi, affine_mat] = img2roi_coord(transform_params, x_full, y_full, quantize)
% Get coordinates in the region-of-interest space.
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
% x_full are x-coordinates in the 'full' space.
%
% y_full are y-coordinates in the 'full' space.
%
% quantize: If true, this function will return integer coordinates which
% represent which binned pixel the initial pixel is a part of. If quantize
% is false, this function returns coordinates using the conventions
% described above.
%
% RETURNS
%
% x_roi are x-coordinates in the region-of-interest space
%
% y_roi are y-coordinates in the region-of-interest space

if ~exist('quantize', 'var') || isempty(quantize)
    quantize = 0;
end

x_roi_start = transform_params(1);
y_roi_start = transform_params(2);
nbin_x = transform_params(3);
nbin_y = transform_params(4);

% should I use the ceil function here? It seems like it depends on if I'm
% only allowing values at exactly integer values. This is like only
% allowing values at each pixel of an image...
% TODO: Maybe I should write a more general affine transformation function
% and this function would be a wrapper/possibly with additional ceil
% function

% This method is simpler to write down, and was initial approach. Thought
% maybe affine transform approach would unify my many disparate coordinate
% transformations.
% x_roi = ceil( (x_full - (x_roi_start - 1)) / nbin_x );
% y_roi = ceil( (y_full - (y_roi_start - 1)) / nbin_y );

% first translate our coordinates to be referenced to the roi. This
% transformation means that (x_roi_start, y_roi_start) -> (1, 1), as we
% would expect.
[x_intermediate, y_intermediate, affine_xform1] = affine_transform([1 0; 0 1], [-(x_roi_start - 1); -(y_roi_start - 1)], x_full, y_full);
% do scaling transformation (binning)
[x_roi, y_roi, affine_xform2] = affine_transform([1/nbin_x 0; 0 1/nbin_y], [0; 0], x_intermediate, y_intermediate);

% kind of obnoxious that matlab uses affine transformation for row vectors
% instead of column vectors.
affine_mat = transpose(affine_xform1.T * affine_xform2.T);

if quantize
    x_roi = ceil(x_roi);
    y_roi = ceil(y_roi);
end

end