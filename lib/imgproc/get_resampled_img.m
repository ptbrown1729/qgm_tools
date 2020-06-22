function [resampled_img, dx, dy, std, sem, unc] = get_resampled_img(img, num_pix, x_roi, y_roi, mode, img_unc)
%   Resample a region of interest of an image at a new number of pixels. If
%   the number of pixels is smaller than the initial number this
%   corresponds to 'binning' the image to decrease the noise. If the number
%   of pixels is larger than the initial image this is oversampling.
%
%   arguments:
%   ------------
%   img: A 2d array to be resampled
%
%   num_pix: [num_pix_x, num_pix_y]. Number of pixels in the resampled
%   array.
%
%   x_roi: [x_start, x_end], corresponding to the edges of the resampled
%   image. We regard the initial image as spanning the coordinates [0.5,
%   size(img, 2) + 0.5], i.e. the centers of the pixels at positions 1, ...
%   size(img, 2). If no x_roi is supplied, the default value is [0.5,
%   size(img, 2) + 0.5].
%
%   y_roi:
%
%   mode: may be either 'mean' or 'sum'. If sum, then the overall weight of
%   the initial image is preserved. If 'mean', then the density of the
%   image is preserved.
%
%   return:
%   ------------
%   resampled_img: an num_pix(2) x num_pix(1) image
%
% TODO: add weighted uncertainty propogation to this ... have to think
% about how to do that appropriately.

% ROI arguments
if ~exist('img_unc', 'var') || isempty(img_unc)
    img_unc = ones( size(img) );
end

if ~exist('x_roi', 'var') || isempty(x_roi)
    x_roi = [0.5, size(img, 2) + 0.5];
end

if ~exist('y_roi', 'var') || isempty(y_roi)
    y_roi = [0.5, size(img, 1) + 0.5];
end

if x_roi(1) >= x_roi(2) || y_roi(1) >= y_roi(2)
    error('roi coordinates must be increasing');
end

if x_roi(1) < 0.5 || x_roi(2) > size(img, 2) + 0.5 ...
   || y_roi(1) < 0.5 || y_roi(2) > size(img, 1) + 0.5
    error('roi goes outside of image. Image runs from [x1, x2] x [y1, y2] = [0.5, lx + 0.5], [0.5, ly + 0.5].');
end

% mode argument
if ~exist('mode', 'var') || isempty(mode)
    mode = 'mean';
end

allowed_modes = {'mean', 'sum'};
fn = @(c) strcmp(c, mode);
if ~any(cellfun(fn, allowed_modes))
    error('mode direction argument was not one of the allowed keywords.');
end

nimgs = size(img, 3);
if nimgs > 1
    % if we have multiple images, call our function on each individual
    % image
    resampled_img = zeros([num_pix(1), num_pix(2), nimgs]);
    std = zeros([num_pix(1), num_pix(2), nimgs]);
    sem = zeros([num_pix(1), num_pix(2), nimgs]);
    unc = zeros([num_pix(1), num_pix(2), nimgs]);
    
    for ii = 1 : nimgs
        [resampled_img(:, :, ii), dx, dy, std(:, :, ii), sem(:, :, ii), unc(:, :, ii)] ...
         = get_resampled_img(img(:, :, ii), num_pix, x_roi, y_roi, mode);
    end
    
else
    % coordinates of initial image
    [xx_i, yy_i] = meshgrid( 1 : size(img, 2), 1:size(img, 1));

    % coordinates of new image
    dx = (x_roi(2) - x_roi(1)) / num_pix(1);
    dy = (y_roi(2) - y_roi(1)) / num_pix(2);
    [xx_f, yy_f] = meshgrid( x_roi(1) + dx/2 : dx : x_roi(2) - dx/2,...
                             y_roi(1) + dy/2 : dy : y_roi(2) - dy/2); 
                         
    resampled_img = zeros( size(xx_f) );
    std = zeros( size(xx_f) );
    unc = zeros( size(xx_f) );
    for ii = 1 : numel(resampled_img)
        % weighted mean
        overlaps = get_pixel_overlap(xx_i, 1, xx_f(ii), dx) .* ...
                  get_pixel_overlap(yy_i, 1, yy_f(ii), dy);
        % note that sum(weights) = dx * dy
        resampled_img(ii) = sum(sum( img .*  overlaps )) / (dx * dy);
                  
        % calculate weighted uncertainty
        % This has several convenient properties: if all weights are equal
        % then it reduces to uncertainty propogation, i.e. unc = sqrt( sum( unc.^2) );
        % On the other hand, if all sigmas are equal, it reduces to the
        % standard deviation of the mean, unc = sigma / sqrt( \sum w_i),
        % where the total weight is the area being summed over (i.e. number
        % of pixels, which may be fractional in this case).
        
        % ignore points with zero uncertainty for calculating uncertainty
        unc_weights = 1 ./ img_unc.^2;
        unc_weights( isinf(unc_weights) ) = 0;
        
        % this form corrects issue where uncertainty would explode if we
        % had a resampled point that had very small overlap (e.g. 1%) with
        % some pixel. i.e. if unc = 1 and w_i = 0.01, then without the
        % second term in the product below we would get 
        % unc = 1 / sqrt( 0.01) = 10. Whereas we should expect the
        % uncertainty to decrease. 
        unc(ii) = 1 ./ sqrt( sum(sum( overlaps .* unc_weights  )) ) .* sqrt( sum(overlaps(unc_weights ~= 0)) / dx /dy); 
        % handle case where none of the points in the ROI have finite
        % uncertainty.
        unc( isinf(unc) ) = 0;
        unc( isnan(unc) ) = 0;
        
        % weighted sample variance
        std(ii) = sqrt( sum(sum( (img - resampled_img(ii)).^2 .* overlaps )) / (dx * dy - 1) );
    end
    % this should at least agree with the sem in the case where we are
    % binning pixels. Maybe should think more about the weighted
    % average case if thsi still makes sense.
    sem = std / sqrt(dx * dy);
    
    if strcmp(mode, 'sum')
        resampled_img = resampled_img * ( dx * dy);
    end
end
    
end