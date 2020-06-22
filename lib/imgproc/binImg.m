function [img_binned, img_bin_unc, img_bin_std, img_bin_sem] = binImg(img, nbin_x, nbin_y, output_size_mode, img_sem)
%   Bins an image into square bins with size Nbin.
%
% [img_binned] = binImg(img, nbin_x, nbin_y, output_size_mode)
%
% TODO: accept arbitrary dimensional arrays for img
%
% arguments:
% ------------------------
%   
%   img: An Ny x Nx x Nimgs array
%
%   nbin_x:
%
%   nbin_y:
%
%   output_size_mode: May be 'Same' or 'Normal'. By default is 'Same'. If
%   'Same', img_binned is the same size as img. Points in the same bin are
%   represented by several pixels which are all assigned the same value. 
%   If 'Normal', points in the same bin are represented by a single pixel
%   and img_binned is smaller than img.
%
%   img_sem: the standard error of the mean associated with each point int
%   he image. This is assuming that the image is already produced by some
%   sort of averaging (e.g. if the image is an average of several images).
%   If so, the standard error of the mean of the value of each pixel is a
%   measure of the uncertainty of the value of that pixel. The sem's of all
%   the pixels in a bin can be combined as unc_I = 1 / sqrt( sum( 1 / sigma_i^2 ) )
%   where the i's index all of the pixels in bin I
%
%   outputs:
% ------------------------
%
%   img_binned: image after binning
%
%   img_unc: weighted uncertainty of the image. These values are calculated
%   from img_sem according to 1 / sqrt( sum( 1 / sigma^2 ) )
%
%   img_std: the standard deviation of values in each bin. This includes
%   'Bessels correction', i.e. we divide by 1/N-1 instead of 1/N
%
%   img_sem: the standard error of the mean for values in each bin. This is
%   1/sqrt(nbin_x * nbin_y) * img_std and is provided for convenience.

if ~exist('img_sem', 'var') || isempty(img_sem)
    img_sem = zeros( size(img) );
end

ny = size(img, 1);
nx = size(img, 2);

if ~exist('output_size_mode', 'var')
    output_size_mode = 'Same';
end

if (mod(nbin_x, 1) ~= 0) || (mod(nbin_y, 1) ~= 0)
    error('nbin_x or nbin_y was not an integer.');
end

if (nbin_x == 0) || (nbin_y==0)
    error('nbin_x or nbin_y was 0.');
end

if ny < nbin_y || nx < nbin_x
    error('Ny was less than NBinY, or Nx was less than NBinX')
end

if (nbin_x == 1) && (nbin_y == 1)
    img_binned = img;
    return;
end

if ndims(img) > 3
    % for arrays of many dimensions, loop over 3 dimensional sections. No
    % efficiency optimization here, just provide this feature for
    % convenience.
    
    % get first image to determine sizes
    [img_binned_first, img_bin_unc_first, img_bin_std_first, img_bin_sem_first]...
        = binImg(img(:, :, :, 1), nbin_x, nbin_y, output_size_mode, img_sem(:, :, :, 1));
    dims = size(img);
    dims_binned = size(img_binned_first);
    
    % declare arrays to store results
    img_binned = zeros([dims_binned(1:3), dims(4:end)]);
    img_bin_unc = zeros([dims_binned(1:3), dims(4:end)]);
    img_bin_std = zeros([dims_binned(1:3), dims(4:end)]);
    img_bin_sem = zeros([dims_binned(1:3), dims(4:end)]);
    
    % store first values
    img_binned(:, :, :, 1) = img_binned_first;
    img_bin_unc(:, :, :, 1) = img_bin_unc_first;
    img_bin_std(:, :, :, 1) = img_bin_std_first;
    img_bin_sem(:, :, :, 1) = img_bin_sem_first;
    
    % get other values
    for ii = 2 : numel(img(1, 1, 1, :))
        [img_binned(:, :, :, ii), img_bin_unc(:, :, :, ii),...
         img_bin_std(:, :, :, ii), img_bin_sem(:, :, :, ii)]...
            = binImg(img(:, :, :, ii), nbin_x, nbin_y, output_size_mode, img_sem(:, :, :, ii));
    end
    
else
    % for smaller arrays, to binning
    
    if strcmp(output_size_mode, 'Same')
    % bn image, keeping the same number of pixels by duplicating.

        mean_fun = @(block_struct) mean2(block_struct.data) * ones(size(block_struct.data));
        std_fun = @(block_struct) std2(block_struct.data) * ones(size(block_struct.data));
        sum_fun = @(block_struct) sum(sum(block_struct.data)) * ones( size(block_struct.data) );

        img_binned = zeros( size(img) );
        img_bin_unc = zeros( size(img) );
        img_bin_std = zeros( size(img) );
        img_bin_sem = zeros( size(img) );

    %     for ii = 1 : size(img, 3)
    %         img_binned(:, :, ii) = blockproc(img(:, :, ii), [nbin_y nbin_x], mean_fun);
    %         img_bin_std(:, :, ii) = sqrt( blockproc( (img(:, :, ii) - img_binned(:, :, ii)).^2, [nbin_y, nbin_x], sum_fn) / (nbin_y * nbin_x - 1) );
    %         img_bin_sem(:, :, ii) = img_bin_std / sqrt(nbin_y * nbin_x);
    %         img_bin_unc(:, :, ii) = 1 ./ sqrt(blockproc( 1 ./ img_sem(:, :, ii).^2, [nbin_y, nbin_x], sum_fn));
    %     end
    %     
    elseif strcmp(output_size_mode, 'Normal')
    % bin image, reducing number of pixels.
        mean_fun = @(block_struct) mean2(block_struct.data);
        std_fun = @(block_struct) std2(block_struct.data);
        sum_fun = @(block_struct) sum(sum(block_struct.data));

        %Check that can resize correctly.
        if (mod(ny, nbin_y) ~= 0) || (mod(nx, nbin_x) ~= 0)
            warning('Image size was not divisble by bin size. Rounded image size down.');
            extra_y_pts = mod(ny, nbin_y);
            extra_x_pts = mod(nx, nbin_x);

            img = img(1 : end - extra_y_pts, 1 : end - extra_x_pts, :);
            img_sem = img_sem(1 : end - extra_y_pts, 1 : end - extra_x_pts, :);

            ny = size(img, 1);
            nx = size(img, 2);
        end

        img_binned = zeros(ny / nbin_y, nx / nbin_x, size(img, 3));
        img_bin_unc = zeros(ny / nbin_y, nx / nbin_x, size(img, 3));
        img_bin_std = zeros(ny / nbin_y, nx / nbin_x, size(img, 3));
        img_bin_sem = zeros(ny / nbin_y, nx / nbin_x, size(img, 3));

    else
        error('Incorrect option given for output_size_mode in binImg.');
    end

    % do processing
    for ii = 1 : size(img,3)
        img_binned(:, :, ii) = blockproc(img(:, :, ii), [nbin_y nbin_x], mean_fun);
        img_bin_std(:, :, ii) = blockproc( img(:, :, ii), [nbin_y, nbin_x], std_fun);
        img_bin_sem(:, :, ii) = img_bin_std(:, :, ii) / sqrt( nbin_y * nbin_x);
        img_bin_unc(:, :, ii) = 1 ./ sqrt(blockproc( 1 ./ img_sem(:, :, ii).^2, [nbin_y nbin_x], sum_fun));
    end

end

end

