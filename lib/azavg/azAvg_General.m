function [distances, distances_sdm, azavg, sdm, wunc, npts_bin, mask_stack] = ...
    azAvg_General(img, weights, distance_grid, bin_edges)
% [Distances, Distances_SDM, Img_Avg, Img_SDM, Weighted_Unc, NPts, Mask_Stack] = ...
%     azAvg_General(img, weights, distance_grid, bin_edges)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Img: is an n-D array to be averaged. nans are ignored during the averaging.
%
% DistGrid: is an n-D array the same size as Img
% which gives the 'distances' used to choose points to average.
% Points are sorted into bins based on these distances and the BinEndPts.
% So, e.g., if you have an Ny x Nx x NImgs stack of pictures and you want to
% azimuthally average each one, provide a DistGrid of size Ny x Nx.
% azAvg_General will average over the first two dimensions only...i.e. it
% will take the azimuthal average of each individually.
% On the other hand, consider the same Ny x Nx x NImgs stack. But provide
% instead a stack of identical DistGrids, which has total size Ny x Nx x NImgs.
% azAvg_General will then simultaneously average over pictures and take and
% azimuthal average.
%
% Weights: wi = 1/(sigma_i^2), where sigma_i is typically the standard
% deviation of the mean of a given point in Img. This array must have the
% same size as Img.
% If all the points are equally weighted, set Weights = []. The function
% will treat this the same as if you had entered all ones.
%
% BinEdges: gives the bin edges.
% e.g. if BinEdges = [0,2,5,10], azAvg_General will search for all points
% in Img where the corresponding entry in DistGrid is greater than or equal
% to zero, and less than 2. These points are averaged together.
%
% Return values:
%
% Distances: is an NBins x 1 array containing the average value of DistGrid 
%over the points in each bin.
%
% BinMeanDistance_SDM: is an Nbins x 1 array containing the standard deviation 
% of the mean of these points corresponding uncertainty.
%
% Img_AzAvg: is an NBins x [size of dimensions not averaged over] matrix, which
% holds the average of all points in Img where the corresponding
% entry in DistGrid falls between points in BinEdges.
%
% Img_AzAvg_SDM: is the standard deviation of the mean of these points. This
% is the appropriate uncertainty to quote in some cases.
%
% Weighted_Unc: treats the azimuthal average process as a weighted
% average, and calculates the uncertainty of each bin from the input weights.
% That is, for each bin, it returns 1/sqrt(sum_PtsInBin(Weights))
%
% NPtsInBin: is an NBins x 1 size array containing the number of points found 
% in each bin. 
% i.e. the number of points in DistanceGrid falling between the two appropriate end points.
%
% Mask_Stack: is a 'stack' of masks describing each bin. If DistGrid is an
% array of size N1 x N2 x ... Nk, then Mask_Stack has 
% size N1 x N2 x ... x Nk x NBins.
% Each k-D slice of the stack is an array of zeros and ones, if a given point
% is in or out of a given bin respectively.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check arguments make sense
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(weights)
    weights = ones(size(img));
end

if sum(isinf(weights(:))) > 0
    weights(isinf(weights)) = 0;
    warning('Encountered Infinite Weights. Set these Weights to zero.')
end

if sum(isnan(weights(:))) > 0
    weights(isnan(weights)) = 0;
    warning('Encountered NaN Weights. Set these Weights to zero.')
end

if ~isequal(size(img), size(weights))
    error('Img and Weights were not the same size.')
end

ImgSizes = size(img);
if ~isequal(ImgSizes(1:ndims(distance_grid)), size(distance_grid))
    error('Img and DistanceGrid dimensions do not agree.')
end

if sum(isnan(img(:))) > 0
    warning('NaN values in Img for azAvg_General. NaNs will be ignored.')
end

bin_edges = sort(bin_edges);
nbins = length(bin_edges) - 1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parsing the sizes of Img and DistGrid...
% Deal with case where Img has more dimensions than DistGrid...
% We will average over the first ndim(DistGrid) dimensions, and
% leave the others untouched.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%

img_full_size = size(img);
if ndims(img) > ndims(distance_grid)
    img_extra_dims_lens = img_full_size(ndims(distance_grid) + 1:end);
else
    img_extra_dims_lens = 1;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define variables for output information
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% store information about our bins...and bin error bars.
npts_bin = zeros([nbins, img_extra_dims_lens]);
distances = zeros([nbins, img_extra_dims_lens]);
distances_sdm = zeros([nbins, img_extra_dims_lens]);
mask_stack = [];
%variables to store Azimuthal average.
azavg = zeros([nbins, img_extra_dims_lens]);
sdm = zeros([nbins, img_extra_dims_lens]);
wunc = zeros([nbins, img_extra_dims_lens]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Perform Azimuthal Average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1 : nbins
    bin_mins = bin_edges(ii);
    bin_maxs = bin_edges(ii + 1);
    
    %create a mask for a single bin.
    mask = ones(size(distance_grid));
    mask(distance_grid > bin_maxs) = 0;
    mask(distance_grid <= bin_mins) = 0;
    mask( isnan(distance_grid) ) = 0;
    %store these in a stack, in case we want them later.
    mask_stack = cat(ndims(distance_grid) + 1, mask_stack, mask);
    
    %we want our mask to be the same size as our array...so repeat it over
    %any extra array dimensions...
    mask_expanded = repmat(mask, [ones(1, ndims(distance_grid)), img_extra_dims_lens]);
    masked_grid = repmat(mask .* distance_grid, [ones(1, ndims(distance_grid)), img_extra_dims_lens]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Perform weighted mean over several array dimensions. Could not find a
    %built in matlab function to sum over several dimensions at once, so we
    %have a loop.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dist_helper = masked_grid .* (~isnan(img));
    BinMeanSquareDistanceHelper = masked_grid.^2 .* (~isnan(img));
    n_azavg_pts = mask_expanded .* (~isnan(img));
    imgsum_weighted = mask_expanded .* img .* weights;
    imgsum_sqr_weighted = mask_expanded .* img .^2 .* weights;
    sum_weights = mask_expanded .* weights .* (~isnan(img));
    for jj = 1:ndims(distance_grid)
        %TODO: This section has problems with singleton dimensions.
        %Want to remove first dimension each time, but maybe not with squeeze,
        %because this can strip off later singleton dimensions we want
        %to keep! Also sum does not operate on singleton dimensions by
        %default.
        dist_helper = squeeze(sum(dist_helper));
        
        BinMeanSquareDistanceHelper = squeeze(sum(BinMeanSquareDistanceHelper));
        
        n_azavg_pts = squeeze(sum(n_azavg_pts));
        
        imgsum_weighted = squeeze(nansum(imgsum_weighted));
        
        imgsum_sqr_weighted = squeeze(nansum(imgsum_sqr_weighted));
        sum_weights = squeeze(sum(sum_weights));
        %TODO weights sum can be wrong if we have nans around.
    end
    
     %TODO with NaN handling...the number of points may differ for each
    %slice...NPts should have size [NBins x NExtraDims]
    
    %get mean bin distance and its uncertainty...
    npts_bin(ii,:) = n_azavg_pts(:);
    distances(ii,:) = dist_helper(:) ./ n_azavg_pts(:);
    BinMeanSquareDistance = BinMeanSquareDistanceHelper./n_azavg_pts;
    %needed to correct for factor of (n-1) in SD...
    %also handle the case where there is only a single point in a bin.
    NMinusOnePts = n_azavg_pts - 1;
    NMinusOnePts(NMinusOnePts==0) = 1;
    distances_sdm(ii,:) = sqrt(n_azavg_pts(:)./NMinusOnePts(:)) ...
                       .* sqrt(BinMeanSquareDistance(:) ...
                        - permute(distances(ii, :), [2, 1]) .^2) ./ sqrt(n_azavg_pts(:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Assign the average values for this bin to function output parameters.
    %In the weighted case, the SDM of the mean calculated here may not make
    %sense...it is really intended for the case where 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    azavg(ii,:) = imgsum_weighted(:) ./ sum_weights(:);
    Img_Squares_AzAvg = imgsum_sqr_weighted ./ sum_weights;
    %permute because Img_AzAvg(ii,:) has size 1 x N, but Img_AzAvg(:) has size N x 1...
    sdm(ii,:) = sqrt(n_azavg_pts(:)./NMinusOnePts(:)) ...
                 .* sqrt(Img_Squares_AzAvg(:) - permute(azavg(ii,:),[2,1]).^2) ...
                 ./ sqrt(n_azavg_pts(:));
    wunc(ii,:) = 1 ./ sqrt(sum_weights(:));
end