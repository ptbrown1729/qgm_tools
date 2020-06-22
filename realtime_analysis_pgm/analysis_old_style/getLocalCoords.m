function LocalCoords = getLocalCoords(GlobalCoords, ROIs, NBin)
%
%   TODO: replace this function with get_roi_coords
%
%LocalCoords = getLocalCoords(GlobalCoords,ROIs,NBin)
%%%arguments%%%
%GlobalCoordinates = [X1_global, X2_global, ..., Xn_global]is the list of
%coordinates to transform.
%ROI = [[X1_start, X1_end], [X2_start, X2_end], ..., [Xn_start, Xn_end]] is a list
%of regions of interest that correspond to each coordinate. If the same ROI
%is to be used for all coordinates, you can specify only a single ROI =
%[XStart, XEnd]
%NBin = [NBin1, NBin2, ..., NBin_n] is a list of binning values to be used for
%each coordinate. This is the number of global coordinates per local
%coordinate. If the same binning is to be used for all coordinates, you only
%need to specify that single number.
%%%Return values%%%
%LocalCoords = [X1_local, X2_local, ...., Xn_local]

%%%Argument checking
if length(ROIs)==2 && length(GlobalCoords)>1
    %if we are always using the same ROI, it is more convenient to enter it
    %only once, so handle that case.
    ROIs = repmat(ROIs, [1, length(GlobalCoords)]);
end

if 2*length(GlobalCoords) ~= length(ROIs)
    errorStruct.message = 'ROIs was not twice as long as Coordinates.';
    errorStruct.identifier = 'getGlobalCoordsFromROI:samelen';
    error(errorStruct)
end

if sum(isnan(GlobalCoords)) ~= 0 || sum(isinf(GlobalCoords)) ~= 0
    errorStruct.message = 'Coordinates contained nans or infs.';
    errorStruct.identifier = 'getGlobalCoordsFromROI:nanOrinf';
    error(errorStruct);
end

if sum(isnan(ROIs)) ~= 0 || sum(isinf(ROIs)) ~= 0
    errorStruct.message = 'ROIs contained nans or infs.';
    errorStruct.identifier = 'getGlobalCoordsFromROI:nanOrinf';
    error(errorStruct);
end

if ~exist('NBin', 'var')
    NBin = ones(length(GlobalCoords),1);
end

if (length(NBin) == 1 && length(GlobalCoords) > 1)
    %so can input single NBin value to be applied to all coordinates.
    NBin = NBin*ones(length(GlobalCoords),1);
end

if sum(isnan(NBin)) ~= 0 || sum(isinf(NBin)) ~= 0
    errorStruct.message = 'NBin contained nans or infs.';
    errorStruct.identifier = 'getGlobalCoordsFromROI:nanOrinf';
    error(errorStruct);
end

if length(GlobalCoords) ~= length(NBin)
    errorStruct.message = 'NBin was different size than Coordinates.';
    errorStruct.identifier = 'getGlobalCoordsFromROI:samelen';
    error(errorStruct)
end

%%%Convert coordinates
LocalCoords = zeros(1, length(GlobalCoords));
for ii = 1:length(GlobalCoords)
    LocalCoords(ii) = ceil( (GlobalCoords(ii) - (ROIs(2 * ii - 1) - 1)) / NBin(ii) );
    %Xl = ceil[(Xg - X1 +1)/NBin]
end

end