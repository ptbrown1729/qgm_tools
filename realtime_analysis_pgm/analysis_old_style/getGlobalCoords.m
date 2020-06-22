function GlobalCoordinates = getGlobalCoords(Coordinates, ROIs, NBin)
%GlobalCoordinate = getGlobalCoords(Coordinates,ROIs,NBin)
%%%Arguments%%%
%Coordinates = [X1_local, X2_local, ....Xn_local] is the list of
%coordinates to transform.
%ROI = [[X1_start,X1_end], [X2_start,X2_end], ..., [Xn_start,Xn_end]] is a list
%of regions of interest that correspond to each coordinate. If the same ROI
%is to be used for all coordinates, you can specify only a single ROI =
%[XStart, XEnd]
%NBin = [NBin1, NBin2, ..., NBin_n] is a list of binning values to be used for
%each coordinate. This is the number of global coordinates per local
%coordinate. If the same binning is to be used for all coordinates, you only
%need to specify that single number.
%%%Return values%%%
%GlobalCoordinates = [X1_global, X2_global, ..., Xn_global]

%%%Argument checking
if length(ROIs) == 2 && length(Coordinates) > 1
    %if we are always using the same ROI, it is more convenient to enter it
    %only once, so handle that case.
    ROIs = repmat(ROIs, [1, length(Coordinates)]);
end

if 2*length(Coordinates) ~= length(ROIs)
    errorStruct.message = 'ROIs was not twice as long as Coordinates.';
    errorStruct.identifier = 'getGlobalCoordsFromROI:samelen';
    error(errorStruct)
end

if sum(isnan(Coordinates)) ~= 0 || sum(isinf(Coordinates)) ~= 0
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
    %handle case if NBin is not supplied as argument.
    NBin = ones(length(Coordinates),1);
end

if (length(NBin)==1 && length(Coordinates)>1)
    %so can input single NBin value to be applied to all coordinates.
    NBin = NBin*ones(length(Coordinates),1);
end

if sum(isnan(NBin))~=0 || sum(isinf(NBin))~=0
    errorStruct.message = 'NBin contained nans or infs.';
    errorStruct.identifier = 'getGlobalCoordsFromROI:nanOrinf';
    error(errorStruct);
end

if length(Coordinates)~=length(NBin)
    errorStruct.message = 'NBin was different size than Coordinates.';
    errorStruct.identifier = 'getGlobalCoordsFromROI:samelen';
    error(errorStruct)
end

%Convert coordinates
GlobalCoordinates = zeros(1, length(Coordinates));
for ii = 1:length(Coordinates)
    %first paranthesis 'undo' binning as much as possible.
    %Then ROI coords added.
    GlobalCoordinates(ii) = ...
        ( Coordinates(ii) * NBin(ii) - floor(NBin(ii) / 2) ) + ROIs(2 * ii - 1) - 1;
    %Xg = (Xl*NBin - NBin/2) + X1 - 1
    %Subtract NBin/2 to get approximately the right pixel. Since Xg
    %oversamples Xl, can't choose an exact pixel.
    %choose floor instead of ceiling becuase this correctly reproduces case
    %where Xg = 1, Xl = 1, NBin = 1.
end

end