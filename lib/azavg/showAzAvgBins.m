function [ImgWithBinEdgesMasked] = showAzAvgBins(Image,DistanceGrid,BinEndPts)
%[ImgWithBinEdgesMasked] = showAzAvgBins(Image,DistanceGrid,BinEndPts)
%Returns an input image where points on the edges of bins used for
%azimuthal averaging are given a large value. This image can be displayed
%as a visual aid to monitor azimuthal averaging.

%Check arguments...
if sum(isnan(BinEndPts)) ~= 0 || sum(isinf(BinEndPts)) ~= 0
    errorStruct.message = 'BinEndPts contained nans or infs.';
    errorStruct.identifier = 'showAzAvgBins:nanOrinf';
    error(errorStruct);
end

    BinEdgeMask = zeros(size(DistanceGrid));
    for ii = 1:length(BinEndPts)
        PtA = BinEndPts(ii)-1;
        PtB = BinEndPts(ii)+1;
        BinEdgeMask((PtA<DistanceGrid)&(DistanceGrid<PtB)) = 1;
    end
    
    ImgWithBinEdgesMasked = Image;
    ImgWithBinEdgesMasked(BinEdgeMask==1)=10;

end