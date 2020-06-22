function binSatTest(M,BinSizeList)
%function binSatTest(M,BinSizeList)

figure('name','Binning Saturation Test');
NRows = 2;
NCols = ceil(length(BinSizeList)/2);
for ii = 1:length(BinSizeList)
    BinnedImg = binImg(M,BinSizeList(ii),BinSizeList(ii),'Normal');
    subplot(NRows,NCols,ii)
    BinnedOD = getOD(BinnedImg);
    imagesc(BinnedOD(:,:,1),[-0.15,1]);
    axis equal
    axis image
    title(sprintf('%i x %i Binning',BinSizeList(ii),BinSizeList(ii)))
end


end

