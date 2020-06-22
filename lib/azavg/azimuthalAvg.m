function [Distances,Distances_SDM,Img_Avg,Img_SDM,Weighted_Unc,NPts,Mask_Stack] = azimuthalAvg(Img,Weights,DistGrid,BinEdges)
%[Distances,Distances_SDM,Img_AzAvg,Img_AzAvg_SDM,Weighted_Unc,NPts,Mask_Stack] = azAvg_General(Img,Weights,DistGrid,BinEdges)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Img is an n-D array to be averaged.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DistGrid is an n-D array the same size as Img
%which gives the 'distances' used to choose points to average.
%Points are sorted into bins based on these distances and the BinEndPts.
%%%%So, e.g., if you have an Ny x Nx x NImgs stack of pictures and you want to
%azimuthally average each one, provide a DistGrid of size Ny x Nx.
%azAvg_General will average over the first two dimensions only...i.e. it
%will take the azimuthal average of each individually.
%%%%On the other hand, consider the same Ny x Nx x NImgs stack. But provide
%instead a stack of identical DistGrids, which has total size Ny x Nx x NImgs.
%azAvg_General will then simultaneously average over pictures and take and
%azimuthal average.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Weights, wi = 1/(sigma_i^2), where sigma_i is typically the standard
%deviation of the mean of a given point in Img. This array must have the
%same size as Img.
%If all the points are equally weighted, set Weights = []. The function
%will treat this the same as if you had entered all ones.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BinEdges gives the bind edges.
%%%%e.g. if BinEdges = [0,2,5,10], azAvg_General will search for all points
%in Img where the corresponding entry in DistGrid is greater than or equal
%to zero, and less than 2. These points are averaged together.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Return Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Distances is an NBins x 1 array containing the average value of DistGrid 
%over the points in each bin.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BinMeanDistance_SDM is an Nbins x 1 array containing the standard deviation 
%of the mean of these points corresponding uncertainty.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Img_AzAvg is an NBins x [size of dimensions not averaged over] matrix, which
%holds the average of all points in Img where the corresponding
%entry in DistGrid falls between points in BinEdges.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Img_AzAvg_SDM is the standard deviation of the mean of these points. This
%is the appropriate uncertainty to quote in some cases.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Weighted_Unc treats the azimuthal average process as a weighted
%average, and calculates the uncertainty of each bin from the input weights.
%That is, for each bin, it returns 1/sqrt(sum_PtsInBin(Weights))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NPtsInBin is an NBins x 1 size array containing the number of points found 
%in each bin. 
%i.e. the number of points in DistanceGrid falling between the two appropriate end points.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mask_Stack is a 'stack' of masks describing each bin. If DistGrid is an
%array of size N1 x N2 x ... Nk, then Mask_Stack has 
%size N1 x N2 x ... x Nk x NBins.
%Each k-D slice of the stack is an array of zeros and ones, if a given point
%is in or out of a given bin respectively.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check arguments make sense
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(Weights)
    Weights = ones(size(Img));
end

if sum(isinf(Weights(:)))>0
    Weights(isinf(Weights)) = 0;
    warning('Encountered Infinite Weights. Set these Weights to zero.')
end

if sum(isnan(Weights(:)))>0
    Weights(isnan(Weights)) = 0;
    warning('Encountered NaN Weights. Set these Weights to zero.')
end

if ~isequal(size(Img),size(Weights))
    error('Img and Weights were not the same size')
end

ImgSizes = size(Img);
if ~isequal(ImgSizes(1:ndims(DistGrid)),size(DistGrid))
    error('Img and DistanceGrid dimensions do not agree.')
end

BinEdges = sort(BinEdges);
NBins = length(BinEdges)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parsing the sizes of Img and DistGrid...
%Deal with case where Img has more dimensions than DistGrid...
%We will average over the first ndim(DistGrid) dimensions, and
%leave the others untouched.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ImgFullSize = size(Img);
if ndims(Img)>ndims(DistGrid)
    ImgExtraDimensionLengths = ImgFullSize(ndims(DistGrid)+1:end);
else
    ImgExtraDimensionLengths = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define variables for output information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store information about our bins...and bin error bars.
NPts = zeros(NBins,1);
Distances = zeros(NBins,1);
Distances_SDM = zeros(NBins,1);
%variables to store Azimuthal average.
Img_Avg = zeros([NBins,ImgExtraDimensionLengths]);
Img_SDM = zeros([NBins,ImgExtraDimensionLengths]);
Weighted_Unc = zeros([NBins,ImgExtraDimensionLengths]);
%masks
Mask_Stack = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Perform Azimuthal Average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:NBins
    BinMin = BinEdges(ii);
    BinMax = BinEdges(ii+1);
    
    %create a mask for a single bin.
    Mask = ones(size(DistGrid));
    Mask(DistGrid>BinMax) = 0;
    Mask(DistGrid<=BinMin) = 0;
    %store these in a stack, in case we want them later.
    %Mask_Stack(:,:,ii) = Mask;
    Mask_Stack = cat(ndims(DistGrid)+1,Mask_Stack,Mask);
    
    NAzAvgPts = sum(Mask(:));
    NPts(ii) = NAzAvgPts;
    %get mean bin distance and its uncertainty...
    MaskedGrid = Mask.*DistGrid;
    Distances(ii) = sum(MaskedGrid(:))/NAzAvgPts;
    BinMeanSquareDistance = sum(MaskedGrid(:).^2)/NAzAvgPts;
    Distances_SDM(ii) = sqrt(BinMeanSquareDistance - Distances(ii)^2)/sqrt(NAzAvgPts);
    
    %we want our mask to be the same size as our array...so repeat it over
    %any extra array dimensions...
    Mask_Expanded_To_All_Dimensions = repmat(Mask,[ones(1,ndims(DistGrid)),ImgExtraDimensionLengths]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Perform weighted mean over several array dimensions. Could not find a
    %built in matlab function to sum over several dimensions at once, so we
    %have a loop.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Weighted_ImgSum = Mask_Expanded_To_All_Dimensions.*Img.*Weights;
    Weighted_ImgSumSquares = Mask_Expanded_To_All_Dimensions.*Img.^2.*Weights;
    Weights_Sum = Mask_Expanded_To_All_Dimensions.*Weights;
    for jj = 1:ndims(DistGrid)
        Weighted_ImgSum = squeeze(sum(Weighted_ImgSum));
        Weighted_ImgSumSquares = squeeze(sum(Weighted_ImgSumSquares));
        Weights_Sum = squeeze(sum(Weights_Sum));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Assign the average values for this bin to function output parameters.
    %In the weighted case, the SDM of the mean calculated here may not make
    %sense...it is really intended for the case where 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Img_Avg(ii,:) = Weighted_ImgSum(:)./Weights_Sum(:);
    Img_Squares_AzAvg = Weighted_ImgSumSquares./Weights_Sum;
    %permute because Img_AzAvg(ii,:) has size 1 x N, but Img_AzAvg(:) has
    %size N x 1...
    Img_SDM(ii,:) = sqrt(Img_Squares_AzAvg(:)-permute(Img_Avg(ii,:),[2,1]).^2)/sqrt(NAzAvgPts);
    Weighted_Unc(ii,:) = 1./sqrt(Weights_Sum(:));

end
