function ROI = getROI(Img, XStart, XSize, YStart, YSize)
%getROI return a region of interest from a matrix.
%
%   ROI = getROI(obj,Img,XStart,XSize,YStart,YSize) Returns a
%   subregion of the matrix Img as ROI =
%   Img(YStart:YStart+YSize-1,XStart:XStart+XSize-1). If the
%   region of interest extends outside the initial matrix,
%   these regions are treated as zeros.
%
%   TODO: make work for stacks of images.
%   TODO: probably should treat regions outside image as NANS!!!

ROI = zeros([YSize, XSize]);
ShiftImg = getShiftedMat(Img, -XStart + 1, -YStart + 1, 0);

if XSize <= size(ShiftImg, 2) && YSize <= size(ShiftImg, 1)
    ROI = ShiftImg(1:YSize, 1:XSize);
elseif XSize > size(ShiftImg, 2)
    ROI(:,1:size(ShiftImg, 2)) = ShiftImg(1:YSize, :);
elseif YSize > size(ShiftImg, 1)
    ROI(1:size(ShiftImg, 1), :) = ShiftImg(:, 1:XSize);
else
    ROI(1:size(ShiftImg, 1), 1:size(ShiftImg, 2)) = ShiftImg;
end

end