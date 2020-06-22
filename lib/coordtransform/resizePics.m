function [ResPics, ROIStart_X, ROIStart_Y] = resizePics(Pics, FinalSize)
%resizePics Resize pictures, keeping center fixed as best as
%possible.
%
%   [ResPics,ROIStart_X,ROIStart_Y] = resizePics(Pics,FinalSize)
%   Given a Ny x Nx x NImgs array, Pics, and a vector FinalSize
%   specifying the first two dimensions of the final array,
%   return a cropped version of Pics which is approximately
%   centered. Also returns ROIStart_X and ROIStart_Y, which
%   give the coordinate of the upper left most point in the
%   region of interest in terms of the initial picture
%   coordinates. The region of interest is the same for all
%   slices in Pics.
%
%   TODO: add option on how to resize.
ResPics = zeros(size(Pics));
if size(Pics,1) > FinalSize(1)
    SizeToCut = size(Pics,1) - FinalSize(1);
    SizeToCutTop = ceil(SizeToCut / 2);
    SizeToCutBottom = SizeToCut - SizeToCutTop + 1;
    ResPics = Pics(SizeToCutTop:end - SizeToCutBottom, :, :);
    ROIStart_Y = SizeToCutTop;
elseif size(Pics,1) < FinalSize(1)
    SizeToPad = abs(size(Pics, 1) - FinalSize(1));
    SizeToPadTop = ceil(SizeToPad / 2);
    SizeToPadBottom = SizeToPad - SizeToPadTop + 1;
    ResPics = padarray(ResPics, [SizeToPadTop, 0, 0], 0, 'pre');
    ResPics = padarray(ResPics, [SizeToPadBottom, 0, 0], 0, 'post');
    ROIStart_Y = -SizeToPadTop + 1;
else
    ResPics = Pics;
    ROIStart_Y = 1;
end

if size(Pics,2) > FinalSize(2)
    SizeToCut = size(Pics,2) - FinalSize(2);
    SizeToCutLeft = ceil(SizeToCut/2);
    SizeToCutRight = SizeToCut-SizeToCutLeft+1;
    ResPics = ResPics(:, SizeToCutLeft:end - SizeToCutRight, :);
    ROIStart_X = SizeToCutLeft;
elseif size(Pics,2) < FinalSize(2)
    SizeToPad = abs(size(Pics, 2) - FinalSize(2));
    SizeToPadLeft = ceil(SizeToPad / 2);
    SizeToPadRight = SizeToPad - SizeToPadLeft + 1;
    ResPics = padarray(ResPics, [0, SizeToPadLeft, 0], 'pre');
    ResPics = padarray(ResPics, [0, SizeToPadRight, 0], 'post');
    ROIStart_X = -SizeToPadLeft + 1;
else
    ResPics = ResPics;
    ROIStart_X = 1;
end
end