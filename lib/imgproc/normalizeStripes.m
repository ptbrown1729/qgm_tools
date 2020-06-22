function [ImgNorm] = normalizeStripes(Img,NormRegion,Orientation)
%function [ImgNorm] = normalizeStripes(Img,NormRegion,Orientation)
%NORMALIZESTRIPES removes stripes due to image readout noise.
%Norm Region should have the same number of rows as Img.
%Can do this for a stack of images.
%Must be very careful about the normalization region you choose. This
%region should be small enough that the light level is relatively constant.
%Otherwise this function will change the global light level average.

NFrames = size(Img,3);
ImgNorm = zeros(size(Img));

if ischar(Orientation)
    if strcmp(Orientation,'Horizontal')
        Orientation = 2;
    elseif strcmp(Orientation,'Vertical')
        Orientation = 1;
    end
end

if Orientation == 2
    for ii = 1:NFrames
        Frame = Img(:,:,ii);
        FrameNormRegion = NormRegion(:,:,ii);
        Avg = mean(FrameNormRegion(:));
        RowAvg = mean(FrameNormRegion,2);
        [~,SizeX] = size(Frame);
        Correction = kron(ones(1,SizeX),RowAvg-Avg);
        ImgNorm(:,:,ii) = Frame-Correction;
    end
elseif Orientation == 1
    for ii = 1:NFrames
        Frame = Img(:,:,ii);
        FrameNormRegion = NormRegion(:,:,ii);
        Avg = mean(FrameNormRegion(:));
        ColAvg = mean(FrameNormRegion,1);
        [SizeY,~] = size(Frame);
        Correction = kron(ones(SizeY,1),ColAvg-Avg);
        ImgNorm(:,:,ii) = Frame-Correction;
    end
end




end

