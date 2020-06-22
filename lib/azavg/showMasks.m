function [DisplayMask_Img] = showMasks(MaskStack)

DisplayMask_Img = zeros(size(MaskStack,1),size(MaskStack,2));

NMasks = size(MaskStack,3);
for ii = 1:NMasks
    CurrMask = MaskStack(:,:,ii);
    CurrMask(CurrMask==1) = ii;
    DisplayMask_Img = DisplayMask_Img + CurrMask;
end
end