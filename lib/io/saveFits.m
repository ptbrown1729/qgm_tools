function [] = saveFits(ImageStack, SavePath)
%[] = saveFits(ImageStack,SavePath)
%save matrix as a .fits image file
%TODO: add ability to save metadata

Images = flip(ImageStack, 1); %so orientation matches andor solis program.
fitswrite(Images, SavePath);

end