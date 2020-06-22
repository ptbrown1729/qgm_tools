function [OD] = getOD_Alpha(Images,Alpha,PixelArea,ImagingDuration,Wavelength,QE,CountsPerPhoton,Isat) 
%[OD] = getOD(Images)
%Get absorption image from array of images. First layer is atoms, second
%layer is beam, third layer is dark.
%Other common sizes return a stack of OD images.
%Change function so that 1. can accept a stack of [Ny,Nx,3,NImgs] and
%process correctly.
%For different cases, instead of rewriting the code, reformulate so that
%array size is [Ny,Nx,3,NImgs]

atoms = Images(:,:,2);
beam = Images(:,:,3);
dark = Images(:,:,4);

divided = (atoms-dark)./(beam-dark);

divided(isnan(divided))=1;
divided(divided<=0)=1;
divided(isinf(divided))=1;

IntAtoms = getIntensity(atoms,ImagingDuration,Wavelength,QE,CountsPerPhoton,PixelArea);
IntBeam = getIntensity(beam,ImagingDuration,Wavelength,QE,CountsPerPhoton,PixelArea);

OD = -Alpha*log(divided)+(IntBeam-IntAtoms)/Isat;
