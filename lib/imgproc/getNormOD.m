function [NormImages] = getNormOD(Images,NormRegions)
%[NormImages] = getNormOD(Images,NormRegions) 
%takes atoms, beam, and dark image and rescales beam image
%to give zero background.

Atoms = Images(:,:,1);
Beam = Images(:,:,2);
Dark = Images(:,:,3);

AtomsNormRegion = NormRegions(:,:,1);
BeamNormRegion = NormRegions(:,:,2);
DarkNormRegion = NormRegions(:,:,3);

%Dark region should not be part of the scaling. This correction becomes
%imiportant for low counts in atom/beam pictures.
AtomAvg = mean(mean(AtomsNormRegion-DarkNormRegion));
BeamAvg = mean(mean(BeamNormRegion-DarkNormRegion));
BeamNorm = (Beam-Dark)*AtomAvg/BeamAvg+Dark;

NormImages = cat(3,Atoms,BeamNorm,Dark);








end

