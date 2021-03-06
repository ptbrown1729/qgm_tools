CurrentFolder = pwd;
Experiment = 'gaussianExptUpdate'; %'gaussianExptUpdate';
%Hardware Properties.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Camera = 'Zyla4.2';
PictureType = 'Absorption';
ActivePictures = 2:4;
PixelSize = 6.5e-6;
QuantumEfficiency = 0.75;
CountsPerPhoton = 1/0.48;
ImagingDuration = 10e-6;
ImagingAxis = 'Vertical';
Magnification = 750/25;
HardwareBinSizeH = 4;
HardwareBinSizeV = 4;

%Software Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Binning
SoftwareBinSizeH = 4;
SoftwareBinSizeV = 4;

%Normalize OD
NormalizeOD = 0;
NormRegionLimits = [100,400,300,350];

%ROI
WindowLimits = [1,512,1,512];

%Remove Stripes
RemoveStripes = 0;
StripesOrientation = 'Horizontal';
StripesRegionLimits = [1,50,130,450];

%Full range for Neo: [1,2560,1,2160]
%Fitting Properties.
%[Gauss,TF] %[CxG,CyG,SxG,SyG,AmpG,ThetaG,Bg]
Fitting = 1;
FitInitialParameters = [354, 250, 24, 21 ,0.5, 0, 0];
FitFixedParameters = [0, 0, 0, 0, 0, 1, 0];
FitLowerLims = [-inf, -inf, 10, 10, -1, -inf, -0.5];
FitUpperLims = [inf, inf, 250, 250, 2, inf, 0.5];
UseFitAutoGuess = 0;

%Display Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Image Display
ColorMap = 'bone';
ImageMultiplicity = [1];
FitMultiplicity = [1];
ResidualBool = [1];
%Locations
ImageLocations =[1,1];
FitLocations = [2,1];
ResidualLocations = [1,2];
%Limits
Lims = [-0.15,1];
ImageLimits = [Lims];
FitLimits = [Lims];
ResidualLimits = [-0.25,0.2];

%Markers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MarkerLocations = [336,236]; %[311,218]; % [185,198]; %[178,195]; % [204.6,192];
%Debug Mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DebugMode = 0;

