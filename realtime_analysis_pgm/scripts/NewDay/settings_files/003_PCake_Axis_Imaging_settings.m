CurrentFolder = pwd;
Experiment = 'gaussianExptUpdate';
%Hardware Properties.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Camera = 'Zyla4.2';
PictureType = 'Absorption';
ActivePictures = 2:4;
PixelSize = 6.5e-6;
QuantumEfficiency = 0.75;
CountsPerPhoton = 1/0.48;
ImagingDuration = 10e-6;
ImagingAxis = 'Pancake';
Magnification = 750/150;
HardwareBinSizeH = 1;
HardwareBinSizeV = 1;

%Software Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Binning
SoftwareBinSizeH = 8;
SoftwareBinSizeV = 8;

%Normalize OD
NormalizeOD = 1;
NormRegionLimits = [500,700,800,900];

%ROI
WindowLimits = [600,1400,800,1500];

%Remove Stripes
RemoveStripes = 0;
StripesOrientation = 'Horizontal';
StripesRegionLimits = [1,50,130,450];

%Full range for Neo: [1,2560,1,2160]
%Fitting Properties.
%[Gauss,TF] %[CxG,CyG,SxG,SyG,AmpG,ThetaG,Bg]
Fitting = 1;
FitInitialParameters = [1050,1150,100,10,0.5,0,0];
FitFixedParameters = [0,0,0,0,0,0,0];
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
Lims = [-0.15,0.4];
ImageLimits = [Lims];
FitLimits = [Lims];
ResidualLimits = [-0.25,0.2];

%Markers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MarkerLocations = [1050,1165]; %[310,213]; %[311,218]; % [185,198]; %[178,195]; % [204.6,192];
%Debug Mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DebugMode = 0;

