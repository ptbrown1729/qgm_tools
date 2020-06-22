CurrentFolder = pwd;
Experiment = 'gaussianExpt';
%Hardware Properties.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Camera = 'Guppy F-080';
Manufacturer = 'AVT';
PixelSize = 4.65e-6;
QuantumEfficiency = 0.225; %guppy
CountsPerPhoton = 1;
ImagingDuration = 10e-6;
ImagingAxis = 'Vertical';
Magnification = 30/250;
HardwareBinSizeH = 1;
HardwareBinSizeV = 1;

%Software Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Binning
SoftwareBinSizeH = 2;
SoftwareBinSizeV = 2;

%Normalize OD
NormalizeOD = 0;
NormRegionLimits = [100,400,300,350];

%ROI
WindowLimits = [500,800,200,500];

%Remove Stripes
RemoveStripes = 0;
StripesOrientation = 'Horizontal';
StripesRegionLimits = [1,50,130,450];

%Full range for Neo: [1,2560,1,2160]
%Fitting Properties.
%[Gauss,TF] %[CxG,CyG,SxG,SyG,AmpG,ThetaG,Bg]
Fitting = 1;
FitInitialParameters = [640,332,4,21,0.5,0,0];
FitFixedParameters = [0,0,0,0,0,0,0];
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
MarkerLocations = [314,215]; %[311,218]; % [185,198]; %[178,195]; % [204.6,192];
%Debug Mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DebugMode = 0;

