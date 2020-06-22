CurrentFolder = pwd;
Experiment = 'gaussianExptUpdate';
%Hardware Properties.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Camera = 'Zyla4.2';
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
WindowLimits = [200,470,100,350];

%Remove Stripes
RemoveStripes = 0;
StripesOrientation = 'Horizontal';
StripesRegionLimits = [1,50,130,450];

%Full range for Neo: [1,2560,1,2160]
%Fitting Properties.
%[Gauss,TF] %[CxG,CyG,SxG,SyG,AmpG,ThetaG,Bg]
Fitting = 1;
FitInitialParameters = [311,222,31,29,0.5,-2.98,0];
FitFixedParameters = [0,0,1,1,0,0,0];
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
MarkerLocations = [309,221]; %[310,213]; %[311,218]; % [185,198]; %[178,195]; % [204.6,192];
%Debug Mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DebugMode = 0;

