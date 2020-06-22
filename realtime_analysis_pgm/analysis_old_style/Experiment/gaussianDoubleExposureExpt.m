function [MoreVals,MoreKeys,ImgCrop,ImgFit] = gaussianDoubleExposureExpt(Images,Vals,Keys,Settings,Constants)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Process Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set main ROI
X1 = Settings.WindowLimits(1);
X2 = Settings.WindowLimits(2);
Y1 = Settings.WindowLimits(3);
Y2 = Settings.WindowLimits(4);
ImagesCrop = Images(Y1:Y2,X1:X2,:);

if Settings.DebugMode
   displayROIs(Images,Settings);
end

%Remove Stripes
tic;
if Settings.RemoveStripes
    StrX1 = Settings.StripesRegionLimits(1);
    StrX2 = Settings.StripesRegionLimits(2);
    StrY1 = Settings.StripesRegionLimits(3);
    StrY2 = Settings.StripesRegionLimits(4);
    StripesNormRegion = Images(StrY1:StrY2,StrX1:StrX2,:);
    ImagesCrop = normalizeStripes(ImagesCrop,StripesNormRegion,Settings.StripesOrientation);  
end
TRemoveStripes = toc;
fprintf('Removing Stripes took %d \n',TRemoveStripes);

%Remove OD offsets...maybe could do this after binning?
tic;
if Settings.NormalizeOD
    X1Norm = Settings.NormRegionLimits(1);
    X2Norm = Settings.NormRegionLimits(2);
    Y1Norm = Settings.NormRegionLimits(3);
    Y2Norm = Settings.NormRegionLimits(4);
    ProcImages1 = getNormOD(ImagesCrop(:,:,[2,5,7]),Images(Y1Norm:Y2Norm,X1Norm:X2Norm,[2,5,7]));
    ProcImages2 = getNormOD(ImagesCrop(:,:,[3,6,7]),Images(Y1Norm:Y2Norm,X1Norm:X2Norm,[3,6,7]));
else
    ProcImages1 = ImagesCrop(:,:,[2,5,7]);
    ProcImages2 = ImagesCrop(:,:,[3,6,7]);
end
TNormOD = toc;
fprintf('Normalizing OD took %d \n',TNormOD);

%software bin
ProcImages1 = binImg(ProcImages1,Settings.SoftwareBinSizeH,Settings.SoftwareBinSizeV);
ProcImages2 = binImg(ProcImages2,Settings.SoftwareBinSizeH,Settings.SoftwareBinSizeV);

tic;
ImgCrop1 = getOD(ProcImages1);
ImgCrop2 = getOD(ProcImages2);
TOD = toc;
fprintf('Getting OD took %d \n',TOD);

tic;
DifferenceCrop = ImgCrop2-ImgCrop1;
ImgCrop = cat(3,ImgCrop1,ImgCrop2,DifferenceCrop);
TBinCrop = toc;
fprintf('Binning and cropping took %d \n',TBinCrop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Functions = {'gaussian2D'};
InitParams = Settings.FitInitialParameters;
FixedParams = Settings.FitFixedParameters;

if Settings.Fitting
try
    tic;
    Tol = 1e-5;
    Xs = X1:X2;
    Ys = Y1:Y2;
    [X,Y] = meshgrid(Xs,Ys);
    [FParams1,PNames1,FittedFunctionHandle1] = fit2D(X,Y,ImgCrop1,[],Functions,InitParams,FixedParams,Tol);
    CxG1 = FParams1(1); CyG1 = FParams1(2); SxG1 = FParams1(3); SyG1 = FParams1(4);
    AG1 = abs(FParams1(5)); ThetaG1 = FParams1(6); Bg1 = FParams1(7);
%     [SizeY,SizeX]=size(ImgCrop1);
%     [X,Y] = meshgrid(1:SizeX,1:SizeY);
    ImgFit1 = FittedFunctionHandle1(X,Y);
    
    [FParams2,PNames2,FittedFunctionHandle2] = fit2D(X,Y,ImgCrop2,[],Functions,InitParams,FixedParams,Tol);
    CxG2 = FParams2(1); CyG2 = FParams2(2); SxG2 = FParams2(3); SyG2 = FParams2(4);
    AG2 = abs(FParams2(5)); ThetaG2 = FParams2(6); Bg2 = FParams2(7);
%     [SizeY,SizeX]=size(ImgCrop2);
%     [X,Y] = meshgrid(1:SizeX,1:SizeY);
    ImgFit2 = FittedFunctionHandle2(X,Y);
    Tfit = toc;
    fprintf('Fitting took %d \n',Tfit);
    
    DifferenceFit = ImgFit2-Bg2-ImgFit1+Bg1;
    ImgFit = cat(3,ImgFit1,ImgFit2,DifferenceFit);
    
    
catch Exception
    MsgString = getReport(Exception);
    ImgFit1 = zeros(size(ImgCrop1));
    FParams1 = [0,0,0,0,0,0,0];
    PNames1 = repmat({'P'},1,length(FParams1));
    
    ImgFit2 = zeros(size(ImgCrop2));
    FParams2 = [0,0,0,0,0,0,0];
    PNames2 = repmat({'P'},1,length(FParams2));
    
    DifferenceFit = zeros(size(ImgCrop1));
    
    ImgFit = cat(3,ImgFit1,ImgFit2,DifferenceFit);
    disp('Error performing fit.')
    disp(MsgString);
end
else
    ImgFit1 = zeros(size(ImgCrop1));
    FParams1 = [0,0,0,0,0,0,0];
    PNames1 = repmat({'P'},1,length(FParams1));
    
    ImgFit2 = zeros(size(ImgCrop2));
    FParams2 = [0,0,0,0,0,0,0];
    PNames2 = repmat({'P'},1,length(FParams2));
    
    DifferenceFit = zeros(size(ImgCrop1));
    
    ImgFit = cat(3,ImgFit1,ImgFit2,DifferenceFit);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis From Fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
try
    %These ignore software binning...
    AtomAvg1 = mean(mean(ImagesCrop(:,:,2)));
    AtomAvg2 = mean(mean(ImagesCrop(:,:,3)));
    BeamAvg1 = mean(mean(ImagesCrop(:,:,5)));
    BeamAvg2 = mean(mean(ImagesCrop(:,:,6)));
    BeamSD1 = std2(ImagesCrop(:,:,3));
    BeamSD2 = std2(ImagesCrop(:,:,5));
    DarkAvg = mean(mean(ImagesCrop(:,:,7)));
    
    PixelArea = Settings.getEffPixelArea;
    Lambda = Constants.lambda_D2;
    QE = Settings.QuantumEfficiency;
    CountsPerPhoton = Settings.CountsPerPhoton;
    ImgT = Settings.ImagingDuration;
    IoverIsat1 = getIntensity(BeamAvg1,ImgT,Lambda,QE,CountsPerPhoton,PixelArea)/Constants.Isat_D2;
    IoverIsat2 = getIntensity(BeamAvg2,ImgT,Lambda,QE,CountsPerPhoton,PixelArea)/Constants.Isat_D2;
    
    CountsAtFitOD1 = (BeamAvg1-DarkAvg)*exp(-AG1);
    CountsAtFitOD2 = (BeamAvg2-DarkAvg)*exp(-AG2);
    
    NSum1 = getAtomNumber(sum(sum(ImgCrop1)),0,PixelArea,Constants);
    NSum2 = getAtomNumber(sum(sum(ImgCrop2)),0,PixelArea,Constants);
        
    NFit2 = getAtomNumber(abs(2*pi*SxG2*SyG2*AG2),0,PixelArea,Constants);
    NFit1 = getAtomNumber(abs(2*pi*SxG1*SyG1*AG1),0,PixelArea,Constants);
        
    NPix = Size(ImgCrop1,2)*Size(ImgCrop1,1);
    NSumNoBg1 = getAtomNumber(sum(sum(ImgCrop1))-NPix*Bg1,0,PixelArea,Constants);
    NSumNoBg2 = getAtomNumber(sum(sum(ImgCrop2))-NPix*Bg2,0,PixelArea,Constants);

catch Exception
    if(~exist('IoverIsat1','var'))
        IoverIsat1 = 0;
    end
    if(~exist('IoverIsat2','var'))
        IoverIsat2 = 0;
    end    
    MsgString = getReport(Exception);
    disp(MsgString);
end

    
try
    TOF = Vals(strcmp(Keys,'TOF')==1);
    if isempty(TOF)
        throw(MException('TOF:isempty',''));
    end
    Tx1 = getTempSingle(SxG1,TOF,PixelArea,Constants);
    Ty1 = getTempSingle(SyG1,TOF,PixelArea,Constants);
    Tx2 = getTempSingle(SxG2,TOF,PixelArea,Constants);
    Ty2 = getTempSingle(SyG2,TOF,PixelArea,Constants);
    
catch Exception
    if(~exist('NFit1','var'))
        NFit1 = 0;
    end
    if(~exist('NSum1','var'))
        NSum1 = 0;
    end
    if(~exist('NFit2','var'))
        NFit2 = 0;
    end
    if(~exist('NSum2','var'))
        NSum2 = 0;
    end
    if(~exist('NSumNoBg1','var'))
        NSumNoBg1 = 0;
    end
    if(~exist('NSumNoBg2','var'))
        NSumNoBg2 = 0;
    end
    if(~exist('CountsAtFitOD1','var'))
        CountsAtFitOD1 = 0;
    end
    if(~exist('CountsAtFitOD2','var'))
        CountsAtFitOD2 = 0;
    end
    
    
    
    
    Tx1 = 0; Ty1 = 0; Tx2 = 0; Ty2 = 0;
    MsgString = getReport(Exception);
    disp('Analysis failed.')
    disp(MsgString);
end

try
    Polarization = -(NFit1-NFit2)/(NFit1+NFit2);
catch Exception
    Polarization = 0;
    MsgString = getReport(Exception);
    disp('Analysis failed.')
    disp(MsgString);
end

TInfo = toc;
fprintf('Getting atom# etc. took %d \n',TInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Organize Results and Return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MoreVals = [FParams1,NSum1,NSumNoBg1,NFit1,Tx1,Ty1,FParams2,NSum2,NSumNoBg2,NFit2,Tx2,Ty2,...
    Polarization,AtomAvg1,AtomAvg2,BeamAvg1,BeamAvg2,IoverIsat1,IoverIsat2,BeamSD1,BeamSD2,CountsAtFitOD1,CountsAtFitOD2,DarkAvg];
MoreKeys = cat(2,PNames1,{'NSum1','NSumNoBG1','NFit1','Tx1','Ty1'},PNames2,{'NSum2','NSumNoBg2','NFit2','Tx2','Ty2'},...
    {'Polarization','Atoms1 AvgCts','Atoms2 AvgCts','Beam1 AvgCts','Beam2 AvgCts','IoverIsat1','IoverIsat2','BeamSD1','BeamSD2','Estimated Counts At Fit OD1','Estimated Counts At Fit OD2','Dark AvgCts'});
end

