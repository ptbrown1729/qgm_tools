function [MoreVals,MoreKeys,ImgCrop,ImgFit] = doubleExposureCondFrac(Images,Vals,Keys,Settings,Constants)
%[MoreVals,MoreKeys,ImgCrop,ImgFit] = doubleExposureCondFrac(Images,Vals,Keys,Settings,Constants)

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

%Remove OD offsets
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


tic;
ImgCrop1NoBin = getOD(ProcImages1);
ImgCrop2NoBin = getOD(ProcImages2);
TOD = toc;
fprintf('Getting OD took %d \n',TOD);


%software bin
tic;
ImgCrop1 = binImg(ImgCrop1NoBin,Settings.SoftwareBinSizeH,Settings.SoftwareBinSizeV);
ImgCrop2 = binImg(ImgCrop2NoBin,Settings.SoftwareBinSizeH,Settings.SoftwareBinSizeV);
DifferenceCrop = ImgCrop1-ImgCrop2;
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
    %2D Fitting
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
%     [X,Y] = meshgrid(1:SizeX,1:SizeY);g
    ImgFit2 = FittedFunctionHandle2(X,Y);
    Tfit = toc;
    fprintf('Fitting took %d \n',Tfit);
    
    DifferenceFit = ImgFit1-Bg1-ImgFit2+Bg2;
    ImgFit = cat(3,ImgFit1,ImgFit2,DifferenceFit);
    
    %1D Fitting
    Profile1 = sum(ImgCrop1,1);
    XPts1 = 1:length(Profile1);
    Profile2 = sum(ImgCrop2,1);
    XPts2 = 1:length(Profile2);
    InitParams1 = [length(Profile1)/2,length(Profile1)/4,0.3,0,length(Profile1)/2,length(Profile1)/4,0.3,0];
    [FOneDParams1,PNamesOneD1,FittedFunctionHandleOneD1] = fit1D(XPts1,Profile1,{'gaussian1D','thomasFermi1D'},InitParams1,zeros(1,8),1e-5);    
    SG1 = FOneDParams1(2); AG1 = FOneDParams1(3);
    RTF1 = FOneDParams1(6); ATF1 = FOneDParams1(7);
    InitParams2 = [length(Profile2)/2,length(Profile2)/4,0.3,0,length(Profile2)/2,length(Profile2)/4,0.3,0];
    [FOneDParams2,PNamesOneD2,FittedFUnctionHandleOneD2] = fit1D(XPts2,Profile2,{'gaussian1D','thomasFermi1D'},InitParams2,zeros(1,8),1e-5);
    SG2 = FOneDParams2(2); AG2 = FOneDParams2(3);
    RTF2 = FOneDParams2(6); ATF2 = FOneDParams2(7);    

    %plot condensate fraction parts
    FigureName = 'Condensate Fraction';
    fhandle = findobj('type','figure','name',FigureName);
    if isempty(fhandle)
        figure('name',FigureName);
    elseif length(fhandle)>1
        close(fhandle)
        figure('name',FigureName);
    else
        figure(fhandle);
    end
    
    subplot(2,1,1)
    plot(XPts1,Profile1,'b')
    hold on;
    plot(XPts1,FittedFunctionHandleOneD1(XPts1),'r')
    plot(XPts1,gaussian1D(FOneDParams1(1:4),XPts1,0),'black')
    plot(XPts1,thomasFermi1D(FOneDParams1(5:8),XPts1,0),'g')
    hold off;
    subplot(2,1,2)
    plot(XPts2,Profile2,'b')
    hold on;
    plot(XPts2,FittedFUnctionHandleOneD2(XPts2),'r')
    plot(XPts2,gaussian1D(FOneDParams2(1:4),XPts2,0),'black')
    plot(XPts2,thomasFermi1D(FOneDParams2(5:8),XPts2,0),'g')
    hold off;
    
    
catch Exception
    MsgString = getReport(Exception);
    ImgFit1 = zeros(size(ImgCrop1));
    FParams1 = zeros(1,7);
    PNames1 = repmat({'P'},1,length(FParams1));
    
    ImgFit2 = zeros(size(ImgCrop2));
    FParams2 = zeros(1,7);
    PNames2 = repmat({'P'},1,length(FParams2));
    
    DifferenceFit = zeros(size(ImgCrop1));
    
    ImgFit = cat(3,ImgFit1,ImgFit2,DifferenceFit);
    
    FOneDParams1 = zeros(1,8);
    PNamesOneD1 = repmat({'P'},1,length(FOneDParams1));
    
    FOneDParams2 = zeros(1,8);
    PNamesOneD2 = repmat({'P'},1,length(FOneDParams2));
    
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

%Analysis From Fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
try
    AtomAvg1 = mean(mean(ImagesCrop(:,:,2)));
    AtomAvg2 = mean(mean(ImagesCrop(:,:,3)));
    BeamAvg1 = mean(mean(ImagesCrop(:,:,5)));
    BeamAvg2 = mean(mean(ImagesCrop(:,:,6)));
    DarkAvg = mean(mean(ImagesCrop(:,:,7)));
    
    PixelArea = Settings.getEffPixelArea;
    
    NSum1 = getAtomNumber(sum(sum(ImgCrop1)),0,PixelArea,Constants);
    NSum2 = getAtomNumber(sum(sum(ImgCrop2)),0,PixelArea,Constants);
        
    NFit2 = getAtomNumber(abs(2*pi*SxG2*SyG2*AG2),0,PixelArea,Constants);
    NFit1 = getAtomNumber(abs(2*pi*SxG1*SyG1*AG1),0,PixelArea,Constants);
        
    NPix = size(ImgCrop1,2)*size(ImgCrop1,1);
    NSumNoBg1 = getAtomNumber(sum(sum(ImgCrop1))-NPix*Bg1,0,PixelArea,Constants);
    NSumNoBg2 = getAtomNumber(sum(sum(ImgCrop2))-NPix*Bg2,0,PixelArea,Constants);
    
    G1Area = sqrt(2*pi)*SG1*AG1;
    TF1Area = 4/3*ATF1*RTF1;
    CondFrac1 = TF1Area/(G1Area+TF1Area);
    G2Area = sqrt(2*pi)*SG2*AG2;
    TF2Area = 4/3*ATF2*RTF2;
    CondFrac2 = TF2Area/(G2Area+TF2Area);
    
    if NFit1 == NFit2
        Polarization = 0;
    else
        Polarization = abs((NFit1-NFit2)/(NFit1+NFit2));
    end
    
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
    
%     NSumNoBg1 = 0;
%     NSumNoBg2 = 0;
%     Polarization = 0;
    
    MsgString = getReport(Exception);
    Tx1 = 0;
    Ty1 = 0;
    Tx2 = 0;
    Ty2 = 0;
    disp('Analysis failed.')
    disp(MsgString);
end
TInfo = toc;
fprintf('Getting atom# etc. took %d \n',TInfo);

%Organize Results and Return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MoreVals = [FParams1,NSum1,NSumNoBg1,NFit1,Tx1,Ty1,FParams2,NSum2,...
    NSumNoBg2,NFit2,Tx2,Ty2,Polarization,AtomAvg1,AtomAvg2,BeamAvg1,...
    BeamAvg2,DarkAvg,FOneDParams1,FOneDParams2,CondFrac1,CondFrac2];
MoreKeys = cat(2,PNames1,{'NSum1','NSumNoBG1','NFit1','Tx1','Ty1'}...
    ,PNames2,{'NSum2','NSumNoBg2','NFit2','Tx2','Ty2'},...
    {'Polarization','Atoms1 AvgCts','Atoms2 AvgCts','Beam1 AvgCts',...
    'Beam2 AvgCts','Dark AvgCts'},PNamesOneD1,PNamesOneD2,{'CondFrac1','CondFrac2'});
end


