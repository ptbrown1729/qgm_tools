function [MoreVals,MoreKeys,ImgCrop,ImgFit] = gaussianTFDoubleExposureExpt(Images,Vals,Keys,Settings,Constants)

X1 = Settings.WindowLimits(1);
X2 = Settings.WindowLimits(2);
Y1 = Settings.WindowLimits(3);
Y2 = Settings.WindowLimits(4);

%crop image
ImagesCrop = Images(Y1:Y2,X1:X2,:);

%remove stripes
if Settings.RemoveStripes
    ImagesCrop = normalizeStripes(ImagesCrop,Images(Y1:Y2,300:400,:));
    %ImagesCrop = normalizeStripes(ImagesCrop,Images(Y1:Y2,X2:X2+80,:));  
end
%remove OD offsets
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

%software bin
ProcImages1 = binImg(ProcImages1,Settings.SoftwareBinSizeH,Settings.SoftwareBinSizeV);
ProcImages2 = binImg(ProcImages2,Settings.SoftwareBinSizeH,Settings.SoftwareBinSizeV);

ImgCrop1 = getOD(ProcImages1);
ImgCrop2 = getOD(ProcImages2);


DifferenceCrop = ImgCrop1-ImgCrop2;
ImgCrop = cat(3,ImgCrop1,ImgCrop2,DifferenceCrop);

Functions = {'gaussian2D','thomasFermi2D'};
InitParams = Settings.FitInitialParameters;
FixedParams = Settings.FitFixedParameters;

%Fitting
try
    
    Tol = 1e-5;
    Xs = X1:X2;
    Ys = Y1:Y2;
    [X,Y] = meshgrid(Xs,Ys);
    [FParams1,PNames1,FittedFunctionHandle1] = fit2D(X,Y,ImgCrop1,[],Functions,InitParams,FixedParams,Tol);
    CxG1 = FParams1(1); CyG1 = FParams1(2); SxG1 = FParams1(3); SyG1 = FParams1(4);
    AG1 = abs(FParams1(5)); ThetaG1 = FParams1(6);
    CxTF1 = FParams1(8); CyTF1 = FParams1(9); RxTF1 = FParams1(10); RyTF1 = FParams1(11);
    ATF1 = abs(FParams1(12)); ThetaTF1 = FParams1(13);    
    Bg1 = FParams1(7) + FParams1(14);
%     [SizeY,SizeX]=size(ImgCrop1);
%     [X,Y] = meshgrid(1:SizeX,1:SizeY);
    ImgFit1 = FittedFunctionHandle1(X,Y);
    
    [FParams2,PNames2,FittedFunctionHandle2] = fit2D(X,Y,ImgCrop2,[],Functions,InitParams,FixedParams,Tol);
    CxG2 = FParams2(1); CyG2 = FParams2(2); SxG2 = FParams2(3); SyG2 = FParams2(4);
    AG2 = abs(FParams2(5)); ThetaG2 = FParams2(6);
    CxTF2 = FParams2(8); CyTF2 = FParams2(9); RxTF2 = FParams2(10); RyTF2 = FParams2(11);
    ATF2 = abs(FParams2(12)); ThetaTF2 = FParams2(13);    
    Bg2 = FParams2(7) + FParams2(14);
%     [SizeY,SizeX]=size(ImgCrop2);
%     [X,Y] = meshgrid(1:SizeX,1:SizeY);
    ImgFit2 = FittedFunctionHandle2(X,Y);
    
    DifferenceFit = ImgFit1-ImgFit2;
    ImgFit = cat(3,ImgFit1,ImgFit2,DifferenceFit);
    
catch Exception
    MsgString = getReport(Exception);
    ImgFit1 = zeros(size(ImgCrop1));
    FParams1 = [0,0,0,0,0,0,0];
    PNames1 = repmat({'P'},1,length(FParams1));
    
    ImgFit2 = zeros(size(ImgCrop2));
    FParams2 = [0,0,0,0,0,0,0];
    PNames2 = repmat({'P'},1,length(FParams2));
    
    ImgFit = cat(3,ImgFit1,ImgFit2);
    disp('Error performing fit.')
    disp(MsgString);
end

%Analysis
try
    AtomAvg1 = mean(mean(ImagesCrop(:,:,2)));
    AtomAvg2 = mean(mean(ImagesCrop(:,:,3)));
    BeamAvg1 = mean(mean(ImagesCrop(:,:,5)));
    BeamAvg2 = mean(mean(ImagesCrop(:,:,6)));
    DarkAvg = mean(mean(ImagesCrop(:,:,7)));
    
    PixelArea = Settings.getEffPixelArea;
    NFitG1 = getAtomNumber(abs(2*pi*SxG1*SyG1*AG1),0,PixelArea,Constants);
    NFitTF1 = getAtomNumber(0.5*pi*RxTF1*RyTF1*ATF1,0,PixelArea,Constants);
    NSum1 = getAtomNumber(sum(sum(ImgCrop1)),0,PixelArea,Constants);
    NFitG2 = getAtomNumber(abs(2*pi*SxG2*SyG2*AG2),0,PixelArea,Constants);
    NFitTF2 = getAtomNumber(0.5*pi*RxTF2*RyTF2*ATF2,0,PixelArea,Constants);
    NSum2 = getAtomNumber(sum(sum(ImgCrop2)),0,PixelArea,Constants);
    
    if NFitTF1 == NFitTF2
        Polarization = 0;
    else
        Polarization = abs((NFitTF1-NFitTF2)/(NFitTF1+NFitTF2));
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
%     if(~exist('NFit1','var'))
%         NFitG1 = 0;
%     end
    if(~exist('NSum1','var'))
        NSum1 = 0;
    end
%     if(~exist('NFit2','var'))
%         NFitG2 = 0;
%     end
    if(~exist('NSum2','var'))
        NSum2 = 0;
    end
    
    MsgString = getReport(Exception);
    Tx1 = 0;
    Ty1 = 0;
    Tx2 = 0;
    Ty2 = 0;
%     disp('Analysis failed.')
%     disp(MsgString);
end

MoreVals = [FParams1,NSum1,NFitG1,NFitTF1,Tx1,Ty1,FParams2,NSum2,NFitG2,NFitTF2,Tx2,Ty2,Polarization,AtomAvg1,AtomAvg2,BeamAvg1,BeamAvg2,DarkAvg];
MoreKeys = cat(2,PNames1,{'NSum1','NFitG1','NFitTF1','Tx1','Ty1'},PNames2,{'NSum2','NFitG2','NFitTF2','Tx2','Ty2'},{'Polarization','Atoms1 AvgCts','Atoms2 AvgCts','Beam1 AvgCts','Beam2 AvgCts','Dark AvgCts'});
end

