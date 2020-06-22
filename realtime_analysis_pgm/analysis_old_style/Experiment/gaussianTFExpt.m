function [MoreVals,MoreKeys,ImgCrop,ImgFit] = gaussianTFExpt(Images,Vals,Keys,Settings,Constants)

X1 = Settings.WindowLimits(1);
X2 = Settings.WindowLimits(2);
Y1 = Settings.WindowLimits(3);
Y2 = Settings.WindowLimits(4);

ImagesCrop = Images(Y1:Y2,X1:X2,:);
ImagesCrop = binImg(ImagesCrop,Settings.SoftwareBinSizeH,Settings.SoftwareBinSizeV);
ImgCrop = getOD(ImagesCrop);


Functions = {'gaussian2D','thomasFermi2D'};
InitParams = Settings.FitInitialParameters;
FixedParams = Settings.FitFixedParameters;

%Fitting
try
    
    Tol = 1e-5;
    Xs = X1:X2;
    Ys = Y1:Y2;
    [X,Y] = meshgrid(Xs,Ys);
    [FParams,PNames,FittedFunctionHandle] = fit2D(X,Y,ImgCrop,[],Functions,InitParams,FixedParams,Tol);
    CxG = FParams(1); CyG = FParams(2); SxG = FParams(3); SyG = FParams(4); 
    AG = abs(FParams(5)); ThetaG = FParams(6);
    CxTF = FParams(8); CyTF = FParams(9); RxTF = FParams(10); RyTF = FParams(11);
    ATF = abs(FParams(12)); ThetaTF = FParams(13);    
    Bg = FParams(7) + FParams(14);
%     [SizeY,SizeX]=size(ImgCrop);
%     [X,Y] = meshgrid(1:SizeX,1:SizeY);
    ImgFit = FittedFunctionHandle(X,Y);        
    
    
catch Exception
    MsgString = getReport(Exception);
    ImgFit = zeros(size(ImgCrop));
    FParams = [[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]];
    PNames = repmat({'P'},1,length(FParams));
    disp('Error performing fit.')
    disp(MsgString);
end

%Analysis
try
    AtomAvg = mean(mean(ImagesCrop(:,:,1)));
    BeamAvg = mean(mean(ImagesCrop(:,:,2)));
    DarkAvg = mean(mean(ImagesCrop(:,:,3)));
    
    PixelArea = Settings.getEffPixelArea;
    NFitG = getAtomNumber(2*pi*SxG*SyG*AG,0,PixelArea,Constants);
    NFitTF = getAtomNumber(0.5*pi*RxTF*RyTF*ATF,0,PixelArea,Constants);
    NSum = getAtomNumber(sum(sum(ImgCrop)),0,PixelArea,Constants);
    TOF = Vals(strcmp(Keys,'TOF')==1);
    if isempty(TOF)
        throw(MException('TOF:isempty',''));
    end
    TxG = getTempSingle(SxG,TOF,PixelArea,Constants);
    TyG = getTempSingle(SyG,TOF,PixelArea,Constants);
    
catch Exception
    if(~exist('NFitG','var'))
        NFitG = 0;
    end
    if(~exist('NSum','var'))
        NSum = 0;
    end
    if(~exist('NFitTF','var'))
        NFitTF = 0;
    end
    
    MsgString = getReport(Exception);
    TxG = 0;
    TyG = 0;
    disp('Analysis failed.')
    disp(MsgString);
end

MoreVals = [FParams,NSum,NFitTF, NFitG,TxG,TyG,AtomAvg,BeamAvg,DarkAvg];
MoreKeys = cat(2,PNames,{'NSum','NFitTF','NFitG','Tx','Ty','Atom AvgCts','Beam AvgCts','Dark AvgCts'});
end

