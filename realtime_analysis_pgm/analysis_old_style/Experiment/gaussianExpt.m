function [MoreVals,MoreKeys,ImgCrop,ImgFit] = gaussianExpt(Images,Vals,Keys,Settings,Constants)


X1 = Settings.WindowLimits(1);
X2 = Settings.WindowLimits(2);
Y1 = Settings.WindowLimits(3);
Y2 = Settings.WindowLimits(4);
ImagesCrop = Images(Y1:Y2,X1:X2,:);

if(size(ImagesCrop,3)==3)
    Indices = [1,2,3];
elseif (size(ImagesCrop,3)==4)
    Indices = [2,3,4];
else
    disp('Unsupported image size in gaussianExpt. Only accepts 3 or 4 images.');
end

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
    ProcImages = getNormOD(ImagesCrop(:,:,Indices),Images(Y1Norm:Y2Norm,X1Norm:X2Norm,Indices));
else
    ProcImages = ImagesCrop(:,:,Indices);
end
TNormOD = toc;
fprintf('Normalizing OD took %d \n',TNormOD);

%software bin
ProcImages = binImg(ProcImages,Settings.SoftwareBinSizeH,Settings.SoftwareBinSizeV);

tic;
ImgCrop = getOD(ProcImages);
TOD = toc;
fprintf('Getting OD took %d \n',TOD);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
if Settings.Fitting
    
    Functions = {'gaussian2D'};
    InitParams = Settings.FitInitialParameters;
    FixedParams = Settings.FitFixedParameters;
    if Settings.UseFitAutoGuess
       [CxAuto,CyAuto,SxAuto,SyAuto] = guessCenter(ImgCrop);
       AutoGuesses = [CxAuto,CyAuto,SxAuto,SyAuto];
       for ii = 1:4
        if ~FixedParams(ii)
            InitParams(ii) = AutoGuesses(ii);
        end
       end
    end
    
    
    %Fitting
    try
        Tol = 1e-5;
        Xs = X1:X2;
        Ys = Y1:Y2;
        [X,Y] = meshgrid(Xs,Ys);
        [FParams,PNames,FittedFunctionHandle] = fit2D(X,Y,ImgCrop,[],Functions,InitParams,FixedParams,Tol);
%         [SizeY,SizeX]=size(ImgCrop);
%         [X,Y] = meshgrid(1:SizeX,1:SizeY);
        ImgFit = FittedFunctionHandle(X,Y);
        
    catch Exception
        MsgString = getReport(Exception);
        ImgFit = zeros(size(ImgCrop));
        FParams = [0,0,0,0,0,0,0];
        PNames = repmat({'P'},1,length(FParams));
        disp('Error performing fit.')
        disp(MsgString);
    end
else
    ImgFit = zeros(size(ImgCrop));
    FParams = [0,0,0,0,0,0,0];
    PNames = repmat({'P'},1,length(FParams));
end
%Name variables for later use.
CxG = FParams(1); CyG = FParams(2); SxG = FParams(3); SyG = FParams(4);
AG = FParams(5); ThetaG = FParams(6); Bg = FParams(7);
Tfit = toc;
fprintf('Fitting took %d \n',Tfit);
%Analysis
try
    %Use unprocessed images for beam params.
    if(size(ImagesCrop,3)==4)
        AtomAvg = mean(mean(ImagesCrop(:,:,2)));
        BeamAvg = mean(mean(ImagesCrop(:,:,3)));
        DarkAvg = mean(mean(ImagesCrop(:,:,4)));
        BeamSD = std2(ImagesCrop(:,:,3));
    elseif(size(ImagesCrop,3)==3)
        AtomAvg = mean(mean(ImagesCrop(:,:,1)));
        BeamAvg = mean(mean(ImagesCrop(:,:,2)));
        DarkAvg = mean(mean(ImagesCrop(:,:,3)));
        BeamSD = std2(ImagesCrop(:,:,2));
    end
        
    
    PixelArea = Settings.getEffPixelArea;
    Lambda = Constants.lambda_D2;
    QE = Settings.QuantumEfficiency;
    CountsPerPhoton = Settings.CountsPerPhoton;
    ImgT = Settings.ImagingDuration;
    IoverIsat = getIntensity(BeamAvg-DarkAvg,ImgT,Lambda,QE,CountsPerPhoton,PixelArea)/Constants.Isat_D2;
    CountsAtFitOD = (BeamAvg-DarkAvg)*exp(-AG);
    
    NFit = getAtomNumber(2*pi*SxG*SyG*AG,0,PixelArea,Constants);
    NSum = getAtomNumber(sum(sum(ImgCrop)),0,PixelArea,Constants);
     
    NPix = size(ImgCrop,1)*size(ImgCrop,2);
    NSumNoBg = getAtomNumber(sum(sum(ImgCrop))-NPix*Bg,0,PixelArea,Constants);
    
catch Exception
    if(~exist('IoverIsat','var'))
        IoverIsat = 0;
    end
    if(~exist('CountsAtFitOD','var'))
        CountsAtFitOD = 0;
    end
    if(~exist('NFit','var'))
        NFit = 0;
    end
    if(~exist('NSum','var'))
        NSum = 0;
    end
    if(~exist('NSumNoBg','var'))
        NSumNoBg = 0;
    end
    
    MsgString = getReport(Exception);
    disp(MsgString);
end

MoreVals = [FParams,NSum,NSumNoBg,NFit,AtomAvg,BeamAvg,IoverIsat,BeamSD,CountsAtFitOD,DarkAvg];
MoreKeys = cat(2,PNames,{'NSum','NSumNoBg','NFit','Atom AvgCts','Beam AvgCts','IOverISat','BeamSD','Estimated Counts At Fit OD','Dark AvgCts'});
end

