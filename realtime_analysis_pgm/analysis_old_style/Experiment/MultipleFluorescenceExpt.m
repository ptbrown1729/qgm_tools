function [MoreVals,MoreKeys,ImgCrop,ImgFit] = FluorescenceExpt(Images,Vals,Keys,Settings,Constants)


X1 = Settings.WindowLimits(1);
X2 = Settings.WindowLimits(2);
Y1 = Settings.WindowLimits(3);
Y2 = Settings.WindowLimits(4);
ImgCrop = Images(Y1:Y2,X1:X2,:);
NImgs = size(ImgCrop,3);
ImgCrop = binImg(ImgCrop,Settings.SoftwareBinSizeH,Settings.SoftwareBinSizeV);

if Settings.DebugMode
   displayROIs(Images,Settings);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
if Settings.Fitting
    
    Functions = {'gaussian2D'};
    InitParams = Settings.FitInitialParameters;
    FixedParams = Settings.FitFixedParameters;
    
    %Fitting
    for ii = 1:NImgs
        try
            Tol = 1e-5;
            Xs = X1:X2;
            Ys = Y1:Y2;
            [X,Y] = meshgrid(Xs,Ys);
            [FParams,PNames,FittedFunctionHandle] = fit2D(X,Y,ImgCrop(:,:,ii),[],Functions,InitParams,FixedParams,Tol);
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
    BeamAvg = mean(mean(ImgCrop));
    BeamSD = std2(ImgCrop);
    CountSum = sum(sum(ImgCrop));
    CountSumNoBg = CountSum - size(ImgCrop,1)*size(ImgCrop,2)*Bg;
    CountFit = 2*pi*SxG*SyG*AG;
    
    PixelArea = Settings.getEffPixelArea;
    Lambda = Constants.lambda_D2;
    QE = Settings.QuantumEfficiency;
    CountsPerPhoton = Settings.CountsPerPhoton;
    ImgT = Settings.ImagingDuration;
    IoverIsat = getIntensity(BeamAvg,ImgT,Lambda,QE,CountsPerPhoton,PixelArea)/Constants.Isat_D2;

    
catch Exception
    if(~exist('IoverIsat','var'))
        IoverIsat = 0;
    end

  
    MsgString = getReport(Exception);
    disp(MsgString);
end

MoreVals = [FParams,CountSum,CountSumNoBg,CountFit,BeamAvg,IoverIsat,BeamSD];
MoreKeys = cat(2,PNames,{'CountSum','CountSumNoBg','CountFit','Beam AvgCts','IOverISat','BeamSD'});
end

