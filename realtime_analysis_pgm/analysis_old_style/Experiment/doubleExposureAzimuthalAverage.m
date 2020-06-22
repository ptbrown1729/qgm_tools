function [MoreVals,MoreKeys,ImgCrop,ImgFit] = doubleExposureAzimuthalAverage(Images,Vals,Keys,Settings,Constants)
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

PixelArea = Settings.getEffPixelArea;

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
% size(ProcImages1)
TNormOD = toc;
fprintf('Normalizing OD took %d \n',TNormOD);

%software bin
tic;
BinImg1 = getOD(binImg(ProcImages1,Settings.SoftwareBinSizeH,Settings.SoftwareBinSizeV,'Normal'));
BinImg2 = getOD(binImg(ProcImages2,Settings.SoftwareBinSizeH,Settings.SoftwareBinSizeV,'Normal'));
% figure
% imagesc(BinImg1)

ImgCrop1 = getOD(binImg(ProcImages1,Settings.SoftwareBinSizeH,Settings.SoftwareBinSizeV));
ImgCrop2 = getOD(binImg(ProcImages2,Settings.SoftwareBinSizeH,Settings.SoftwareBinSizeV));

ImgCropNoBin1 = getOD(ProcImages1);
ImgCropNoBin2 = getOD(ProcImages2);


DifferenceCrop = ImgCrop1-ImgCrop2;
ImgCrop = cat(3,ImgCrop1,ImgCrop2,DifferenceCrop);
TBinCrop = toc;
fprintf('Binning, cropping and getting OD took %d \n',TBinCrop);

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
%         [SizeY,SizeX]=size(ImgCrop1);
%         [X,Y] = meshgrid(1:SizeX,1:SizeY);
        ImgFit1 = FittedFunctionHandle1(X,Y);
        
        [FParams2,PNames2,FittedFunctionHandle2] = fit2D(X,Y,ImgCrop2,[],Functions,InitParams,FixedParams,Tol);
        CxG2 = FParams2(1); CyG2 = FParams2(2); SxG2 = FParams2(3); SyG2 = FParams2(4);
        AG2 = abs(FParams2(5)); ThetaG2 = FParams2(6); Bg2 = FParams2(7);
%         [SizeY,SizeX]=size(ImgCrop2);
%         [X,Y] = meshgrid(1:SizeX,1:SizeY);
        ImgFit2 = FittedFunctionHandle2(X,Y);
        Tfit = toc;
        fprintf('Fitting took %d \n',Tfit);
        
        DifferenceFit = ImgFit1-Bg1-ImgFit2+Bg2;
        ImgFit = cat(3,ImgFit1,ImgFit2,DifferenceFit);
        %
        %     %plot sum
        %     FigureName = 'Sum';
        %     fhandle = findobj('type','figure','name',FigureName);
        %     if isempty(fhandle)
        %         figure('name',FigureName);
        %     elseif length(fhandle)>1
        %         close(fhandle)
        %         figure('name',FigureName);
        %     else
        %         figure(fhandle);
        %     end
        %
        %     subplot(2,1,1)
        %     plot(sum(ImgCropNoBin1,1));
        %     subplot(2,1,2)
        %     plot(sum(ImgCropNoBin2,1));
        %
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
        
        disp('Error performing fit.')
        disp(MsgString);
    end
    
    try
        %Make 2D Grid
        [SizeY,SizeX]=size(BinImg1);
        Xvec=1:1:SizeX;
        Yvec=1:1:SizeY;
        Grid1=zeros(SizeY,SizeX);
        Grid2=zeros(SizeY,SizeX);
        GridParams = Settings.GridParams;
        GridFixedParams = Settings.GridFixedParams;
        if GridFixedParams(1)
            CxG1=GridParams(1); CxG2=GridParams(1);
        end
        if GridFixedParams(2)
            CyG1=GridParams(2); CyG2=GridParams(2);
        end
        if GridFixedParams(3)
            ThetaG1=GridParams(3); ThetaG2=GridParams(3);
        end
        %     CxG1=180;CyG1=150;CxG2=180;CyG2=150;
        %     ThetaG1=0;ThetaG2=0;
        for i=1:SizeX
            for j = 1:SizeY
                Grid1(j,i)=sqrt(((Xvec(i)-CxG1)*cos(ThetaG1)+(Yvec(j)-CyG1)*sin(ThetaG1))^2....
                    +(-(Xvec(i)-CxG1)*sin(ThetaG1)+(Yvec(j)-CyG1)*cos(ThetaG1))^2);
                Grid2(j,i)=sqrt(((Xvec(i)-CxG2)*cos(ThetaG2)+(Yvec(j)-CyG2)*sin(ThetaG2))^2....
                    +(-(Xvec(i)-CxG2)*sin(ThetaG2)+(Yvec(j)-CyG2)*cos(ThetaG2))^2);
            end
        end
        
        %1D Azimuthal Average
        AzInitParams = Settings.AzimuthalAverageInitialParameters;
        AzFixedParams = Settings.AzimuthalAverageFixedParameters;
        
        [bins1,AzAvg1,AzAvgErr1,AzParams1,AzPNames1,Comb1DFit1,GFit1,TF1Fit] = AzimuthalAverage(BinImg1,Grid1,GridParams(4),GridParams(5),AzInitParams,AzFixedParams);
        [bins2,AzAvg2,AzAvgErr2,AzParams2,AzPNames2,Comb1DFit2,GFit2,TF2Fit] = AzimuthalAverage(BinImg2,Grid2,GridParams(4),GridParams(5),AzInitParams,AzFixedParams);
        
        %     [bins1,AzAvg1,AzAvgErr1,AzParams1,AzPNames1,Comb1DFit1,GFit1,TF1Fit] = AzimuthalAverage(ImgCrop1,Grid1,GridParams(4),GridParams(5),AzInitParams,AzFixedParams);
        %     [bins2,AzAvg2,AzAvgErr2,AzParams2,AzPNames2,Comb1DFit2,GFit2,TF2Fit] = AzimuthalAverage(ImgCrop2,Grid2,GridParams(4),GridParams(5),AzInitParams,AzFixedParams);
        
        %     [FPDbleG1,PG1,FhDG1] = fit1Derrorbars(bins1,AzAvg1,AzAvgErr1,{'gaussian1D','gaussian1D'},[0,10,0.3,0,0,50,0.8,0],[1,0,0,0,1,0,0,1],1e-5);
        %     [FPDbleG2,PG2,FhDG2] = fit1Derrorbars(bins2,AzAvg2,AzAvgErr2,{'gaussian1D','gaussian1D'},[0,100,0.5,0,0,50,0.8,0],[1,0,0,0,1,0,0,1],1e-5);
        %
        %   plot condensate fraction parts
        FigureName = 'Azimuthal Average';
        fhandle = findobj('type','figure','name',FigureName);
        if isempty(fhandle)
            figure('name',FigureName);
        elseif length(fhandle)>1
            close(fhandle)
            figure('name',FigureName);
        else
            figure(fhandle);
        end
        
        ylims = [-0.1,1.2];
        xlims = [0,GridParams(5)];
        subplot(1,2,1)
        errorbar(bins1,AzAvg1,AzAvgErr1)
        hold on;
        plot(bins1,Comb1DFit1,'r')
        plot(bins1,TF1Fit,'green')
        plot(bins1,GFit1,'black')
        ylim(ylims);
        xlim(xlims)
        hold off;
        subplot(1,2,2)
        errorbar(bins2,AzAvg2,AzAvgErr2)
        hold on;
        plot(bins2,Comb1DFit2,'r')
        plot(bins2,TF2Fit,'green')
        plot(bins1,GFit2,'black')
        ylim(ylims);
        xlim(xlims)
        hold off;
        
        %      subplot(2,2,3)
        %      errorbar(bins1,AzAvg1,AzAvgErr1);
        %      hold on;
        % %      plot(bins1,FhDG1(bins1),'r')
        % %      plot(bins1,gaussian1D(FPDbleG1(1:4),bins1,0),'green')
        % %      plot(bins1,gaussian1D(FPDbleG1(5:end),bins1,0),'black')
        %      ylim(lims);
        %      hold off;
        %      subplot(2,2,4)
        %      errorbar(bins2,AzAvg2,AzAvgErr2);
        %      hold on;
        % %      plot(bins2,FhDG2(bins2),'r')
        % %      plot(bins2,gaussian1D(FPDbleG2(1:4),bins2,0),'green')
        % %      plot(bins2,gaussian1D(FPDbleG2(5:end),bins2,0),'black')
        %      ylim(lims);
        %      hold off;
        
    catch Exception
        MsgString = getReport(Exception);
        AzParams1 = zeros(1,8);
        AzPNames1 = repmat({'P'},1,length(AzParams1));
        
        AzParams2 = zeros(1,8);
        AzPNames2 = repmat({'P'},1,length(AzParams2));
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
    BeamSD1 = std2(ImagesCrop(:,:,3));
    BeamSD2 = std2(ImagesCrop(:,:,5));
    BeamAvg2 = mean(mean(ImagesCrop(:,:,6)));
    DarkAvg = mean(mean(ImagesCrop(:,:,7)));
    CountsAtFitOD1 = (BeamAvg1-DarkAvg)*exp(-AG1);
    CountsAtFitOD2 = (BeamAvg2-DarkAvg)*exp(-AG2);
    
    
    
    NSum1 = getAtomNumber(sum(sum(ImgCrop1)),0,PixelArea,Constants);
    NSum2 = getAtomNumber(sum(sum(ImgCrop2)),0,PixelArea,Constants);
    
    NFit2 = getAtomNumber(abs(2*pi*SxG2*SyG2*AG2),0,PixelArea,Constants);
    NFit1 = getAtomNumber(abs(2*pi*SxG1*SyG1*AG1),0,PixelArea,Constants);
    
    NPix = SizeY*SizeX;
    NSumNoBg1 = getAtomNumber(sum(sum(ImgCrop1))-NPix*Bg1,0,PixelArea,Constants);
    NSumNoBg2 = getAtomNumber(sum(sum(ImgCrop2))-NPix*Bg2,0,PixelArea,Constants);
    
    PixelAreaBinned=Settings.SoftwareBinSizeH*Settings.SoftwareBinSizeV*PixelArea;
    
    NGauss1 = getAtomNumber(2*pi*AzParams1(3)^2*abs(AzParams1(4)),0,PixelAreaBinned,Constants);
    NTotal1 = getAtomNumber(2*pi*AzParams1(3)^2*abs(AzParams1(4))+2*pi*AzParams1(1)^2*abs(AzParams1(2)),0,PixelAreaBinned,Constants);
    CondFrac1 = (NTotal1-NGauss1)/NTotal1;
    
    NGauss2 = getAtomNumber(2*pi*AzParams2(3)^2*abs(AzParams2(4)),0,PixelAreaBinned,Constants);
    NTotal2 = getAtomNumber(2*pi*AzParams2(3)^2*abs(AzParams2(4))+2*pi*AzParams2(1)^2*abs(AzParams2(2)),0,PixelAreaBinned,Constants);
    CondFrac2 = (NTotal2-NGauss2)/NTotal2;
    
    PolAzAvg = -(NTotal1-NTotal2)/(NTotal1+NTotal2);
    
    %     NDbleGaussA1 =FPDbleG1(2)^2*FPDbleG1(3);
    %     NDbleGaussB1 = FPDbleG1(6)^2*FPDbleG1(7);
    %     NTotal1 = NDbleGaussA1+NDbleGaussB1;
    %     if abs(FPDbleG1(2))<abs(FPDbleG1(6))
    %         CondFracDG1 = NDbleGaussA1/NTotal1;
    %     else
    %         CondFracDG1 = NDbleGaussB1/NTotal1;
    %     end
    %
    %     NDbleGaussA2 =FPDbleG2(2)^2*FPDbleG2(3);
    %     NDbleGaussB2 = FPDbleG2(6)^2*FPDbleG2(7);
    %     NTotal2 = NDbleGaussA2+NDbleGaussB2;
    %     if abs(FPDbleG2(2))<abs(FPDbleG2(6))
    %         CondFracDG2 = NDbleGaussA2/NTotal2;
    %     else
    %         CondFracDG2 = NDbleGaussB2/NTotal2;
    %     end
    
    if NFit1 == NFit2
        Polarization = 0;
    else
        Polarization = -((NFit1-NFit2)/(NFit1+NFit2));
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
    BeamSD1,CountsAtFitOD1,BeamAvg2,CountsAtFitOD2,BeamSD2,DarkAvg,...
    AzParams1,NTotal1,AzParams2,NTotal2,CondFrac1,CondFrac2,PolAzAvg];%,CondFracDG1,CondFracDG2];
MoreKeys = cat(2,PNames1,{'NSum1','NSumNoBG1','NFit1','Tx1','Ty1'}...
    ,PNames2,{'NSum2','NSumNoBg2','NFit2','Tx2','Ty2'},...
    {'Polarization','Atoms1 AvgCts','Atoms2 AvgCts',...
    'Beam1 AvgCts','Beam1 SD','Estimated Counts At Fit OD1',...
    'Beam2 AvgCts','Beam2 SD','Estimated Counts At Firt OD2','Dark AvgCts'},...
    AzPNames1,{'NAzAvg1'},AzPNames2,{'NAzAvg2'},{'CondFrac1','CondFrac2','Polarization AzAvg'});%,'CondFracDG1','CondFracDG2'});
end
