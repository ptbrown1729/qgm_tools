function [ImageData] = gaussianExpt(ImageDataIn, Settings, Constants, fig_handle)
%[ImageData] = gaussianExpt(ImageDataIn,Settings,Constants,FHandle)
%%%Arguments%%%
%%%ImageDataIn is a ImageDataClass instance which contains the image files
%and other associated data.
%Settings is a SettingsClass instance which contains information necessary
%to process the image data. For example, the region of interest
%coordinates, initial fitting parameters, camera quantum efficiency, etc.
%%%Constants is a ConstantsClass instance which contains physical constants
%such as hbar and the saturation intensity of the lithium-6 D2 line.
%%%FigHandleAlt is a figure handle. Any extra information that needs to be
%displayed graphically should be displayed on this figure. Note that basic
%fitting information is automatically displayed in another way (i.e. ProgramClass()
%takes care of this directly based on the ImageData output.
%%%Output%%%
%ImageData is an ImageDataClass instance which contains image files and
%additional information after processing these files. The most important
%information added is stored in ImageData.Vals and ImageData.Keys, which
%are ultimately written to a log file by ProgramClass()

ImageData = ImageDataIn;

%get information from Settings and ImageDataIn instances.
imgs = ImageDataIn.OriginalImages;
vals = ImageDataIn.Vals;
keys = ImageDataIn.Keys;

%set up ROI
if isempty(Settings.WindowLimits)
    X1 = 1;
    X2 = size(imgs, 2);
    Y1 = 1;
    Y2 = size(imgs, 1);
else
    %make sure these don't exceed the dimensions of the image.
    %These are in 'local' coordinates.
    X1 = Settings.WindowLimits(1);
    X2 = min([Settings.WindowLimits(2), size(imgs,2)]);
    Y1 = Settings.WindowLimits(3);
    Y2 = min([Settings.WindowLimits(4), size(imgs,1)]);
end

%Crop image based on settings file values
imgs_crop = imgs(Y1:Y2, X1:X2, :);

%Produce ROI mask
roi_mask = zeros(size(imgs(:, :, 1)));
roi_mask(Y1:Y2, X1:X2) = ones(Y2-Y1+1, X2-X1+1);

%handle both the case where we take a blank picture first, and where we
%only take the three images necessary to calculate optical density.
if(size(imgs_crop, 3) == 3)
    indices = [1, 2, 3];
elseif (size(imgs_crop, 3) == 4)
    indices = [2, 3, 4];
else
    indices = [1, 2, 3];
    error('Unsupported image size in gaussianExpt. Image had %d dimensions. Only accepts 3 or 4 images. Using first three images.',...
        ndims(imgs_crop));
end

% if we are having problems, display more information
if Settings.DebugMode
    displayROIs(imgs, Settings);
end

% remove camera stripes
if Settings.RemoveStripes
    tic;
    
    StrX1 = Settings.StripesRegionLimits(1);
    StrX2 = Settings.StripesRegionLimits(2);
    StrY1 = Settings.StripesRegionLimits(3);
    StrY2 = Settings.StripesRegionLimits(4);
    StripesNormRegion = imgs(StrY1 : StrY2, StrX1 : StrX2, :);
    imgs_crop = normalizeStripes(imgs_crop, StripesNormRegion, Settings.StripesOrientation);
    
    t_remove_stripes = toc;
    fprintf('Removing Stripes took %d \n', t_remove_stripes);
end

% remove OD offsets
if Settings.NormalizeOD
    tic;
    
    X1Norm = Settings.NormRegionLimits(1);
    X2Norm = Settings.NormRegionLimits(2);
    Y1Norm = Settings.NormRegionLimits(3);
    Y2Norm = Settings.NormRegionLimits(4);
    imgs_proc = getNormOD( imgs_crop(:, :, indices),...
        imgs(Y1Norm:Y2Norm, X1Norm:X2Norm, indices) );
    
    TNormOD = toc;
    fprintf('Normalizing OD took %d \n', TNormOD);
else
    imgs_proc = imgs_crop(:, :, indices);
end

% software binning
imgs_proc = binImg(imgs_proc,...
    Settings.SoftwareBinSizeH, Settings.SoftwareBinSizeV, 'Normal');

% calculate OD
tic;
imgs_od = getOD(imgs_proc);
TOD = toc;
fprintf('Getting OD took %d \n', TOD);

%
xs = X1 : X2;
ys = Y1 : Y2;
[xxs, yys] = meshgrid(xs, ys);
ImageData.XcoordsOriginal = xxs; 
ImageData.YcoordsOriginal = yys;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
if Settings.Fitting
    
    functions = {'gaussian2D'};
    init_params = Settings.FitInitialParameters;
    fixed_params = Settings.FitFixedParameters;
    ubs = Settings.FitUpperLims;
    lbs = Settings.FitLowerLims;
    
    if Settings.UseFitAutoGuess
        % TODO: don't think this works... remove?
        [CxAuto, CyAuto, SxAuto, SyAuto] = guessCenter(imgs_od);
        AutoGuesses = [CxAuto, CyAuto, SxAuto, SyAuto];
        for ii = 1:4
            if ~fixed_params(ii)
                init_params(ii) = AutoGuesses(ii);
            end
        end
    end
     
    %Fitting
    transform_params = [X1, Y1, Settings.SoftwareBinSizeH, Settings.SoftwareBinSizeV];
    try
        [xs_local, ys_local] = meshgrid(1:size(imgs_od, 2), 1:size(imgs_od, 1));
        [xxs, yys] = roi2img_coord(transform_params, xs_local, ys_local);

        ImageData.XcoordsProc = xxs; 
        ImageData.YcoordsProc = yys;
        
        [fit_params, pnames, fn_handle] = fit2D(xxs, yys, imgs_od,...
            [], functions, init_params, fixed_params, lbs, ubs);
        imgs_fit = fn_handle(xxs, yys);
        
    catch Exception
        MsgString = getReport(Exception);
        imgs_fit = zeros(size(imgs_od));
        fit_params = zeros(1, 7); %[0,0,0,0,0,0,0];
        pnames = repmat({'P'}, 1, length(fit_params));
        disp('Error performing fit.')
        disp(MsgString);
    end
else
    imgs_fit = zeros(size(imgs_od));
    fit_params = zeros(1, 7); %[0,0,0,0,0,0,0];
    pnames = repmat({'P'}, 1, length(fit_params));
end

%Name variables for later use. These are in global coordinates
CxG = fit_params(1); 
CyG = fit_params(2); 
SxG = fit_params(3); 
SyG = fit_params(4);
AG = fit_params(5); 
ThetaG = fit_params(6); 
Bg = fit_params(7);
%Also keep track of the local coordinates.
[CxGLocal, CyGLocal] = img2roi_coord(transform_params, CxG, CyG, 1);
% CxGLocal = getLocalCoords(CxG, [X1, X2], Settings.SoftwareBinSizeH);
% CyGLocal = getLocalCoords(CyG, [Y1, Y2], Settings.SoftwareBinSizeV);
Tfit = toc;
fprintf('Fitting took %d \n', Tfit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    %Use unprocessed images for beam params.   
    atom_avg = mean(mean( imgs_crop(:, :, indices(1)) ));
    beam_avg = mean(mean( imgs_crop(:, :, indices(2)) ));
    dark_avg = mean(mean( imgs_crop(:, :, indices(3)) ));
    beam_sd = std2( imgs_crop(:, :, indices(2)) );
    
    %also interested in information about the center of the cloud.
    SigmaX_Eff = sqrt( 1 ./ ( cos(ThetaG)^2 / SxG^2 + sin(ThetaG)^2 / SyG^2) );
    WidthX = SigmaX_Eff;
%     WidthX = 1.5 * abs(SxG * cos(ThetaG) + SyG * sin(ThetaG)); %since SxG is the width in the rotated coordinates...
    CenterX1 = max([floor(CxGLocal - WidthX / 2), 1]); %don't want to exceed picture size...
    CenterX2 = min([ceil(CxGLocal + WidthX / 2), size(imgs_crop, 2)]);
    
    SigmaY_Eff = sqrt( 1 ./ ( cos(ThetaG)^2 / SyG^2 + sin(ThetaG)^2 / SxG^2) );
    WidthY = SigmaY_Eff;
%     WidthY = 1.5 * abs(SyG * cos(ThetaG) + SxG * sin(ThetaG));
    CenterY1 = max([floor(CyGLocal - WidthY / 2), 1]);
    CenterY2 = min([ceil(CyGLocal + WidthY / 2), size(imgs_crop,1)]);
    
    AtomAvgCenter = mean(mean( imgs_crop(CenterY1:CenterY2, CenterX1:CenterX2, indices(1)) ));
    BeamAvgCenter = mean(mean( imgs_crop(CenterY1:CenterY2, CenterX1:CenterX2, indices(2)) ));
    DarkAvgCenter = mean(mean( imgs_crop(CenterY1:CenterY2, CenterX1:CenterX2, indices(3)) ));
    BeamSDCenter = std2( imgs_crop(CenterY1:CenterY2, CenterX1:CenterX2, indices(2)) );
    
    fprintf('X1 = %0.2f, X2 = %0.2f, Y1 = %0.2f, Y2 = %0.2f \n',...
        CenterX1, CenterX2, CenterY1, CenterY2);
    
    pixel_area = Settings.getEffPixelArea;
    lambda = Constants.lambda_D2;
    QE = Settings.QuantumEfficiency;
    counts_per_photon = Settings.CountsPerPhoton;
    image_time = Settings.ImagingDuration;
    IoverIsatROI = getIntensity(beam_avg - dark_avg, image_time, lambda, QE, counts_per_photon, pixel_area) / Constants.Isat_D2;
    CountsAtFitOD = (beam_avg - dark_avg) * exp(-AG);
    
    IoverIsatCenter = getIntensity(BeamAvgCenter - DarkAvgCenter,...
        image_time, lambda, QE, counts_per_photon, pixel_area) / Constants.Isat_D2;
%     NFit = getAtomNumber(2 * pi * SxG * SyG * AG, 0, pixel_area, Constants);
    NFit = getAtomNumber(2 * pi * SxG * SyG * AG, Constants.lambda_D2, pixel_area);
    n_sum = getAtomNumber(sum(sum(imgs_od)), Constants.lambda_D2, pixel_area);
    
%     PeakAtomDensity = getAtomNumber(AG, 0, pixel_area, Constants) / pixel_area;
    PeakAtomDensity = getAtomNumber(AG, Constants.lambda_D2, pixel_area);
    
    n_pix = size(imgs_od,1) * size(imgs_od, 2);
    
%     NSumNoBg = getAtomNumber(sum(sum(ImgCrop)) - n_pix * Bg, 0, pixel_area, Constants);
    NSumNoBg = getAtomNumber(sum(sum(imgs_od)) - n_pix * Bg, Constants.lambda_D2, pixel_area);
    
catch Exception
    if(~exist('IoverIsatROI','var'))
        IoverIsatROI = 0;
    end
    
    if(~exist('IoverIsatCenter','var'))
        IoverIsatCenter = 0;
    end
    
    if(~exist('CountsAtFitOD','var'))
        CountsAtFitOD = 0;
    end
    
    if(~exist('NFit','var'))
        NFit = 0;
    end
    
    if(~exist('n_sum','var'))
        n_sum = 0;
    end
    
    if(~exist('NSumNoBg','var'))
        NSumNoBg = 0;
    end
    
    MsgString = getReport(Exception);
    disp(MsgString);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MoreVals = [fit_params, n_sum, NSumNoBg, NFit, atom_avg, beam_avg, dark_avg,...
            beam_sd, CountsAtFitOD, IoverIsatROI, IoverIsatCenter, PeakAtomDensity];
MoreKeys = cat(2,pnames,{'NSum', 'NSumNoBg', 'NFit', 'Atom AvgCts',...
                         'Beam AvgCts', 'Dark AvgCts', 'BeamSD',...
                         'Estimated Counts At Fit OD', 'IOverISat',...
                         'IOverISatCenter', 'PeakAtomDensity'});

%Everything ultimately stored in ImageData
ImageData.ProcessedImages = imgs_od;
ImageData.Fits = imgs_fit;
ImageData.ROIMask = roi_mask;
% TODO: need to modify this so that it doesn't change the order of the
% fitting parameters if their are other vals present, i.e. if cicero is
% connected vs not.
ImageData.Vals = cat(2, vals, MoreVals);
ImageData.Keys = cat(2, keys, MoreKeys);
ImageData.Settings = Settings;

% plot results
% TODO: would like to do this here, instead of separate function...

end

