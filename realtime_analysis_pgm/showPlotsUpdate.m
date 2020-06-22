function FigHandle = showPlots(ImgDataClass, Settings, FigHandle)
% showPlots(ImgDataClass, Settings, FigHandle)
% Imgs and Fits must both be the same size. They are stacks of images of size [Ny,Nx,NImgs].
% Uses information from the settings file to decide on which plots to show
% in what locations etc. Also converts from 'local' ROI coordinates to 'global' coordinates
% of the original uncropped picture.

% TODO:
% I think it would be much better to get rid of this function and instead
% accept a figure directly to show. 

%make the appropriate figure active
figure(FigHandle)

Imgs = ImgDataClass.ProcessedImages;
Fits = ImgDataClass.Fits;
Vals = ImgDataClass.Vals;
Keys = ImgDataClass.Keys;
Stamp = ImgDataClass.TimeStamp;

ImgMult = Settings.ImageMultiplicity;
ImgLocs = Settings.ImageLocations;
ImgLims = Settings.ImageLimits;
FitMult = Settings.FitMultiplicity;
FitLocs = Settings.FitLocations;
FitLims = Settings.FitLimits;
ResBool = Settings.ResidualBool;
ResLocs = Settings.ResidualLocations;
ResLims = Settings.ResidualLimits;

%Handle markers
MarkerLocations = Settings.MarkerLocations;
if ~isempty(MarkerLocations)
    XROI = Settings.WindowLimits(1:2);
    YROI = Settings.WindowLimits(3:4);
    x_roi_start = XROI(1);
    y_roi_start = YROI(1);
    nx_bin = Settings.SoftwareBinSizeH;
    ny_bin = Settings.SoftwareBinSizeV;
    transform_params = [x_roi_start, y_roi_start, nx_bin, ny_bin];
    [XMarkerPts, YMarkerPts] = img2roi_coord(transform_params, MarkerLocations(:, 1), MarkerLocations(:, 2), 1);
    MarkerSize = 200;
else
    XMarkerPts = 1;
    YMarkerPts = 1;
    MarkerSize = 0;
end
%Get rid of warnings about setting marker size to zero.
warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode')

%Find number of desired plots.
PlotsIndices = max([ImgLocs; FitLocs; ResLocs], [], 1);
NPlotsY = PlotsIndices(1);
NPlotsX = PlotsIndices(2);
PlotIndices_To_SingleIndex = @(Y,X) (Y-1) * NPlotsX + X ;

%Generate image indices with multiplicity.
%Need to handle case where different numbers of Fits and Images
ImageIndices = [];
FitIndices = [];
for ii = 1:length(ImgMult)
    ImageIndices = cat(2, ImageIndices, repmat(ii, 1, ImgMult(ii)));
    FitIndices = cat(2, FitIndices, repmat(ii, 1, FitMult(ii)));
end
ResIndices = (1:length(ResBool)).*ResBool;
ResIndices = ResIndices(ResIndices~=0);

%display images
colormap(Settings.ColorMap);
if (sum(ImgMult) == 1) &&...
        (sum(FitMult) == 0) &&...
        (sum(ResBool) == 0) &&...
        (PlotIndices_To_SingleIndex(ImgLocs(1, 1), ImgLocs(1, 2)) == 1)
    
    %if there is only a single image, display full screen
    ImageIndex = ImageIndices(ImgMult(ImgMult ~= 0));
    imagesc(Imgs(:, :, ImageIndex), ImgLims(1, :))
    colorbar
    axis equal;
    axis image;
    hold on;
    ax = gca;
    localToGlobalTicMarks(Settings, ax);
    if XMarkerPts >=1 &&...
            XMarkerPts <= size(Imgs(:,:,ImageIndex), 2) &&...
            YMarkerPts >=1 &&...
            YMarkerPts <= size(Imgs(:,:,ImageIndex), 1)
        
        scatter(XMarkerPts,YMarkerPts,10,MarkerSize,'r','x')
    else
        disp('Scatter point is outside of picture');
    end
    hold off;
else
    for jj=1:sum(ImgMult)
        I = PlotIndices_To_SingleIndex(ImgLocs(jj, 1), ImgLocs(jj, 2));
        ax = subplot(NPlotsY, NPlotsX, I);
        ImageIndex = ImageIndices(jj);
        imagesc( Imgs(:, :, ImageIndex), ImgLims(jj, :) );
        colorbar;
        axis equal;
        axis image;
        hold on;
        localToGlobalTicMarks(Settings, ax);
        if XMarkerPts >=1 &&...
                XMarkerPts <= size(Imgs(:, :, ImageIndex), 2) &&...
                YMarkerPts >=1 &&...
                YMarkerPts <= size(Imgs(:, :, ImageIndex), 1)
            
            scatter(XMarkerPts, YMarkerPts, 10, MarkerSize, 'r', 'x')
        else
            disp('Scatter point is outside of picture');
        end
        hold off;
    end
end

%display plots
for kk = 1:sum(FitMult)
    I = PlotIndices_To_SingleIndex(FitLocs(kk, 1), FitLocs(kk, 2));
    ax = subplot(NPlotsY, NPlotsX, I);
    FitIndex = FitIndices(kk);
    imagesc(Fits(: ,:, FitIndex), FitLims(kk, :));
    colorbar;
    axis equal;
    axis image;
    hold on;
    localToGlobalTicMarks(Settings, ax);
    if XMarkerPts >=1 &&...
            XMarkerPts <= size(Fits(:, :, FitIndex), 2) &&...
            YMarkerPts >=1 &&...
            YMarkerPts <= size(Fits(:, :, FitIndex), 1)
        
        scatter(XMarkerPts, YMarkerPts, 10, MarkerSize, 'r', 'x')
    else
        disp('Scatter point is outside of picture');
    end
    hold off;
end

%display residuals
for ll = 1:length(ResIndices)
    I = PlotIndices_To_SingleIndex(ResLocs(ll, 1), ResLocs(ll, 2));
    ax = subplot(NPlotsY, NPlotsX, I);
    ResIndex = ResIndices(ll);
    imagesc(Imgs(:, :, ResIndex) - Fits(:, :, ResIndex), ResLims(ll, :));
    colorbar;
    axis equal;
    axis image;
    hold on;
    localToGlobalTicMarks(Settings, ax);
    if XMarkerPts >=1 &&...
            XMarkerPts <= size(Imgs(:, :, ResIndex), 2) &&...
            YMarkerPts >=1 &&...
            YMarkerPts <= size(Imgs(:, :, ResIndex), 1)
        
        scatter(XMarkerPts, YMarkerPts, 10, MarkerSize, 'r', 'x')
    else
        disp('Scatter point is outside of picture');
    end
    hold off;
end

%Change data tip so that it displays coordinates for the full image,
%instead of coordinates for the cropped image. This way, if you change your
%region of interest (ROI) the same coordinates will still point to the same
%location.
dcm_obj = datacursormode(gcf);
updaterFun = @(obj, event_obj) update_data_tip_fn(Settings, obj, event_obj);
set(dcm_obj, 'UpdateFcn', updaterFun);

%Create title for all plots.
suptitle(char(Stamp));

warning('on', 'MATLAB:handle_graphics:exceptions:SceneNode');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function localToGlobalTicMarks(Settings, hAxes)
% Converts the tic marks of a given figure from ROI coordinates
% to global coordinates (with respect to the original image).

XROI = Settings.WindowLimits(1:2);
YROI = Settings.WindowLimits(3:4);
x_local = hAxes.XTick;
y_local = hAxes.YTick;

xform_params = [XROI(1), YROI(1), Settings.SoftwareBinSizeH, Settings.SoftwareBinSizeV];
[x_global, ~] = roi2img_coord(xform_params, x_local, zeros(size(x_local)));
[~, y_global] = roi2img_coord(xform_params, zeros(size(y_local)), y_local);

hAxes.XTickLabels = num2cell(x_global);
hAxes.YTickLabels = num2cell(y_global);
end
