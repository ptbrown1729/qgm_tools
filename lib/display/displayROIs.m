function displayROIs(Images,Settings)
%displayROIs(Images,Settings)
%Takes stack of 7 images...displays regions of interest used for
%de-stripping and normalizing od images.

    figure('Name','Regions of Interest')
    %display whole area
    OD = getOD(Images(:,:,[2,5,7]));
    imagesc(OD,[Settings.ImageLimits(1,1),Settings.ImageLimits(1,2)])
    colorbar
    axis equal;
    axis image;
    hold on;
    %main ROI
    X1 = Settings.WindowLimits(1);
    X2 = Settings.WindowLimits(2);
    Y1 = Settings.WindowLimits(3);
    Y2 = Settings.WindowLimits(4);
    LegROI = plot([X1,X2],[Y1,Y1],'black','LineWidth',2);
    plot([X1,X2],[Y2,Y2],'black','LineWidth',2);
    plot([X1,X1],[Y1,Y2],'black','LineWidth',2);
    plot([X2,X2],[Y1,Y2],'black','LineWidth',2);
    
    %ROI used for stripe removal.
    StrX1 = Settings.StripesRegionLimits(1);
    StrX2 = Settings.StripesRegionLimits(2);
    StrY1 = Settings.StripesRegionLimits(3);
    StrY2 = Settings.StripesRegionLimits(4);
    LegStrROI = plot([StrX1,StrX2],[StrY1,StrY1],'b','LineWidth',2);
    plot([StrX1,StrX2],[StrY2,StrY2],'b','LineWidth',2);
    plot([StrX1,StrX1],[StrY1,StrY2],'b','LineWidth',2);
    plot([StrX2,StrX2],[StrY1,StrY2],'b','LineWidth',2);
    
    %ROI used for background removal.
    X1Norm = Settings.NormRegionLimits(1);
    X2Norm = Settings.NormRegionLimits(2);
    Y1Norm = Settings.NormRegionLimits(3);
    Y2Norm = Settings.NormRegionLimits(4);
    LegNormROI = plot([X1Norm,X2Norm],[Y1Norm,Y1Norm],'r','LineWidth',2);
    plot([X1Norm,X2Norm],[Y2Norm,Y2Norm],'r','LineWidth',2);
    plot([X1Norm,X1Norm],[Y1Norm,Y2Norm],'r','LineWidth',2);
    plot([X2Norm,X2Norm],[Y1Norm,Y2Norm],'r','LineWidth',2);
    
    legend([LegROI,LegStrROI,LegNormROI],{'ROI','Stripes ROI','Norm ROI'})


end

