function [Img] = displayFile(FName)
%[Img] = displayFile(FName)
%displays all images in file.

Img = readimg(FName);
Nimgs = size(Img,3);

if Nimgs == 3
    %OD Files
    OD = getOD(Img);
    
    subplot(2,2,1)
    imagesc(Img(:,:,1),[0,3500])
    colorbar
    axis equal
    
    subplot(2,2,2)
    imagesc(Img(:,:,2),[0,3500])
    colorbar
    axis equal
    
    subplot(2,2,3)
    imagesc(Img(:,:,3),[0,1000])
    colorbar
    axis equal
    
    subplot(2,2,4)
    imagesc(OD,[2.3,3.8])
    colorbar
    axis equal
    
else
    %Other files
    Nplots = ceil(Nimgs/2)*2;
    subplot(2,Nplots/2,1)
    imagesc(Img(:,:,1))
    [Cmin,Cmax] = caxis(gca);
    for ii = 2:Nimgs
        subplot(2,Nplots/2,ii);
        imagesc(Img(:,:,ii),[Cmin,Cmax]);
    end
    
end


end

