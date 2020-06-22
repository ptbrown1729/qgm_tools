% RecFilePath = 'Z:\IMAGING DATA\2017\03\08\026_CDW_WithReservoir_From150-250_DMDPix\2017_03_08_026_001_002_reconstr.mat';
% ImgFilePath = 'Z:\IMAGING DATA\2017\03\08\026_CDW_WithReservoir_From150-250_DMDPix\17_21_51_fl.fits';

ImgFilePath = '\\128.112.86.75\lithium\IMAGING DATA\2017\03\13\016_Doublon-Doublon-Corr\16_51_41_fl.fits';
RecFilePath = '\\128.112.86.75\lithium\IMAGING DATA\2017\03\13\016_Doublon-Doublon-Corr\2017_03_13_016_001_002_reconstr.mat';

border_pix = 5;

Img = readimg(ImgFilePath);
Img = Img(border_pix:end-border_pix, border_pix:end-border_pix, 2);


Rec = load(RecFilePath);
RecImg = Rec.occupationsRounded;
gp = Rec.coordTrafo.gridParameters;
IndexOffset = Rec.coordTrafo.indexOffset;

%center of reconstruction is roughly
%(188,157) -> (26.5,-4.5) %LHS is index, RHS is relative to center
%center of picture is roughly
%(450,471) -> (49,70)

Theta1 = gp.theta1;
Theta2 = gp.theta2;
Lambda1 = gp.lambda(1);
Lambda2 = gp.lambda(2);
Phi1 = gp.phi1;
Phi2 = gp.phi2;
Offset1 = IndexOffset(1);
Offset2 = IndexOffset(2);

%The confusing about these functions is...
%the image coordinates are ordered as a function i.e. (x,y)
%where as the reconstruction coords are ordered as a matrix i.e. (y,x)
gp.latticeToImage([172,149], IndexOffset)
% gp.latticeToImageSimple([152,188],IndexOffset)

gp.imageToLattice([330,397], IndexOffset) %=[m,n] which is in matrix index order.


x = 1:size(RecImg, 2);
y = 1:size(RecImg, 1);
[Xrec,Yrec] = meshgrid(x,y);




[Xlatt,Ylatt] = meshgrid(1:size(RecImg,2),1:size(RecImg,1));
[Ximg,Yimg] = meshgrid(1:size(Img,2),1:size(Img,1));

%STILL NOT QUITE WORKING
RecInterp = @(Xl,Yl) interp2(Xlatt,Ylatt,RecImg,Xl,Yl);
% RecImgSp = @(Xi,Yi) RecInterp((- Xi*sin(Theta1) + Yi*cos(Theta1))*1/Lambda2*1/sin(Theta2-Theta1) + Phi2 - Offset2,(Xi*sin(Theta2) - Yi*cos(Theta2))*1/Lambda1*1/sin(Theta2-Theta1) + Phi1 - Offset1);
RecImgSp = @(Xi,Yi) RecInterp((-Xi*cos(Theta1) + Yi*sin(Theta1))*1/Lambda2*1/sin(Theta1-Theta2) + Phi2 - Offset2,(Xi*cos(Theta2) - Yi*sin(Theta2))*1/Lambda1*1/sin(Theta1-Theta2) + Phi1 - Offset1);

ImgInterp = @(Xi,Yi) interp2(Ximg,Yimg,Img,Xi,Yi);
ImgLattSp = @(Xl,Yl) ImgInterp((Xl+Offset2-Phi2)*Lambda2*sin(Theta2) + (Yl+Offset1-Phi1)*Lambda1*sin(Theta1),(Xl+Offset2-Phi2)*Lambda2*cos(Theta2) + (Yl+Offset1-Phi1)*Lambda1*cos(Theta1));




figure('name','Rec Transform')
NRows = 3;
NCols = 2;
subplot(NRows,NCols,1)
    imagesc(Img,[144,250]);
    axis equal; axis image;
    title('Raw Image')
subplot(NRows,NCols,2)
    imagesc(RecImg,[0,1]);
    axis equal; axis image;
    title('Reconstruction')
subplot(NRows,NCols,3)
    imagesc(RecImgSp(Ximg,Yimg),[0,1]);
    axis equal; axis image;
    title('Rec, Img Space')
subplot(NRows,NCols,4)
    imagesc(ImgLattSp(Xlatt,Ylatt),[144,250]);
    axis equal; axis image;
    title('Image, Latt Space')
subplot(NRows,NCols,5)
    imagesc(RecImgSp(Ximg,Yimg)-(Img-144)/250,[-1,1])
    axis equal; axis image;
    title('Rec-Img, Img Space')
subplot(NRows,NCols,6)
    imagesc(RecImg-(ImgLattSp(Xlatt,Ylatt)-144)/250,[-1,1]);
    axis equal; axis image;
    title('Rec-Img, Latt Space')