df = DataFolder();

[xx, yy] = meshgrid(1:100, 1:100);
theta = 10 * pi/180;
period = 10;
xx_rot = xx * cos(theta) - yy * sin(theta);
yy_rot = xx * sin(theta) + yy * cos(theta);
stripes = cos(2*pi*yy_rot/period);
stripes(stripes > 0) = 1;
stripes(stripes < 0) = 0;

n_stack = zeros(100, 100, 2);
n_stack(:, :, 1) = stripes;
n_stack(:, :, 2) = ~stripes;


AzAvgGrid = zeros(100, 100);
BinEndPts = [-1, 1];
NumNeighbors = 10;
RestrictSameBin = 1;

[nn, nnsdm, ni, nisdm, nj, njsdm, npts] = ...
                    df.getCorrMoments(n_stack, AzAvgGrid, BinEndPts, NumNeighbors, RestrictSameBin);
[nnc, nncunc] = df.getCorrWithErr(nn, ni, nj, npts);

df.getColorMap();
figure;

subplot(2, 2, 1)
imagesc(n_stack(:, :, 1));
axis equal; 
axis image;
title('1st image');

subplot(2, 2, 2)
imagesc(n_stack(:, :, 2));
axis equal; 
axis image;
title('2nd image');

subplot(2, 2, 4);
imagesc(squeeze(nnc), [-0.25, 0.25]);
axis equal;
axis image;
colormap(df.ColorMap);
colorbar;