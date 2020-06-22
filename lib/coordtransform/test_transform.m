theta1 = -0.964508807972182;
theta2 = 0.605806785333545;
phi1 = 0.307810141410327;
phi2 = 0.669335387203412;
lambda1 = 3.521786000000000;
lambda2 = 3.424710000000000;
offset1 = -184;
offset2 = 0;
latt_xform_params = [theta1, theta2, phi1, phi2, lambda1, lambda2, offset1, offset2];

fig_handle = figure;

% define square in reconstruction space
% CCW from lower left
xs_latt = [30, 40, 40, 30];
ys_latt = [30, 30, 40, 40];

subplot(2, 2, 1);
scatter(xs_latt, ys_latt);
title('latt space (atoms/reconstruction)');

% get parallelogram in image space
% [xs_img, ys_img] = latticeToImg(xs_latt, ys_latt, latt_xform_params);
[xs_img, ys_img, ~] = latt2img_coord(latt_xform_params, xs_latt, ys_latt);

subplot(2, 2, 2);
scatter(xs_img, ys_img);
title('img space (atoms/fl)');


% get parallelogram in object space
cx_img = 500;
cy_img = -203;
scale = 0.5;
angle = 0;
obj_xform_params = [cx_img, cy_img, angle, scale];

[xs_obj, ys_obj] = img2obj_coord(obj_xform_params, xs_img, ys_img);

subplot(2, 2, 3);
scatter(xs_obj, ys_obj);
title('obj space (DMD)');

%% determine transforms...

% DMD space [mirror_x1, mirror_x2] x [mirror_y1, mirror_y2] = [500, 550] x [360, 410]
% DMD coordinates [499.5, 550.5] x [359.5, 410.5]
% Using center of first pixel in upper left corner as (x, y) = (1, 1)

folders = 95:99;
xstarts_dmd = [500, 400, 500, 550, 600];
xends_dmd = [550, 450, 550, 600, 650];
xcom_dmd = 0.5 * (xstarts_dmd + xends_dmd);
ystarts_dmd =[360, 360, 460, 310, 380];
yends_dmd = [410, 410, 510, 310, 430];
ycom_dmd = 0.5 * (ystarts_dmd + yends_dmd);

x_guesses = [47,55,63,35,52];
y_guesses = [49,59,47,46,35];
len = 8;

% fit square to find center of square
mask_fn = @(x, y) heaviside(x + 0.5) .* heaviside(0.5 - x) .* heaviside(y + 0.5) .* heaviside( 0.5 - y);
fit_fn = @(P, X, Y) mask2d_2(P, X, Y, mask_fn);

xcom_latt = zeros(size(xcom_dmd));
ycom_latt = zeros(size(ycom_dmd));
for ii = 1:length(folders)
    df = DataFolder({2018 1 23 folders(ii) 1 1}, []);
    
    % reconstruction image is cropped
    x_crop_start = df.XStarts_ROI(1);
    y_crop_start = df.YStarts_ROI(1);
    
    img_crop = df.Occs_ImgAvg;
    x_crop_guess = x_guesses(ii);
    y_crop_guess = y_guesses(ii);
    img_crop_roi = img_crop(y_crop_guess - len : y_crop_guess + len, x_crop_guess - len : x_crop_guess + len);
    x_roi_start = x_crop_guess - len;
    y_roi_start = y_crop_guess - len;
    
    [xx_roi, yy_roi] = meshgrid(1:size(img_crop_roi, 2), 1:size(img_crop_roi, 1));
    roi_xform = [x_crop_start + x_roi_start - 1, y_crop_start + y_roi_start - 1, 1, 1];
    [xx_img, yy_img] = roi2img_coord(roi_xform, xx_roi, yy_roi);
    
    % threshold and center of mass
    
    cutoffs = 0.1:0.1:0.8;
    coms_x = zeros(1, length(cutoffs));
    coms_y = zeros(1, length(cutoffs));
    for ii = 1:length(cutoffs)
        img_mod = 1 - img_crop_roi;
        img_mod(img_mod < cutoffs(ii)) = 0;
        [com_x, com_y, w] = get_moment(img_mod, 1);
        coms_x(ii) = com_x;
        coms_y(ii) = com_y;
    end
    
    [xcom_latt(ii), ycom_latt(ii)] = roi2img_coord(roi_xform, coms_x(end), coms_y(end));
    
    % fit square
    init_params = [10, 9, 0, 6];
    fixed_params = [0, 0, 1, 0];
    lb = [-inf, -inf, -inf, -inf];
    ub = [inf, inf, inf ,inf];
    % this is not converging. One idea is that the step size is too small.
    % Here you have to change by some reasonable fraction of an integer to
    % see a difference.
    [fitp, ~, ffh, std_err] = fit2D(xx_roi, yy_roi, img_mod, [], fit_fn,...
        init_params, fixed_params, lb, ub, 4, 'lsqnonlin');
    
    figure;
    subplot(2,2,1);
    plot(cutoffs, coms_x,'o');
    grid on;
    xlabel('cutoff density');
    ylabel('x com');
    
    subplot(2,2,2);
    plot(cutoffs, coms_y,'o');
    grid on;
    xlabel('cutoff density');
    ylabel('y com');
    
    subplot(2,2,3);
    %     imagesc(img_crop_roi);
    imagesc(img_mod);
    axis equal;
    axis image;
    hold on;
    scatter(com_x, com_y, 'rx');
    
    subplot(2,2,4);
    imagesc(ffh(xx_roi, yy_roi));
    axis equal;
    axis image;
    hold on;
    % also can try fitting
    
end

%%
% new data set 01/25/2018
folders = [58, 59];
dfs = [];

for ii = 1:length(folders)
    df = DataFolder();
    DistGrid = [];
    BinEdges = [];
    if ii ~= 1
        df.AzAvgType = 'external';
        df.CenterStyle = 'Fixed';
        df.CroppedPicStartCoords = dfs(1).CroppedPicStartCoords;
        df.Cx_AzAvg = dfs(1).Cx_AzAvg;
        df.Cy_AzAvg = dfs(1).Cy_AzAvg;
        df.AzAvgCentering = 'external';
        DistGrid = dfs(1).DistGrid;
        BinEdges = dfs(1).BinEdges;
    end
    
    df.initialize({2018 01 25 folders(ii) 1 1}, [],...
        BinEdges, df.AzAvgType, DistGrid);
    dfs = vertcat(dfs, df);
end

img = dfs(2).Occs_ImgAvg - dfs(1).Occs_ImgAvg;
figure;

subplot(1, 2, 1);
imagesc(img);
axis equal; 
axis image;

subplot(1, 2, 2);
m = img;
m(m > -0.5) = 0;
imagesc(m);
axis equal;
axis image;
