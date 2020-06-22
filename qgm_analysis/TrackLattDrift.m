Saving = 0;

% DateCell = {2017 08 20};
% Folders = 15:37;

% DateCell = {2017 08 21};
% Folders = 7:41;

% DateCell = {2017 08 22};
% Folders = 6:53;

% DateCell = {2017 08 23};
% Folders = 7:49;

% DateCell = {2017 08 24};
% Folders = 7:53;

DateCell = {2019 04 23};
Folders = 1:100;

Cxs_CoM = [];
Cys_CoM = [];
times = [];
phi1s = [];
phi2s = [];

for ii = 1:length(Folders)
    dset = horzcat(DateCell, Folders(ii), 1, 1);
    df = DataFolder(dset, [], [0, 1]);
%     df.Dataset = dset;
    if ~isnan(df.X_CoM_fullpic)
        CropSize = [df.ImgCropSize, df.ImgCropSize];
        [ImgStack, XStarts_ROI, YStarts_ROI, ~, x_coms, y_coms, ~, ~, ~, FlFPaths] = df.loadPics(df.Dataset, CropSize, df.CenterStyle);
        [~, phi1s_curr, phi2s_curr, sort_index] = df.getRecPhases();
        times_curr = df.getPicTimes(df.FlFPaths);

        % aggregate data
        Cxs_CoM = cat(1, Cxs_CoM, df.X_CoM_fullpic);
        Cys_CoM = cat(1, Cys_CoM, df.Y_CoM_fullpic);
        times = cat(1, times, times_curr);
        phi1s = cat(1, phi1s, phi1s_curr);
        phi2s = cat(1, phi2s, phi2s_curr);
    end
end

[times, I] = sort(times);
Cxs_CoM = Cxs_CoM(I);
Cys_CoM = Cys_CoM(I);
phi1s = phi1s(I);
phi2s = phi2s(I);

NanPos = isnan(Cxs_CoM)|isnan(Cys_CoM);
times = times(~NanPos);
Cxs_CoM = Cxs_CoM(~NanPos);
Cys_CoM = Cys_CoM(~NanPos);
phi1s = phi1s(~NanPos);
phi2s = phi2s(~NanPos);

% unwrap all phases
phi1s = unwrap(2 * pi * phi1s) / (2*pi);
phi2s = unwrap(2 * pi * phi2s) / (2*pi);

% fftx = fftshift(fft(Cxs_CoM-mean(Cxs_CoM)));
% ffty = fftshift(fft(Cys_CoM-mean(Cys_CoM)));

dat = struct();
dat.date = DateCell;
dat.times = times;
dat.x_com = Cxs_CoM;
dat.y_com = Cys_CoM;
dat.phi1s = phi1s;
dat.phi2s = phi2s;



%%
%%%Plot
smooth_size = 1;

x_tic_times = linspace(min(times), max(times), 6); %Times(1:floor(length(Times)/10):length(Times));
XTicLabels = {};
for ii = 1:length(x_tic_times)
    XTicLabels{ii} = datestr(x_tic_times(ii), 'HH:MM');
end

% plot results
figure;
nrows = 2;
ncols = 2;

subplot(nrows, ncols, 1)
plot(times, smooth(Cxs_CoM, smooth_size), 'b.');
grid on;
title(sprintf('Cx = %0.2f +/- %0.2f', mean(Cxs_CoM), std(Cxs_CoM)));

ax = gca;
ax.XTick = x_tic_times;
ax.XTickLabel = XTicLabels;

subplot(nrows, ncols, 2)
plot(times, smooth(Cys_CoM,smooth_size), 'r.');
grid on;
title(sprintf('Cy = %0.2f +/- %0.2f', mean(Cys_CoM), std(Cys_CoM)));
ax = gca;
ax.XTick = x_tic_times;
ax.XTickLabel = XTicLabels;

% phi1s
subplot(nrows, ncols, 3)
plot(times, phi1s, 'b.');
grid on;
title('phi1s');
ax = gca;
ax.XTick = x_tic_times;
ax.XTickLabel = XTicLabels;

% phi2s
subplot(nrows, ncols, 4)
plot(times, phi2s, 'b.');
grid on;
title('phi2s');
ax = gca;
ax.XTick = x_tic_times;
ax.XTickLabel = XTicLabels;

suptitle(sprintf('Lattice Position Vs. Time, %04d-%02d-%02d', DateCell{1}, DateCell{2}, DateCell{3}));

if Saving   
    fname = sprintf('%04d-%02d-%02d_LattDrift_Vs_Time.mat', DateCell{1}, DateCell{2}, DateCell{3});
    save(fname, 'dat', '-struct');
    
    fname = sprintf('%04d-%02d-%02d_LattDrift_Vs_Time.txt', DateCell{1}, DateCell{2}, DateCell{3});
    title_cell = {'time (seconds)', 'x com (sites)', 'y com (sites)', 'phi1 (sites)', 'phi2 (sites)'};
    time_seconds = (dat.times - dat.times(1)) * 24 * 60 * 60;
    
    data_array = horzcat(time_seconds, dat.x_com, dat.y_com, dat.phi1s, dat.phi2s);
    save_data_file(data_array, title_cell, '\t', '', '', fname);
end

%%
DriftFiles = dir('*LattDrift_Vs_Time.mat');

times = [];
Cx = [];
Cy = [];

for ii = 1:length(DriftFiles)
    disp(ii);
    a = load(DriftFiles(ii).name);
    times = cat(1,times,a.Times);
    Cx = cat(1,Cx,a.Cxs_CoM);
    Cy = cat(1,Cy,a.Cys_CoM);
end

%%%Sort
[times,I] = sort(times);
Cx = Cx(I);
Cy = Cy(I);

% fftx = fftshift(fft(Cxs_CoM-mean(Cxs_CoM)));
% ffty = fftshift(fft(Cys_CoM-mean(Cys_CoM)));

x_tic_times = linspace(min(times),max(times),6); %Times(1:floor(length(Times)/10):length(Times));
XTicLabels = {};
for ii = 1:length(x_tic_times)
    XTicLabels{ii} = datestr(x_tic_times(ii),'mm/DD-HH:MM');
end

%%%Plot
smooth_size = 10;

figure;
subplot(1,2,1)
plot(times,smooth(Cx,smooth_size),'b.');
grid on;
title(sprintf('Cx = %0.2f +/- %0.2f',mean(Cx),std(Cx)));

ax = gca;
ax.XTick = x_tic_times;
ax.XTickLabel = XTicLabels;

subplot(1,2,2)
plot(times,smooth(Cy,smooth_size),'r.');
grid on;
title(sprintf('Cy = %0.2f +/- %0.2f',mean(Cy),std(Cy)));
ax = gca;
ax.XTick = x_tic_times;
ax.XTickLabel = XTicLabels;

suptitle('Lattice Position Vs. Time');