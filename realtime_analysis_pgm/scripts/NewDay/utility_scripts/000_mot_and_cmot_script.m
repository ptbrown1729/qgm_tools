saving = 1;
folder = 0;
close all
% load log file data
fname = fullfile(getDataSetPath(datenum(2019, 1, 4), folder), 'log.txt');
[vals, keys] = parseLog(fname);

if isempty(keys)
    error();
end

% order of columns in log file can change if server is connected to
% simplicio. Find columns by name.
match_fn = @(c1, c2) strcmp(c1, c2);

peakod_index = find( cellfun(@(c) match_fn(c, 'AG'), keys) );
atomnum_index = find( cellfun(@(c) match_fn(c, 'NFit'), keys) );
cx_index = find( cellfun(@(c) match_fn(c, 'CxG'), keys) );
cy_index = find( cellfun(@(c) match_fn(c, 'CyG'), keys) );
sx_index = find( cellfun(@(c) match_fn(c, 'SxG'), keys) );
sy_index = find( cellfun(@(c) match_fn(c, 'SyG'), keys) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CMOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofs_cmot = [0.2:0.2:3.2]*1e-3;

indices_cmot = 1:11;
indices_cmot = indices_cmot(1: min(size(vals, 1) - min(indices_cmot) + 1, length(indices_cmot)));
peak_ods_cmot = vals(indices_cmot, peakod_index);
% remove any nans
indices_cmot = indices_cmot(~isnan(peak_ods_cmot));
tofs_cmot = tofs_cmot(1:length(peak_ods_cmot));
tofs_cmot = tofs_cmot(~isnan(peak_ods_cmot));
peak_ods_cmot = peak_ods_cmot(~isnan(peak_ods_cmot));
% process other files
atom_num_cmot = vals(indices_cmot, atomnum_index);
cx_cmot = vals(indices_cmot, cx_index);
cy_cmot = vals(indices_cmot, cy_index);
sigmax_cmot = abs(vals(indices_cmot, sx_index));
sigmay_cmot = abs(vals(indices_cmot, sy_index));


min_t_cmot = -inf;
max_t_cmot = inf;
atom_num_cmot = atom_num_cmot(tofs_cmot < max_t_cmot & tofs_cmot > min_t_cmot);
peak_ods_cmot = peak_ods_cmot(tofs_cmot < max_t_cmot & tofs_cmot > min_t_cmot);
cx_cmot = cx_cmot(tofs_cmot < max_t_cmot & tofs_cmot > min_t_cmot);
cy_cmot = cy_cmot(tofs_cmot < max_t_cmot & tofs_cmot > min_t_cmot);
sigmax_cmot = sigmax_cmot(tofs_cmot < max_t_cmot & tofs_cmot > min_t_cmot);
sigmay_cmot = sigmay_cmot(tofs_cmot < max_t_cmot & tofs_cmot > min_t_cmot);
tofs_cmot = tofs_cmot(tofs_cmot < max_t_cmot & tofs_cmot > min_t_cmot);

kb = 1.3806e-23;
m = 9.9883e-27;
pix_size = 4.65e-6;
mag = 30/250;
sigma_fn = @(P, t) sqrt(P(1)^2 + kb / m * P(2) * t.^2 * (mag/pix_size)^2 );

initpx_cmot = [20, 250e-6];
fixedpx_cmot = [0, 0];
lb_cmot = [0, 0];
ub_cmot = [inf, inf];
[fitp_x_cmot, ~, ffh_x_cmot, std_err_x_cmot, ~] = fit1D(tofs_cmot, sigmax_cmot, [], sigma_fn, initpx_cmot, fixedpx_cmot, lb_cmot, ub_cmot);

initpy_cmot = [20, 250e-6];
fixedpy_cmot = [0, 0];
[fitp_y_cmot, ~, ffh_y_cmot, std_err_y_cmot, ~] = fit1D(tofs_cmot, sigmay_cmot, [], sigma_fn, initpy_cmot, fixedpy_cmot, lb_cmot, ub_cmot);

interp_t_cmot = linspace(0, max(tofs_cmot), 300);

% plot
nrows = 2;
ncols = 2;

fighandle_cmot = figure;
subplot(nrows, ncols, 1)
plot(tofs_cmot, peak_ods_cmot, 'o');
grid on;
hold on;
xlabel('TOF (ms)');
ylabel('Peak OD');
% title('CMOT Peak OD vs. Loading Time');
ax = gca;
ax.YLim(1) = 0;
grid on;

subplot(nrows, ncols, 3);
plot(tofs_cmot, atom_num_cmot, 'o');
grid on;
hold on;
xlabel('TOF (ms)');
ylabel('Atom Number');
ax = gca;
ax.YLim(1) = 0;

subplot(nrows, ncols, 2);
plot(tofs_cmot, sigmax_cmot, 'o');
hold on;
plot(interp_t_cmot, ffh_x_cmot(interp_t_cmot));
xlabel('TOF (ms)');
ylabel('sigma x');
grid on;
title(sprintf('sigma_o = %0.0f +/- %0.0f pix, T = %0.0f +/- %0.0f uK',...
    fitp_x_cmot(1), std_err_x_cmot(1), fitp_x_cmot(2)/1e-6, std_err_x_cmot(2)/1e-6));
ax = gca;
ax.YLim(1) = 0;

subplot(nrows, ncols, 4);
plot(tofs_cmot, sigmay_cmot, 'o');
hold on;
plot(interp_t_cmot, ffh_y_cmot(interp_t_cmot));
xlabel('TOF (ms)');
ylabel('sigma y');
grid on;
title(sprintf('sigma_o = %0.0f +/- %0.0f pix, T = %0.0f +/- %0.0f uK',...
    fitp_y_cmot(1), std_err_y_cmot(1), fitp_y_cmot(2)/1e-6, std_err_y_cmot(2)/1e-6));
ax = gca;
ax.YLim(1) = 0;

suptitle(sprintf('%03d cmot temperature', folder));

if saving
    fname = sprintf('%03d_cmot_temperature.fig', folder);
    savefig(fighandle_cmot, fname);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofs_mot = [0.2:0.2:3.2]*1e-3;

indices_mot = 12:22;
indices_mot = indices_mot(1: min(size(vals, 1) - min(indices_mot) + 1, length(indices_mot)));
peak_ods_mot = vals(indices_mot, peakod_index);
% remove any nans
indices_mot = indices_mot(~isnan(peak_ods_mot));
tofs_mot = tofs_mot(1:length(peak_ods_mot));
tofs_mot = tofs_mot(~isnan(peak_ods_mot));
peak_ods_mot = peak_ods_mot(~isnan(peak_ods_mot));
% process other files

atom_num_mot = vals(indices_mot, atomnum_index);
cx_mot = vals(indices_mot, cx_index);
cy_mot = vals(indices_mot, cy_index);
sigmax_mot = abs(vals(indices_mot, sx_index));
sigmay_mot = abs(vals(indices_mot, sy_index));
tofs_mot = tofs_mot(1:length(peak_ods_mot));

min_t_mot = -inf;
max_t_mot = inf;
atom_num_mot = atom_num_mot(tofs_mot < max_t_mot & tofs_mot > min_t_mot);
peak_ods_mot = peak_ods_mot(tofs_mot < max_t_mot & tofs_mot > min_t_mot);
cx_mot = cx_mot(tofs_mot < max_t_mot & tofs_mot > min_t_mot);
cy_mot = cy_mot(tofs_mot < max_t_mot & tofs_mot > min_t_mot);
sigmax_mot = sigmax_mot(tofs_mot < max_t_mot & tofs_mot > min_t_mot);
sigmay_mot = sigmay_mot(tofs_mot < max_t_mot & tofs_mot > min_t_mot);
tofs_mot = tofs_mot(tofs_mot < max_t_mot & tofs_mot > min_t_mot);

kb = 1.3806e-23;
m = 9.9883e-27;
pix_size = 4.65e-6;
mag = 30/250;
sigma_fn = @(P, t) sqrt(P(1)^2 + kb / m * P(2) * t.^2 * (mag/pix_size)^2 );

initpx_mot = [20, 250e-6];
fixedpx_mot = [0, 0];
lb_mot = [0, 0];
ub_mot = [inf, inf];
[fitp_x, ~, ffh_x, std_err_x, ~] = fit1D(tofs_mot, sigmax_mot, [], sigma_fn, initpx_mot, fixedpx_mot, lb_mot, ub_mot);

initpy_mot = [20, 250e-6];
fixedpy_mot = [0, 0];
[fitp_y_mot, ~, ffh_y_mot, std_err_y_mot, ~] = fit1D(tofs_mot, sigmay_mot, [], sigma_fn, initpy_mot, fixedpy_mot, lb_mot, ub_mot);

interp_t_mot = linspace(0, max(tofs_mot), 300);

% plot
nrows = 2;
ncols = 2;

fighandle_mot = figure;
subplot(nrows, ncols, 1)
plot(tofs_mot, peak_ods_mot, 'o');
grid on;
hold on;
xlabel('TOF (ms)');
ylabel('Peak OD');
% title('CMOT Peak OD vs. Loading Time');
ax = gca;
ax.YLim(1) = 0;
grid on;

subplot(nrows, ncols, 3);
plot(tofs_mot, atom_num_mot, 'o');
grid on;
hold on;
xlabel('TOF (ms)');
ylabel('Atom Number');
ax = gca;
ax.YLim(1) = 0;

subplot(nrows, ncols, 2);
plot(tofs_mot, sigmax_mot, 'o');
hold on;
plot(interp_t_mot, ffh_x(interp_t_mot));
xlabel('TOF (ms)');
ylabel('sigma x');
grid on;
title(sprintf('sigma_o = %0.0f +/- %0.0f pix, T = %0.0f +/- %0.0f uK',...
    fitp_x(1), std_err_x(1), fitp_x(2)/1e-6, std_err_x(2)/1e-6));
ax = gca;
ax.YLim(1) = 0;

subplot(nrows, ncols, 4);
plot(tofs_mot, sigmay_mot, 'o');
hold on;
plot(interp_t_mot, ffh_y_mot(interp_t_mot));
xlabel('TOF (ms)');
ylabel('sigma y');
grid on;
title(sprintf('sigma_o = %0.0f +/- %0.0f pix, T = %0.0f +/- %0.0f uK',...
    fitp_y_mot(1), std_err_y_mot(1), fitp_y_mot(2)/1e-6, std_err_y_mot(2)/1e-6));
ax = gca;
ax.YLim(1) = 0;

suptitle(sprintf('%03d mot temperature', folder));

if saving
    fname = sprintf('%03d_mot_temperature.fig', folder);
    savefig(fighandle_mot, fname);
end