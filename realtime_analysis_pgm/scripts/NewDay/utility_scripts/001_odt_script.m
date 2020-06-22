saving = 1;
folder = 1;
close all

HF1 = 5.04;
HF3 = 7.69;
% loading frq
fname = fullfile(getDataSetPath(datenum(2019, 3, 11), folder), 'log.txt');
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


% define indices and dependent variables for our different measurements
% no evaporation
tofs_noevap = [0:20:200]*1e-6;
index_state2 = 2;
index_state3 = 3;
indices_noevap = [4:14];
% evaporation
tofs_evap = [0:200:2000]*1e-6;
index_state3_evap = 26; %
indices_evap = [27:37];
% trapping frequency
frqs_trpfrq = [6:1:15]; %[5:1:15];
indices_trpfrq = 15:24; %[14:24];
% imaging resonance state 1
frqs_state1 = [.95:.01:1.05]*HF1;
indices_state1 = [39:49];
% imaging resonance state 3
frqs_state3 = [.95:.01:1.05]*HF3;
indices_state3 = [50:60];

% constants
kb = 1.3806e-23;
m = 9.9883e-27;
hw_bin_size = 4;
%4 for hardware binning, software binning, but coordinates account for that
pix_size = 6.5e-6 * hw_bin_size; 
mag = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODT, no evaporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state 2 and 3
peak_od_state2 = vals(index_state2, peakod_index);
atom_num_state2 = vals(index_state2, atomnum_index);

peak_od_state3 = vals(index_state3, peakod_index);
atom_num_state3 = vals(index_state3, atomnum_index);

% state 1s
indices_noevap = indices_noevap(1: min(size(vals, 1) - min(indices_noevap) + 1, length(indices_noevap)));
peak_ods = vals(indices_noevap, 5);
% remove any nans
indices_noevap = indices_noevap(~isnan(peak_ods));
tofs_noevap = tofs_noevap(1:length(peak_ods));
tofs_noevap = tofs_noevap(~isnan(peak_ods));
peak_ods = peak_ods(~isnan(peak_ods));
% process other files
atom_num = vals(indices_noevap, atomnum_index);
cx = vals(indices_noevap, cx_index);
cy = vals(indices_noevap, cy_index);
sigmax = abs(vals(indices_noevap, sx_index));
sigmay = abs(vals(indices_noevap, sy_index));
tofs_noevap = tofs_noevap(1:length(peak_ods));

min_t = -inf;
max_t = inf;
atom_num = atom_num(tofs_noevap < max_t & tofs_noevap > min_t);
peak_ods = peak_ods(tofs_noevap < max_t & tofs_noevap > min_t);
sigmax = sigmax(tofs_noevap < max_t & tofs_noevap > min_t);
sigmay = sigmay(tofs_noevap < max_t & tofs_noevap > min_t);
tofs_noevap = tofs_noevap(tofs_noevap < max_t & tofs_noevap > min_t);

sigma_fn = @(P, t) sqrt(P(1)^2 + kb / m * P(2) * t.^2 * (mag/pix_size)^2 );

initpx = [20, 250e-6];
fixedpx = [0, 0];
lb = [0, 0];
ub = [inf, inf];
[fitp_x, ~, ffh_x, std_err_x, ~] = fit1D(tofs_noevap, sigmax, [], sigma_fn, initpx, fixedpx, lb, ub);

initpy = [20, 250e-6];
fixedpy = [0, 0];
[fitp_y, ~, ffh_y, std_err_y, ~] = fit1D(tofs_noevap, sigmay, [], sigma_fn, initpy, fixedpy, lb, ub);

interp_t = linspace(0, max(tofs_noevap), 300);

% plot temperature for ODT, no evaporation
nrows = 2;
ncols = 2;

fighandle = figure;

subplot(nrows, ncols, 1)

plot(tofs_noevap, peak_ods, 'o');
hold on;
plot(0, peak_od_state3, 'o');
plot(0, peak_od_state2, 'o');
grid on;
hold on;
xlabel('TOF (ms)');
ylabel('Peak OD');
ax = gca;
ax.YLim(1) = 0;
grid on;
legend({'state 1', 'state 3', 'state 2'});

subplot(nrows, ncols, 3);
plot(tofs_noevap, atom_num, 'o');
hold on;
plot(0, atom_num_state3, 'o');
% plot(0, atom_num_state2, 'o');
grid on;
xlabel('TOF (ms)');
ylabel('Atom Number');
ax = gca;
ax.YLim(1) = 0;

subplot(nrows, ncols, 2);
plot(tofs_noevap, sigmax, 'o');
hold on;
plot(interp_t, ffh_x(interp_t));
xlabel('TOF (ms)');
ylabel('sigma x');
grid on;
title(sprintf('sigma_o = %0.0f +/- %0.0f pix, T = %0.0f +/- %0.0f uK',...
    fitp_x(1), std_err_x(1), fitp_x(2)/1e-6, std_err_x(2)/1e-6));
ax = gca;
ax.YLim(1) = 0;

subplot(nrows, ncols, 4);
plot(tofs_noevap, sigmay, 'o');
hold on;
plot(interp_t, ffh_y(interp_t));
xlabel('TOF (ms)');
ylabel('sigma y');
grid on;
title(sprintf('sigma_o = %0.0f +/- %0.0f pix, T = %0.0f +/- %0.0f uK',...
    fitp_y(1), std_err_y(1), fitp_y(2)/1e-6, std_err_y(2)/1e-6));
ax = gca;
ax.YLim(1) = 0;

suptitle(sprintf('%03d ODT temperature, no evaporation', folder));

if saving
    fname = sprintf('%03d_odt_temperature_noevap.fig', folder);
    savefig(fighandle, fname);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODT Evaporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% state 3s
peak_od_state3_evap = vals(index_state3_evap, peakod_index);
atom_num_state3_evap = vals(index_state3_evap, atomnum_index);

% state 1s
indices_evap = indices_evap(1: min(size(vals, 1) - min(indices_evap) + 1, length(indices_evap)));
peak_ods_evap = vals(indices_evap, peakod_index);
% remove any nans
indices_evap = indices_evap(~isnan(peak_ods_evap));
tofs_evap = tofs_evap(1:length(peak_ods_evap));
tofs_evap = tofs_evap(~isnan(peak_ods_evap));
peak_ods_evap = peak_ods_evap(~isnan(peak_ods_evap));
% process other files

atom_num_evap = vals(indices_evap, atomnum_index);
cx_evap = vals(indices_evap, cx_index);
cy_evap = vals(indices_evap, cy_index);
sigmax_evap = abs(vals(indices_evap, sx_index));
sigmay_evap = abs(vals(indices_evap, sy_index));
tofs_evap = tofs_evap(1:length(peak_ods_evap));

min_t_evap = -inf;
max_t_evap = 1.9e-3;
atom_num_evap = atom_num_evap(tofs_evap < max_t_evap & tofs_evap > min_t_evap);
peak_ods_evap = peak_ods_evap(tofs_evap < max_t_evap & tofs_evap > min_t_evap);
sigmax_evap = sigmax_evap(tofs_evap < max_t_evap & tofs_evap > min_t_evap);
sigmay_evap = sigmay_evap(tofs_evap < max_t_evap & tofs_evap > min_t_evap);
tofs_evap = tofs_evap(tofs_evap < max_t_evap & tofs_evap > min_t_evap);

sigma_fn = @(P, t) sqrt(P(1)^2 + kb / m * P(2) * t.^2 * (mag/pix_size)^2 );

initpx_evap = [20, 250e-6];
fixedpx_evap = [0, 0];
lb = [0, 0];
ub = [inf, inf];
[fitp_x_evap, ~, ffh_x_evap, std_err_x_evap, ~] = fit1D(tofs_evap, sigmax_evap, [], sigma_fn, initpx_evap, fixedpx_evap, lb, ub);

initpy_evap = [20, 250e-6];
fixedpy_evap = [0, 0];
[fitp_y_evap, ~, ffh_y_evap, std_err_y_evap, ~] = fit1D(tofs_evap, sigmay_evap, [], sigma_fn, initpy_evap, fixedpy_evap, lb, ub);

interp_t_evap = linspace(0, max(tofs_evap), 300);

% plot temperature for ODT after evaporation
nrows = 2;
ncols = 2;

fighandle_evap = figure;

subplot(nrows, ncols, 1)

plot(tofs_evap, peak_ods_evap, 'o');
hold on;
plot(0, peak_od_state3_evap, 'o');
grid on;
xlabel('TOF (ms)');
ylabel('Peak OD');
% title('CMOT Peak OD vs. Loading Time');
ax = gca;
ax.YLim(1) = 0;
grid on;
legend({'state 1', 'state 3'});

subplot(nrows, ncols, 3);
plot(tofs_evap, atom_num_evap, 'o');
hold on;
plot(0, atom_num_state3_evap, 'o');
grid on;
xlabel('TOF (ms)');
ylabel('Atom Number');
ax = gca;
ax.YLim(1) = 0;

subplot(nrows, ncols, 2);
plot(tofs_evap, sigmax_evap, 'o');
hold on;
plot(interp_t_evap, ffh_x_evap(interp_t_evap));
xlabel('TOF (ms)');
ylabel('sigma x');
grid on;
title(sprintf('sigma_o = %0.0f +/- %0.0f pix, T = %0.2f +/- %0.2f uK',...
    fitp_x_evap(1), std_err_x_evap(1), fitp_x_evap(2)/1e-6, std_err_x_evap(2)/1e-6));
ax = gca;
ax.YLim(1) = 0;

subplot(nrows, ncols, 4);
plot(tofs_evap, sigmay_evap, 'o');
hold on;
plot(interp_t_evap, ffh_y_evap(interp_t_evap));
xlabel('TOF (ms)');
ylabel('sigma y');
grid on;
title(sprintf('sigma_o = %0.0f +/- %0.0f pix, T = %0.2f +/- %0.2f uK',...
    fitp_y_evap(1), std_err_y_evap(1), fitp_y_evap(2)/1e-6, std_err_y_evap(2)/1e-6));
ax = gca;
ax.YLim(1) = 0;

suptitle(sprintf('%03d odt temperature, after evaporation', folder));

if saving
    fname = sprintf('%03d_odt_temperature_evap.fig', folder);
    savefig(fighandle_evap, fname);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODT trapping frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indices_trpfrq = indices_trpfrq(1: min(size(vals, 1) - min(indices_trpfrq) + 1, length(indices_trpfrq)));
peak_ods_trpfrq = vals(indices_trpfrq, peakod_index);
% remove any nans
indices_trpfrq = indices_trpfrq(~isnan(peak_ods_trpfrq));
frqs_trpfrq = frqs_trpfrq(1:length(peak_ods_trpfrq));
frqs_trpfrq = frqs_trpfrq(~isnan(peak_ods_trpfrq));
peak_ods_trpfrq = peak_ods_trpfrq(~isnan(peak_ods_trpfrq));
% process other files
peak_ods_trpfrq = vals(indices_trpfrq, peakod_index);
cx = vals(indices_trpfrq, cx_index);
cy = vals(indices_trpfrq, cy_index);
sigmax_trpfrq = abs(vals(indices_trpfrq, sx_index));
sigmay_trpfrq = abs(vals(indices_trpfrq, sy_index));
frqs_trpfrq = frqs_trpfrq(1:length(peak_ods_trpfrq));

init_params_trpfrq = [5, 2, -0.5, 0.75];
fixed_params_trpfrq = [0, 0, 0, 0];
lb_trpfrq = [0, -inf, -inf, -inf];
ub_trpfrq = [inf, inf, 0, inf];
[fp_trpfrq, ~, ffh_trpfrq, stderr_trpfrq, ~] = fit1D(frqs_trpfrq, peak_ods_trpfrq, [], {'lorentzian1D'}, init_params_trpfrq, fixed_params_trpfrq, lb_trpfrq, ub_trpfrq);
interp_frqs_trpfrq = linspace(min(frqs_trpfrq),max(frqs_trpfrq),300);

% plot trapping frequency measurements
nrows = 2;
ncols = 2;

fighandle_trpfrq = figure;
subplot(nrows, ncols, 1)
plot(frqs_trpfrq, peak_ods_trpfrq, 'o');
grid on;
hold on;
plot(interp_frqs_trpfrq, ffh_trpfrq(interp_frqs_trpfrq));
xlabel('Freq (KHz)');
ylabel('Peak OD');

ax = gca;
ax.YLim(1) = 0;
grid on;
title(sprintf('Center = %0.1f +/- %0.1f KHz, HWHM = %0.2f +/- %0.2f KHz',...
    fp_trpfrq(1), stderr_trpfrq(1), fp_trpfrq(2), stderr_trpfrq(2)));

subplot(nrows, ncols, 2);
plot(frqs_trpfrq, sigmax_trpfrq, 'o');
hold on;

xlabel('Freq (KHz)');
ylabel('sigma x');
grid on;

ax = gca;
ax.YLim(1) = 0;

subplot(nrows, ncols, 3);
plot(frqs_trpfrq, sigmay_trpfrq, 'o');
hold on;

xlabel('Freq (KHz)');
ylabel('sigma y');
grid on;

ax = gca;
ax.YLim(1) = 0;

suptitle(sprintf('%03d odt trap frq', folder));

if saving
    fname = sprintf('%03d_odt_trapping_frq.fig', folder);
    savefig(fighandle_trpfrq, fname);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State 1 resonance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% state 1s
indices_state1 = indices_state1(1: min(size(vals, 1) - min(indices_state1) + 1, length(indices_state1)));
peak_ods_state1 = vals(indices_state1, peakod_index);
% remove any nans
indices_state1 = indices_state1(~isnan(peak_ods_state1));
frqs_state1 = frqs_state1(1:length(peak_ods_state1));
frqs_state1 = frqs_state1(~isnan(peak_ods_state1));
peak_ods_state1 = peak_ods_state1(~isnan(peak_ods_state1));

% fit to lorentzian
initp_state1 = [mean(frqs_state1), 0.3, 0.5, 0];
fixedp_state1 = [0, 0, 0, 1];
lb_state1 = [min(frqs_state1), 0, 0, 0];
ub_state1 = [max(frqs_state1), 3, 2, 0];
[fitp_state1, ~, ffh_state1, std_err_state1, ~] = ...
            fit1D(frqs_state1, peak_ods_state1,...
                  [], {'lorentzian1D'}, initp_state1, fixedp_state1, lb_state1, ub_state1);
interp_frq_state1 = linspace(min(frqs_state1), max(frqs_state1), 300);

fighandle_state_res = figure;
plot(frqs_state1, peak_ods_state1, 'o');
hold on;
plot(interp_frq_state1, ffh_state1(interp_frq_state1));
xlabel('Frequency (V)');
ylabel('Peak OD');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State 3 resonance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state 3s
indices_state3 = indices_state3(1: min(size(vals, 1) - min(indices_state3) + 1, length(indices_state3)));
peak_ods_state3 = vals(indices_state3, peakod_index);
% remove any nans
indices_state3 = indices_state3(~isnan(peak_ods_state3));
frqs_state3 = frqs_state3(1:length(peak_ods_state3));
frqs_state3 = frqs_state3(~isnan(peak_ods_state3));
peak_ods_state3 = peak_ods_state3(~isnan(peak_ods_state3));

% fit to lorentzian
initp_state3 = [mean(frqs_state3), 0.3, 0.5, 0];
fixedp_state3 = [0, 0, 0, 1];
lb_state3 = [min(frqs_state3), 0, 0, 0];
ub_state3 = [max(frqs_state3), 3, 2, 0];
[fitp_state3, ~, ffh_state3, std_err_state3, ~] = ...
               fit1D(frqs_state3, peak_ods_state3,...
                     [], {'lorentzian1D'}, initp_state3, fixedp_state3, lb_state3, ub_state3);
interp_frq_state3 = linspace(min(frqs_state3), max(frqs_state3), 300);


plot(frqs_state3, peak_ods_state3, 'o');
hold on;
plot(interp_frq_state3, ffh_state3(interp_frq_state3));

title(sprintf('%03d, state 1 and 3 resonances\n 1s = %0.2f(%.0f)V, 3s = %0.2f(%.0f)V',...
                            folder, fitp_state1(1), round(std_err_state1(1)*1e2),...
                                    fitp_state3(1), round(std_err_state3(1)*1e2)));

if saving
    fname = sprintf('%03d_state1_state3_resonance.fig', folder);
    savefig(fighandle_state_res, fname);
end