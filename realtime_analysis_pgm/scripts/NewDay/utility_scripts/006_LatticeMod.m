%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters and load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Saving = 1;
close all;

df = DataFolder({2017 1 26 8 1 1}, [], [0, 20]);
% df = DataFolder({2019 1 4 6 1 1}, [], [0, 20]);
df.IndependentVariable = [340 : 2.5 :410];
% df.IndependentVariable = [320 : 2 :430];
[frqs, nums, nums_sdm] = df.showAllVsIndVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit atom number loss features to lorentzians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
min_spacing = min( frqs(2:end) - frqs(1:end-1));
offset_guess = mean(nums);
% fit results
init_params = [ [350, 10, -offset_guess, offset_guess],...
                [375, 3, -offset_guess, 0],...
                [390, 10, -offset_guess, 0] ];
fixed_params = [ [0, 0, 0, 0],...
                 [0, 0, 0, 1],...
                 [0, 0, 0, 1] ];
lbs = [ [min(frqs), min_spacing, -inf, -inf],...
        [min(frqs), min_spacing, -inf, -inf],...
        [min(frqs), min_spacing, -inf, -inf] ];
ubs = [ [max(frqs), max(frqs), inf, inf],...
        [max(frqs), max(frqs), 0, inf],...
        [max(frqs), max(frqs), 0, inf] ];

[fit_params, ~, ffh, std_err] = fit1D(frqs, nums, [], {'lorentzian1D', 'lorentzian1D', 'lorentzian1D'},...
                      init_params, fixed_params, lbs, ubs);

% print results
fprintf('Center Frq = %0.1f KHz, HWHM = %0.1f KHz \n', fit_params(1), fit_params(2));
fprintf('Center Frq = %0.1f KHz, HWHM = %0.1f KHz \n', fit_params(5), fit_params(6));
fprintf('Center Frq = %0.1f KHz, HWHM = %0.1f KHz \n\n', fit_params(9), fit_params(10));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit lattice depth from lorentzian positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[peaks, I] = sort([fit_params(1), fit_params(5), fit_params(9)]);

peak_errs = [std_err(1), std_err(5), std_err(9)];
peak_errs = peak_errs(I);

hwhms = [fit_params(2) ,fit_params(6), fit_params(10)];
hwhms = hwhms(I);

hwhm_errs = [std_err(2), std_err(6), std_err(10)];
hwhm_errs = hwhm_errs(I);

print_results = 1;
[latt_fit_params, latt_param_errs] = Fit_Lattice_Depth(10, peaks(1), peaks(2), peaks(3),...
    hwhms(1), hwhms(2), hwhms(3), [5.9, 0.48], [0, 0], [1, 1, 1], print_results);

% plot results
frq_interp = linspace(min(frqs), max(frqs), 300);

fig_handle = figure;
errorbar(frqs, nums, nums_sdm);
hold on;
plot(frq_interp, ffh(frq_interp));
xlabel('modulation frequency (KHz)');
ylabel('filling');

grid on;

ttl = sprintf('%03d Lattice Modulation, @10V Servo = mV\n peaks = %0.1f(%.0f), %0.1f(%.0f), %0.1f(%.0f) KHz \n hwhm = %0.1f(%.0f), %0.1f(%.0f), %0.1f(%.0f) KHz\n Depth = %0.1f(%.0f) Er at 10V, Efield attenuation = %0.2f(%.0f)',...
                df.Dataset{4},...
                peaks(1), peak_errs(1) * 1e1,...
                peaks(2), peak_errs(2) * 1e1,...
                peaks(3), peak_errs(3) * 1e1,...
                hwhms(1), hwhm_errs(1) * 1e1,...
                hwhms(2), hwhm_errs(2) * 1e1,...
                hwhms(3), hwhm_errs(3) * 1e1,...
                10 * latt_fit_params(1), 10 * latt_param_errs(1) * 1e1,...
                latt_fit_params(2), latt_param_errs(2) * 1e2);
title(ttl, 'FontSize', 10);
% ax = gca;
% ax.FontSize = 10;

if Saving
    save_path = fullfile(getOtherDataPath(today) , sprintf('%03d_Lattice_Modulation.fig', df.Dataset{4}));
    savefig(fig_handle, save_path);
end