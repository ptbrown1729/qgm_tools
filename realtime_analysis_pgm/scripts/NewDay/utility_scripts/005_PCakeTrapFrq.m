saving = 1;
folder = 5;
date = {2019 1 4};

% bin_edges = sqrt(linspace(0, 40^2, 40));
bin_edges = [0, 20];

df = DataFolder({date{1} date{2} date{3} folder 1 1}, [], bin_edges);
df.IndependentVariable = [20:2:50];

[frqs_unique, all_dens_unique, all_dens_unc_unique,...
    num_points, std_dev, fig_handle2] = df.showAllVsIndVar;


% fit to lorentzian
init_p = [37, 5, -100, 100];
fixed_p = [0, 0, 0, 0];

delta_frqs = frqs_unique(2:end) - frqs_unique(1:end-1);
lbs = [min(frqs_unique), min(delta_frqs), -inf, -inf];
ubs = [max(frqs_unique), max(frqs_unique) - min(frqs_unique), 0, inf]

[fit_params, ~, ffh, std_errs] = fit1D(frqs_unique, all_dens_unique, 1 ./ all_dens_unc_unique .^2,...
                            {'lorentzian1D'}, init_p, fixed_p, lbs, ubs);
interp_frqs = linspace(min(frqs_unique), max(frqs_unique), 300);

fig_handle = figure;
errorbar(frqs_unique, all_dens_unique, all_dens_unc_unique);
hold on;
plot(interp_frqs, ffh(interp_frqs));
grid on;
xlabel('Modulation frequency (KHz)');
ylabel('Filling');

ylim([0, 1]);

suptitle( sprintf('%s\n trap freq = %0.1f(%.0f) KHz, HWHM = %0.1f(%0.0f) KHz',...
    df.identifier, fit_params(1)/2, std_errs(1)/2 * 1e1,...
    fit_params(2)/2, std_errs(2)/2 * 1e1...
    ));


if saving
    fname = sprintf('%03d_pancake_trap_frequency.fig', df.Dataset{4});
    savefig(fig_handle, fname);
end