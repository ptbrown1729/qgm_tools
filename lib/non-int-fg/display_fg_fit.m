function [fig_handle] = display_fg_fit(fit_struct, mu_struct, fig_handle)
% plot results from fit_dqmc_attractive.m

if ~exist('fig_handle', 'var') || isempty(fig_handle) || ~isvalid(fig_handle)
    fig_handle = figure;
else
    figure(fig_handle);
end

Title = sprintf('T = %0.2f(%.0f) t, mu peak = %0.2f (mu=0 is n=1), h = %0.2f(%.0f) \n Img Fidelity = %0.2f(%.0f), Dbl Fidelity = %0.2f(%.0f), Fit densities n = %0.2f - %0.2f\n expt densities corrected for efficiency',...
    fit_struct.fit_params(1), fit_struct.std_errs(1) * 100,...
    mu_struct.mu_peak,...
    fit_struct.fit_params(2), fit_struct.std_errs(2) * 100,...
    fit_struct.fit_params(3), fit_struct.std_errs(3) * 100,...
    fit_struct.fit_params(4), fit_struct.std_errs(4) * 100,...
    fit_struct.n_cutoffs(1), fit_struct.n_cutoffs(2));

% plot singles density versus density
nrows = 2;
ncols = 4;

% doublons
subplot(nrows, ncols, 1)
ph1 = errorbar(fit_struct.n, fit_struct.nd, fit_struct.nd_err, fit_struct.nd_err, fit_struct.n_err, fit_struct.n_err); 
hold on;
ph2 = plot(fit_struct.n_interp, fit_struct.nd_fit, 'r');
plot(fit_struct.n_interp, fit_struct.nd_fit_maxT, 'r--');
ph3 = plot(fit_struct.n_interp, fit_struct.nd_fit_minT, 'r--');
xlabel('<n>')
ylabel('<d>');
title(sprintf('doubles (fit = %d)', fit_struct.use_param_to_fit(1)));
grid on;
ylim([0, 1]);
legend([ph1, ph2, ph3], {'expt', 'fg' , 'T err'})

% singles
subplot(nrows, ncols, 5)
ph1 = errorbar(fit_struct.n, fit_struct.ns, fit_struct.ns_err, fit_struct.ns_err, fit_struct.n_err, fit_struct.n_err);
hold on;
ph2 = plot(fit_struct.n_interp, fit_struct.ns_fit, 'r');
plot(fit_struct.n_interp, fit_struct.ns_fit_maxT, 'r--');
ph3 = plot(fit_struct.n_interp, fit_struct.ns_fit_minT, 'r--');
xlabel('<n>')
ylabel('<n^s>');
title(sprintf('singles (fit = %d)', fit_struct.use_param_to_fit(2)));
grid on;
ylim([0, 1]);
legend([ph1, ph2, ph3], {'expt', 'fg' , 'T err'})

% plot spin-up correlator versus density
subplot(nrows, ncols, 3)
ph1 = errorbar(fit_struct.n, fit_struct.nup_corr, fit_struct.nup_corr_err, fit_struct.nup_corr_err, fit_struct.n_err, fit_struct.n_err, 'b');
hold on;
ph2 = plot(fit_struct.n_interp, fit_struct.nup_corr_fit, 'b');
plot(fit_struct.n_interp, fit_struct.nup_corr_fit_maxT, 'r--');
plot(fit_struct.n_interp, fit_struct.nup_corr_fit_minT, 'r--');

xlabel('<n>');
ylabel('Corr(0,1)');
title(sprintf('up corr (fit = %d)', fit_struct.use_param_to_fit(3)));
grid on;
legend([ph1, ph2,], {'up corr', 'fg'});
ylim([-0.065, 0.015]);

% spin-dn correlator vs density
subplot(nrows, ncols, 7)
ph1 = errorbar(fit_struct.n, fit_struct.ndn_corr, fit_struct.ndn_corr_err, fit_struct.ndn_corr_err, fit_struct.n_err, fit_struct.n_err, 'b');
hold on;
ph2 = plot(fit_struct.n_interp, fit_struct.ndn_corr_fit, 'r');
% errorbars
plot(fit_struct.n_interp, fit_struct.ndn_corr_fit_maxT, 'r--');
plot(fit_struct.n_interp, fit_struct.ndn_corr_fit_minT, 'r--');

legend([ph1, ph2], {'down corr', 'dqmc'});
xlabel('<n>');
ylabel('Corr(0,1)');
title(sprintf('down corr (fit = %d)', fit_struct.use_param_to_fit(3)));
grid on;
legend([ph1, ph2], {'down corr', 'dqmc'});
ylim([-0.065, 0.015]);

% plot doublon doublon correlator versus density
subplot(nrows, ncols, 2)

ph1 = errorbar(fit_struct.n, fit_struct.nd_corr, fit_struct.nd_corr_err, fit_struct.nd_corr_err, fit_struct.n_err, fit_struct.n_err);
hold on;
ph2 = plot(fit_struct.n_interp, fit_struct.nd_corr_fit, 'r');
plot(fit_struct.n_interp, fit_struct.nd_corr_fit_maxT, 'r--');
plot(fit_struct.n_interp, fit_struct.nd_corr_fit_minT, 'r--');

xlabel('<n>');
ylabel('Corr(0,1)');
title(sprintf('doubles corr (fit = %d)', fit_struct.use_param_to_fit(4)));
grid on;
legend([ph1, ph2,], {'Expt', 'dqmc'});
ylim([-0.065, 0.015]);

% singles correlator versus density
subplot(nrows, ncols, 6)

ph1 = errorbar(fit_struct.n, fit_struct.ns_corr, fit_struct.ns_corr_err, fit_struct.ns_corr_err, fit_struct.n_err, fit_struct.n_err, 'b');
hold on;
ph2 = plot(fit_struct.n_interp, fit_struct.ns_corr_fit, 'r');
ph3 = plot(fit_struct.n_interp, fit_struct.ns_corr_fit_maxT, 'r--');
plot(fit_struct.n_interp, fit_struct.ns_corr_fit_minT, 'r--');

xlabel('<n>');
ylabel('Corr(0,1)');
title(sprintf( 'singles corr (fit = %d)', fit_struct.use_param_to_fit(5)) );
grid on;
legend([ph1, ph2], {'expt', 'fg'});
ylim([-0.065, 0.015]);

% spatially varying chemical potential
ax = subplot(nrows, ncols, 8);
if ~mu_struct.is2d
    display_harmonic_trap_1d(mu_struct, ax);
else
    display_harmonic_trap_2d(mu_struct, ax);
    ax.FontSize = 6;
end

% plot component densities
subplot(nrows, ncols, 4);
ph1 = errorbar(fit_struct.n, fit_struct.nup, fit_struct.nup_err, fit_struct.nup_err, fit_struct.n_err, fit_struct.n_err, 'b');
hold on;
ph2 = plot(fit_struct.n_interp, fit_struct.nup_fit, 'Color', [0.1, 0.5, 1]);
ph3 = errorbar(fit_struct.n, fit_struct.ndn, fit_struct.ndn_err, fit_struct.ndn_err, fit_struct.n_err, fit_struct.n_err, 'r');
ph4 = plot(fit_struct.n_interp, fit_struct.ndn_fit, 'Color', [0.5, 0, 0]);

xlabel('<n>')
ylabel('<n_\sigma>');
title( sprintf('nup and ndn (fit = %d)', fit_struct.use_param_to_fit(6)) );
grid on;
ylim([0, 1]);
legend([ph1, ph2, ph3, ph4], {'nup', 'nup fg', 'ndn', 'ndn fg'})

suptitle(Title);

end