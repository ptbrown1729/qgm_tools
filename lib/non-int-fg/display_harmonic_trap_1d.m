function fig_handle =  display_harmonic_trap_1d(mu_struct, axes)
% display results from fit_harmonic_trap_1d function

if ~exist('axes', 'var')
    fig_handle = figure;
    axes = gca;
else
    fig_handle = axes.Parent;
end

ph1 = errorbar(axes, mu_struct.radius, mu_struct.mus, mu_struct.mus_unc, 'bo');
hold on;
rad_interp = linspace(0, max(mu_struct.radius), 100);
ph2 = plot(axes, rad_interp, mu_struct.fn_mu_parab(rad_interp), 'r');
ph3 = plot(axes, rad_interp, mu_struct.fn_mu_gauss(rad_interp), 'k');
ph4 = plot(axes, rad_interp, mu_struct.fn_mu_gauss_harmonic_approx(rad_interp), 'k--');

grid on;
xlabel('r');
ylabel('mu (t)');

title(axes, sprintf('mu, omega parab = (2pi) %0.2f(%.0f) sqrt(t/h)\n coeff = %0.5f(%0.0f) \n omega gauss = (2pi) %0.2f(%.0f) sqrt(t/h)\n waist gauss = %0.1f(%.0f) a',...
      mu_struct.omega_parab_sqrt_toverh / (2*pi),...
      mu_struct.omega_parab_sqrt_toverh_unc / (2*pi) * 1e2,...
      mu_struct.fit_params_mu_parab(1) ^ 2,...
      2 * abs(mu_struct.fit_params_mu_parab(1)) * (mu_struct.fit_err_mu_parab(1) * 1e5),...
      mu_struct.omega_gauss_sqrt_toverh / (2*pi),...
      mu_struct.omega_gauss_sqrt_toverh_unc / (2*pi) * 1e2,...
      mu_struct.fit_params_mu_gauss(2),...
      mu_struct.fit_err_mu_gauss(2)));

% TODO: figure out why this isn't working
try
legend([ph1, ph2, ph3, ph4], {'mu', 'parabola', 'gaussian', 'gauss harm apprx'});
catch
end

end