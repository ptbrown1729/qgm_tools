function fit_info_struct = fit_harmonic_trap_1d(mus, mus_unc, radii, x_pos, y_pos)

% mus_unc
if ~exist('mus_unc', 'var') || isempty(mus_unc)
    mus_unc = ones(size(mus));
end

if ~exist('x_pos', 'var')
    x_pos = '';
end

if ~exist('y_pos', 'var')
    y_pos = '';
end

mus_unc(mus_unc == 0) = 1;

% define constants
m = 9.9883e-27;
a = 750e-9;
h = 6.626e-34;
% Suppose you have the quantity W obtained from k^2 * (r/a)^2 = V(r/a). Then
% the trapping frequency is given by w = factor * k * sqrt(t/h)
factor = sqrt(2 * h / a^2 /m);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% guess fit parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega_guess = sqrt( (mus(1) - mus(end)) / (radii(1)^2 - radii(end)^2) );

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit trap to 1D parabola
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn_parab = @(P, X) P(2) - P(1)^2 * X.^2;

init_params_parab = [omega_guess, max( mus(:) )];
fixed_params_parab = [0, 0];
lbs_parab = [-inf, -inf];
ubs_parab = [0, inf];

[fit_params_parab, ~, ffh_parab, fit_err_parab, ~] = ...
        fit1D(radii, mus, 1 ./ mus_unc.^2, fn_parab, init_params_parab,...
              fixed_params_parab, lbs_parab, ubs_parab);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit trap to 1D gaussian
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters are trapping frequency and waist
fn_mu_gauss = @(P, X) P(3) +  P(1)^2 * P(2)^2 / 2 * ( exp( - 2 * X.^2 / P(2)^2 ) - 1 );
fn_mu_gauss_harmonic_approx = @(P, X) P(3) - P(1)^2 * X.^2;

waist_guess = 100;
init_params_gauss = [omega_guess, waist_guess, max(mus(:))];
fixed_params_gauss = [0, 0, 0];
lbs_gauss = [-inf, 0, -inf];
ubs_gauss = [inf, inf, inf];

[fit_params_gauss, ~, ffh_gauss, fit_err_gauss, ~] = ...
    fit1D(radii, mus, 1 ./ mus_unc.^2, fn_mu_gauss, init_params_gauss,...
          fixed_params_gauss, lbs_gauss, ubs_gauss);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract trapping frequency information
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trapping frequency from parabolic fit
omega_parab_sqrt_toverh =  factor * abs(fit_params_parab(1));
omega_parab_sqrt_toverh_unc = factor * fit_err_parab(1);
% trapping frequency from gaussian fit
omega_gauss_sqrt_toverh = factor * abs( fit_params_gauss(1) );
omega_gauss_sqrt_toverh_unc = factor * fit_err_gauss(1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% store mu results
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fit_info_struct = struct();
fit_info_struct.is2d = 0;

% mu and spatial information
fit_info_struct.radius = radii;
fit_info_struct.xpos = x_pos;
fit_info_struct.ypos = y_pos;
fit_info_struct.mus = mus;
fit_info_struct.mus_unc = mus_unc;
% TODO: should mu_peak instead by the value extrapolated to the trap
% center?
fit_info_struct.mu_peak = max(mus);

% fit to parabola information
fit_info_struct.fn_mu_parab = ffh_parab;

fit_info_struct.fit_params_mu_parab = fit_params_parab;
fit_info_struct.fit_err_mu_parab = fit_err_parab;

fit_info_struct.omega_parab_sqrt_toverh = omega_parab_sqrt_toverh;
fit_info_struct.omega_parab_sqrt_toverh_unc = omega_parab_sqrt_toverh_unc;

% fit to gaussian information
fit_info_struct.fn_mu_gauss = ffh_gauss;
fit_info_struct.fn_mu_gauss_harmonic_approx = @(X) fn_mu_gauss_harmonic_approx(fit_params_gauss, X);

fit_info_struct.fit_params_mu_gauss = fit_params_gauss;
fit_info_struct.fit_err_mu_gauss = fit_err_gauss;

fit_info_struct.omega_gauss_sqrt_toverh = omega_gauss_sqrt_toverh;
fit_info_struct.omega_gauss_sqrt_toverh_unc = omega_gauss_sqrt_toverh_unc;

end