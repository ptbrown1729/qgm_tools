function fit_info_struct = fit_harmonic_trap_2d(mus, mus_unc, radii, x_pos, y_pos)

% mus_unc
if ~exist('mus_unc', 'var') || isempty(mus_unc)
    mus_unc = ones(size(mus));
end

ignore_nans = 1;

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
cx_guess = mean(x_pos(:));
cy_guess = mean(y_pos(:));
omegax_guess = 0.01;
omegay_guess = 0.01;
theta_guess = 0;
bg_guess = max(mus(:));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit FG to 2D parabola
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters are [cx, cy, omega_x, omega_y, theta, background]
fn_parab = @(P, X, Y) parabola2D([P(1), P(2), -P(3)^2, -P(4)^2, P(5), P(6)], X, Y);

init_params_parab = [cx_guess, cy_guess, omegax_guess, omegay_guess, theta_guess, bg_guess];
fixed_params_parab = [0, 0, 0, 0, 0, 0];
lbs_parab = [min(x_pos(:)), min(y_pos(:)), -inf, -inf, 0, -inf];
ubs_parab = [max(x_pos(:)), max(y_pos(:)), inf, inf, pi, inf];

[fit_params_mu_parab, ~, ffh_parab, fit_err_mu_parab] = ...
        fit2D(x_pos, y_pos, mus, 1 ./ mus_unc.^2, fn_parab,...
              init_params_parab, fixed_params_parab,...
              lbs_parab, ubs_parab, '', '', ignore_nans);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit trap to 2D gaussian
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters are [cx, cy, waist_x, waist_y, amp, theta, bg]
fn_gauss = @(P, X, Y) P(5) * (1 - gaussian2D([P(1), P(2), 0.5 * P(3), 0.5 * P(4), 1, P(6), 0], X, Y)) + P(7);
fn_gauss_harmonic_approx = @(P, X, Y) fn_parab([P(1), P(2), sqrt(2 * P(5) / P(3)^2),...
                                                   sqrt(2 * P(5) / P(4)^2), P(6), P(7)],...
                                                   X, Y);

init_params_gauss = [cx_guess, cy_guess, 100, 100, 1, theta_guess, bg_guess];
fixed_params_gauss = [0, 0, 0, 0, 0, 0, 0];
lbs_gauss = [min(x_pos(:)), min(y_pos(:)), 0, 0, -inf, 0, -inf];
ubs_gauss = [max(x_pos(:)), max(y_pos(:)), 200e-6/a, 200e-6/a, inf, pi, inf];
ignore_nans = 1;

[fit_params_gauss, ~, ffh_gauss, fit_err_gauss] = ...
    fit2D(x_pos, y_pos, mus, 1 ./ mus_unc.^2, fn_gauss,...
          init_params_gauss, fixed_params_gauss,...
          lbs_gauss, ubs_gauss, '', '', ignore_nans);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract trapping frequency information
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% trapping frequency from parabolic fit
omegax_parab_sqrt_toverh =  factor * abs(fit_params_mu_parab(3));
omegax_parab_sqrt_toverh_unc = factor * fit_err_mu_parab(3);

omegay_parab_sqrt_toverh =  factor * abs(fit_params_mu_parab(4));
omegay_parab_sqrt_toverh_unc = factor * fit_err_mu_parab(4);

% trapping frequency from gaussian fit
omegax_gauss_sqrt_toverh = factor * sqrt( 2 * abs(fit_params_gauss(5)) / fit_params_gauss(3)^2 );
omegax_gauss_sqrt_toverh_unc = omegax_gauss_sqrt_toverh * ...
                               sqrt( (fit_err_gauss(5) / fit_params_gauss(5) )^2 ...
                                   + (2 * fit_err_gauss(3) / fit_params_gauss(3) )^2 );
          
omegay_gauss_sqrt_toverh = factor * sqrt( 2 * abs(fit_params_gauss(5)) / fit_params_gauss(4)^2 );
omegay_gauss_sqrt_toverh_unc = omegay_gauss_sqrt_toverh * ...
                               sqrt( (fit_err_gauss(5) / fit_params_gauss(5) )^2 ...
                                   + (2 * fit_err_gauss(4) / fit_params_gauss(4) )^2 );
         

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% store mu results
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fit_info_struct = struct();
fit_info_struct.is2d = 1;


% mu and spatial information
fit_info_struct.radius = radii;
fit_info_struct.xpos = x_pos;
fit_info_struct.ypos = y_pos;
fit_info_struct.mus = mus;
fit_info_struct.mus_unc = mus_unc;
% TODO: should mu_peak instead by the value extrapolated to the trap
% center?
fit_info_struct.mu_peak = max(mus(:));

% fit to parabola information
fit_info_struct.fn_mu_parab = ffh_parab;

fit_info_struct.fit_params_mu_parab = fit_params_mu_parab;
fit_info_struct.fit_err_mu_parab = fit_err_mu_parab;

fit_info_struct.omegax_parab_sqrt_toverh = omegax_parab_sqrt_toverh;
fit_info_struct.omegax_parab_sqrt_toverh_unc = omegax_parab_sqrt_toverh_unc;
fit_info_struct.omegay_parab_sqrt_toverh = omegay_parab_sqrt_toverh;
fit_info_struct.omegay_parab_sqrt_toverh_unc = omegay_parab_sqrt_toverh_unc;

% fit to gaussian information
fit_info_struct.fn_mu_gauss = ffh_gauss;
fit_info_struct.fn_mu_gauss_harmonic_approx = @(X, Y) fn_gauss_harmonic_approx(fit_params_gauss, X, Y);

fit_info_struct.fit_params_mu_gauss = fit_params_gauss;
fit_info_struct.fit_err_mu_gauss = fit_err_gauss;

fit_info_struct.omegax_gauss_sqrt_toverh = omegax_gauss_sqrt_toverh;
fit_info_struct.omegax_gauss_sqrt_toverh_unc = omegax_gauss_sqrt_toverh_unc;
fit_info_struct.omegay_gauss_sqrt_toverh = omegay_gauss_sqrt_toverh;
fit_info_struct.omegay_gauss_sqrt_toverh_unc = omegay_gauss_sqrt_toverh_unc;

end