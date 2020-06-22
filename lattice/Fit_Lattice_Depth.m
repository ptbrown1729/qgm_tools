function [fit_params, std_errs] = fit_lattice_depth(latt_servo_vs, peaks1, peaks2, peaks3,...
    widths1, widths2, widths3, init_params, fixed_params, use_peak, display_results)
% [fit_paramsj, std_errs] = Fit_Lattice_Depth(latt_servo_vs, peaks1, peaks2, peaks3,...
%     widths1, widths2, widths3, init_params, fixed_params, use_peak, display_results)
%
%%%%%%%%%%%%%%%%%%%
% Fits lattice depth and the attenuation for the retro electric field. So the
% attenuation factor for the retro power is this parameter squared. Does
% this using the locations of the band transitions from the ground band to
% the three d-bands.
%
%%% Arguments
%%%%%%%%%%%%%%%%%%%
% latt_servo_vs: is the servo depth in volts. Fit a scaling factor to this, so
% the lattice depth is LatticeServoVoltage*FitParams(1)
%
% peaks1: is a list of peaks for the lowest d-band, one for each LattServo
% voltage.
%
% widths1: is a list of widths for the lowest d-band points. These are used to
% determine the weighting of each peak position in the fit.
%
% init_params: initial guesses for alttice depth fitting parameters,
% [ScalingFactor,Attenuation Factor for Retro E-field]
%
% fixed_params: A list of ones or zeros. A one will fix the appropriate
% parameter to the value provided in InitParams. A zero will allow that
% parameter to be freely fit.
%
% use_peak: specifies whether each peak is to be used in the fit. Default is
% [1,1,1], which directs the fit to use all of the peaks.
%
%%% Return values
%%%%%%%%%%%%%%%%%%%
% fit_params = [ScalingFactor, Attenuation Factor for Retro E-field]

%experimental data.
% LattServos = [10];
% Peaks1 = [351];
% Width1 = [23];
% Peaks2 = [1];
% Width2 = [3];
% Peaks3 = [380];
% Width3 = [9];
% InitParams = [4,0.47];
% FixedParams = [0,0];
% UsePeak = [1,0,1];

if ~exist('fixed_params', 'var') || isempty(fixed_params)
    fixed_params = [0, 0];
end

if ~exist('use_peak', 'var') || isempty(use_peak)
    use_peak = [1, 1, 1];
end

if ~exist('display_results', 'var') || isempty(display_results)
    display_results = 1;
end

latt_data_path = 'LatticeBandData_Interpolation.txt';
n_bands = 8;
e_recoil = 14.66; % recoil energy in KHz

%Read data lattice depth data
latt_data = dlmread(latt_data_path, ',', 1,0);
depths = latt_data(:,1);
efield_retro_atten = latt_data(:,2);

%Reshape these into meshgrid like shapes. This makes interpolation much
%faster.
n_depths = length(unique(depths));
depths_grid = reshape(depths,[length(depths)/n_depths,n_depths]);
retro_atten_grid = reshape(efield_retro_atten,[length(depths)/n_depths,n_depths]);

%reshape the data to have a similar shape, where the third axis is now the
%various different lattice parameters.
latt_data_grid = reshape(latt_data, [length(depths) / n_depths, n_depths,size(latt_data, 2)]);
dband_centers_khz = latt_data_grid(:, :, 13:15) * e_recoil; %Multiply by Er in KHz to get KHz
all_bw_khz = latt_data_grid(:, :, end - n_bands + 1 : end) * e_recoil;
d_bw_khz = all_bw_khz(:, :, 3:5);

%peak functions
d1Fn = @(depth, atten) interp2(depths_grid, retro_atten_grid, dband_centers_khz(:,:,1), depth, atten);
d2Fn = @(depth, atten) interp2(depths_grid, retro_atten_grid, dband_centers_khz(:,:,2), depth, atten);
d3Fn = @(depth, atten) interp2(depths_grid, retro_atten_grid, dband_centers_khz(:,:,3), depth, atten);

%bw functions
dBW1Fn = @(depth, atten) interp2(depths_grid, retro_atten_grid, d_bw_khz(:,:,1), depth, atten);
dBW2Fn = @(depth, atten) interp2(depths_grid, retro_atten_grid, d_bw_khz(:,:,2), depth, atten);
dBW3Fn = @(depth, atten) interp2(depths_grid, retro_atten_grid, d_bw_khz(:,:,3), depth, atten);

%other interpolated functions
txFn = @(depth, atten) interp2(depths_grid, retro_atten_grid, latt_data_grid(:,:,7), depth, atten);
tyFn = @(depth, atten) interp2(depths_grid, retro_atten_grid, latt_data_grid(:,:,8), depth, atten);
tdiagFn = @(depth, atten) interp2(depths_grid, retro_atten_grid, latt_data_grid(:,:,9), depth, atten);
uFn = @(depth, atten) interp2(depths_grid, retro_atten_grid, latt_data_grid(:,:,10), depth, atten);

%fitting
Fit1 = @(P) (d1Fn(P(1) * latt_servo_vs, P(2)) - peaks1) ./ widths1;
Fit2 = @(P) (d2Fn(P(1) * latt_servo_vs, P(2)) - peaks2) ./ widths2;
Fit3 = @(P) (d3Fn(P(1) * latt_servo_vs, P(2)) - peaks3) ./ widths3;


% fit_fn = @(P) [use_peak(1) * Fit1(P .* (1 - fixed_params) + init_params .* fixed_params),...
%                use_peak(2) * Fit2(P .* (1 - fixed_params) + init_params .* fixed_params),...
%                use_peak(3) * Fit3(P .* (1 - fixed_params) + init_params .* fixed_params)];
% fit_params = lsqnonlin(fit_fn, init_params);

fit_fn = @(P) [use_peak(1) * Fit1(P),...
               use_peak(2) * Fit2(P),...
               use_peak(3) * Fit3(P)];

lbs = [0, 0];
ubs = [100, 1];
[fit_params, std_errs, red_chi_sqr, ~] = lsq_fixedp(fit_fn, init_params, fixed_params, lbs, ubs);

%%%%%%%%%%%%%%
%Display Results
%%%%%%%%%%%%%%
tx = txFn(fit_params(1) * 1, fit_params(2));
ty = tyFn(fit_params(1) * 1, fit_params(2));
tdiag = tdiagFn(fit_params(1) * 1, fit_params(2));
u = uFn(fit_params(1) * 1, fit_params(2));

if display_results
    fprintf('Lattice Depth@10V Servo = %0.1f Er \n', fit_params(1) * 10);
    fprintf('Retro EField Attenuation = %0.3f \n\n', fit_params(2));

    fprintf('Estimated parameters at 1V \n');
    fprintf('tx = %0.2f Er = %0.2fHz\n', tx, tx * e_recoil * 1e3);
    fprintf('ty = %0.2f Er = %0.2fHz\n', ty, ty * e_recoil * 1e3);
    fprintf('tdiag = %0.2f Er = %0.2fHz\n', tdiag, tdiag * e_recoil * 1e3);
    fprintf('u = %0.2f Er = %0.2fHz\n', u, u * e_recoil * 1e3);
    fprintf('u/t = %0.2f \n\n', u / (0.5 * (tx + ty)) );
    
    %plot results
    figure('name', 'Lattice Depth Fit');
    if use_peak(1)
        ph = errorbar(latt_servo_vs, peaks1, widths1, 'ro');
    end
    
    hold on;
    
    if use_peak(2)
        ph = errorbar(latt_servo_vs, peaks2, widths2, 'ro');
    end
    
    if use_peak(3)
        ph = errorbar(latt_servo_vs, peaks3, widths3, 'ro');
    end
    
    ph_thry = errorbar(latt_servo_vs, d1Fn(fit_params(1) * latt_servo_vs, fit_params(2)), dBW1Fn(fit_params(1) * latt_servo_vs, fit_params(2)), 'bo');
    ph_thry = errorbar(latt_servo_vs, d2Fn(fit_params(1) * latt_servo_vs, fit_params(2)), dBW2Fn(fit_params(1) * latt_servo_vs, fit_params(2)), 'bo');
    ph_thry = errorbar(latt_servo_vs, d3Fn(fit_params(1) * latt_servo_vs, fit_params(2)), dBW3Fn(fit_params(1) * latt_servo_vs, fit_params(2)), 'bo');
    grid on;
    hold off;
    xlabel('Lattice Servo Voltage');
    ylabel('Lattice Peaks (KHz)');
    title(sprintf('Depth@10V = %0.1f Er, Retro Efield Attenuation = %0.3f', fit_params(1) * 10, fit_params(2)));
    legend([ph, ph_thry], {'experiment', 'theory'});
end

end
