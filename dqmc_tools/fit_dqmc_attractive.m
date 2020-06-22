function [fit_params, std_errs, redchi_sqr, mu_struct, fit_struct] = fit_dqmc(...
    ups, downs, singles, doubles, use_quantity_to_fit, n_cutoffs,...
    init_params, fixed_params, theory_fns_struct, use_az_avg, img_bin_size, fig_handle)
% Fit various densities and correlators versus total density to determine
% Hubbard interaction, temperature, chemical potential imbalance, imaging
% fidelity, and doublon imaging fidelity in a trap independent way.
%
% TODO: add second neighbor correlators as option
%
% Arguments:
%---------------------------
% ups: DataFolder class object representing spin-up density. If one of ups
% or downs is given and the other is omitted, we assume that the gas is
% balanced and use the given DataFolder object for both ups and downs.
%
% downs: DataFolder class object representing spin-down density.
%
% singles: DataFolder class object representing singles density. May be
% safely omitted.
%
% doubles: DataFolder class object representing doubles density. May be
% safely omitted.
%
% use_quantity_to_fit: A six element boolean array, e.g. [1, 1, 1, 1, 1, 1] 
% which selects which of densities and correlators to use during fitting.
% The order is {'<d>, <n^s>, <nup(1)nup(0)>_c, <d(1)d(0)>_c, <n^s(1)n^s(0)>_c,
%                <nup>, <ndn>, <ndn(1)ndn(0)>_c}.
%
% n_cutoffs: [nmin, nmax] limit the densities that are included in the fit
%
% init_params: [U, T, h, imaging fidelity, doublon fidelity]. Note that the
% doublon fidelity is only the transfer efficiency of the doublon process.
% The imaging fidelity is also taken into account for the doublons.
%
% fixed_params: specifies which of the initial parameters are to be fixed. It
% should be the same size as InitParams, with a zero for any parameter that
% is allowed to vary, and a 1 for any parameter that is not allowed to vary.
% If a one is given, the parameter will be fixed at the value provided in
% InitParams.
%
% theory_fn_structure: A .mat file containing interpolating functions
% created from DQMC data. The work flow to create such a dataset is first
% run the DQMC, then aggregate the output files using aggregate_dqmc_output_files(),
% next arrange the datapoints in a grid with convert_dqmc_to_grid(), and
% use this grid to create the interpolating functions with create_dqmc_interpolating_fn().
%
% use_az_avg: boolean. 1 to use azimuthal average data or 0 to use 2D
% density.
%
% img_bin_size: this variable only has an effect if use_az_avg=0. In that
% case, it determines how to bin the 2D density and correlator data from
% ups, downs, singles, and doubles before using it in the fitting
% function.
%
% fig_handles: previously the figure handle where the results of this
% function would be displayed. This style has been deprecated. Now, you can
% use the output structures as arguments to the display_dqmc_fit()
% function.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Argument checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('use_az_avg', 'var') || isempty(use_az_avg)
    use_az_avg = 1;
end

if ~exist('img_bin_size', 'var') || isempty(img_bin_size)
    img_bin_size = 4;
end

if ~exist('use_quantity_to_fit', 'var') || isempty(use_quantity_to_fit)
    use_quantity_to_fit = [1, 1, 1, 1, 1, 1];
end

if length(use_quantity_to_fit) == 4
    % preserve old scripts functionality
    use_quantity_to_fit = horzcat(use_quantity_to_fit, [0, 0]);
elseif length(use_quantity_to_fit) == 5
        use_quantity_to_fit = horzcat(use_quantity_to_fit, 0);
end

if length(use_quantity_to_fit) ~= 6
    error('use_quantity_to_fit must be of length 6');
end

if isempty(doubles)
    use_quantity_to_fit(1) = 0;
    use_quantity_to_fit(4) = 0;
    warning('no doubles file provided, not fitting doubles density or correlations');
end

if isempty(singles)
    use_quantity_to_fit(2) = 0;
    use_quantity_to_fit(5) = 0;
    warning('no singles file provided, not fitting singles density');
end

if length(n_cutoffs) > 2 || ( n_cutoffs(2) <= n_cutoffs(1) )
    error('problem with n_cutoffs in fitDQMC');
end

if ~exist('n_cutoffs', 'var') || isempty(n_cutoffs)
    % n can't exceed fidelity value
    n_cutoffs = [0.1, init_params(2)];
end

if length(init_params) == 4
    init_params = horzcat(init_params, 1);
    fixed_params = horzcat(fixed_params, 1);
end

if length(init_params) ~= 5
    error('init_params must have length 5: [U/t, T/t, h/t, imaging fidelity, doublon fidelity]');
end

if ~exist('fixed_params', 'var') || isempty(fixed_params)
    fixed_params = zeros(size(init_params));
end

if length(fixed_params) ~= 5
    error('fixed_params must have length 5.');
end 

if ischar(theory_fns_struct)
    theory_fns_struct = load(theory_fns_struct);
end

% TODO: finish with getting axes to work
% if ~exist('axes', 'var') || isempty(axes)
%     fighandle = figure;
%     axes = gca;
% else
%     fighandle = axes.Parent;
% end

% account for possibility we may not have all of ones, singles, doubles,
% etc.
if isempty(downs) && ~isempty(ups)
    downs = ups;
elseif isempty(ups) && ~isempty(downs)
    ups = downs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract densities and correlators from DataFolder instances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if use_az_avg
    % Get information from azimuthal average
    [nup, nup_err, nup_corr, nup_corr_err] = get_azavg_quantities_from_DataFolder(ups);
    [ndn, ndn_err, ndn_corr, ndn_corr_err] = get_azavg_quantities_from_DataFolder(downs);

    % singles
    if ~isempty(singles)
        [ns, nserr, ns_corr, ns_corr_err] = get_azavg_quantities_from_DataFolder(singles);
    else
        ns = nan(size(nup));
        nserr = nan(size(ns));
        ns_corr = nan(size(ns));
        ns_corr_err = nan(size(ns));
    end
    % doubles
    if ~isempty(doubles)
        [d, derr, ddcorr, ddcorr_err] = get_azavg_quantities_from_DataFolder(doubles);
    else
        d = nan(size(nup));
        derr = nan(size(d));
        ddcorr = nan(size(d));
        ddcorr_err = nan(size(d));
    end

    % spatial position
    radius_full = ups.RadialPos;
    
    [xx_full, yy_full] = meshgrid( 1 : ups.ImgCropSize, 1 : ups.ImgCropSize ); 
    xx_full = xx_full - ups.Cx_AzAvg;
    yy_full = yy_full - ups.Cy_AzAvg;
    
else
    % Get information from 2D cloud
    % this mode is a work in progress.    
    [nup, nup_err, nup_corr, nup_corr_err] = get_2dbin_quantities_from_DataFolder(ups, img_bin_size);
    [ndn, ndn_err, ndn_corr, ndn_corr_err] = get_2dbin_quantities_from_DataFolder(downs, img_bin_size);

    % singles
    if ~isempty(singles)
        [ns, nserr, ns_corr, ns_corr_err] = get_2dbin_quantities_from_DataFolder(singles, img_bin_size);
    else
        ns = nan(size(nup));
        nserr = nan(size(ns));
        ns_corr = nan(size(ns));
        ns_corr_err = nan(size(ns));
    end
    % doubles
    if ~isempty(doubles)
        [d, derr, ddcorr, ddcorr_err] = get_2dbin_quantities_from_DataFolder(doubles, img_bin_size);
    else
        d = nan(size(nup));
        derr = nan(size(d));
        ddcorr = nan(size(d));
        ddcorr_err = nan(size(d));
    end
    
    % spatial position
    [xx_full, yy_full] = meshgrid( 1 : ups.ImgCropSize, 1 : ups.ImgCropSize ); 
    xx_full = xx_full - ups.Cx_AzAvg;
    yy_full = yy_full - ups.Cy_AzAvg;
    
    [xx_full, ~, xx_std, ~] = binImg(xx_full, img_bin_size, img_bin_size, 'Normal');
    [yy_full, ~, yy_std, ~] = binImg(yy_full, img_bin_size, img_bin_size, 'Normal');
    radius_full = sqrt( xx_full.^2 + yy_full.^2);
    
end

% total density
n = nup + ndn;
nerr = sqrt( nup_err .^ 2 + ndn_err .^ 2);
% magnetization
m = nup - ndn;
merr = sqrt(nup_err .^ 2 + ndn_err .^2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract functions from theory function structure.
% handle possibility some functions may not be there. If they aren't there,
% then don't try and use them for fitting!
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thry = theory_fns_struct;

dqmc_fn_names = {'doubles_fn_u_t_n_h', 'singles_fn_u_t_n_h',...
                'nup_nup_corr_fn_u_t_n_h', 'dd_corr_fn_u_t_n_h',...
                'ns_ns_corr_fn_u_t_n_h', 'nup_fn_u_t_n_h',...
                'ndn_fn_u_t_n_h', 'ndn_ndn_corr_fn_u_t_n_h',...
                'mups_fn_u_t_n_h'};

dqmc_fns = {};
for ii = 1 : length(dqmc_fn_names)
    check_name_fn = @(name) isequal(name, dqmc_fn_names{ii});
    if ~any( cellfun( check_name_fn, fieldnames(thry) ) )
        warning('Theory data did not have a function named %s for use with dqmc fitting function. Set UseFit(ii)=0.', dqmc_fn_names{ii}, ii);
        use_quantity_to_fit(ii) = 0;
        dqmc_fns{ii} = @(P,n) zeros(size(n));
    else
        dqmc_fns{ii} = @(P,n) transpose(squeeze(thry.(dqmc_fn_names{ii})(P(1), P(2), n, P(3))));
    end
end
% rename for ease of use
dd_fn = @(P, n) dqmc_fns{1}(P, n);
ns_fn = @(P, n) dqmc_fns{2}(P, n);
nupcorr_fn = @(P, n) dqmc_fns{3}(P, n);
dd_corr_fn = @(P, n) dqmc_fns{4}(P, n);
ndncorr_fn = @(P, n) nupcorr_fn(P, n); % at some point allow this to be different...
ns_corr_fn = @(P, n) dqmc_fns{5}(P, n);
nup_fn = @(P, n) dqmc_fns{6}(P, n);
ndn_fn = @(P, n) dqmc_fns{7}(P, n);
chem_pot_fn = @(P, n) dqmc_fns{9}(P, n);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use density cutoff
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adjust density cutoffs to maximum and minimum densities in theory data
% set.
if n_cutoffs(1) < thry.n_limits(1)
    n_cutoffs(1) = thry.n_limits(1);
    warning('n_cutoffs lower cutoff was smaller than all DQMC theory data. Increased cutoff');
end

if n_cutoffs(2) > thry.n_limits(2)
    n_cutoffs(2) = thry.n_limits(2);
    warning('n_cutoffs upper cutoff was larger than all DQMC theory data. Increased cutoff');
end

nmin_cutoff = n_cutoffs(1);
nmax_cutoff = n_cutoffs(2);
use_to_fit_index = (nmax_cutoff > n) & (n > nmin_cutoff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define fitting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fixed_params(fixed_params ~= 0) = 1;
lbs = [theory_fns_struct.u_limits(1), theory_fns_struct.t_limits(1),...
       theory_fns_struct.u_limits(1), 0.8, 0.8];
ubs = [theory_fns_struct.u_limits(2), theory_fns_struct.t_limits(2),...
       theory_fns_struct.u_limits(2), 1, 1];

%function for non-int fermi gas
% TODO: correct for imbalance. See fit_fg function
nup_fg_fn = @(P, n) 0.5 * n;
ndn_fg_fn = @(P, n) 0.5 * n;

ns_fg_fn = @(P, n) n - 0.5 * n.^2;
d_fg_fn = @(P, n) (n / 2) .^2;
nup_corr_fg_fn = @(P, n) fg_corr(1 ./ P(2), fg_mu(1./ P(2), 0.5 * n), [0, 1]);
ndn_corr_fg_fn = @(P, n)  nup_corr_fg_fn(P, n);
dd_corr_fg_fn = @(P, n) fg_doubles_corr(1 ./ P(2), fg_mu(1./ P(2), 0.5 * n), fg_mu(1./ P(2), 0.5 * n), [0, 1]);
ns_corr_fg_fn = @(P, n) fg_singles_corr(1 ./ P(2), fg_mu(1./ P(2), 0.5 * n), fg_mu(1./ P(2), 0.5 * n), [0, 1]);
mu_fg_fn = @(P, n) fg_mu(1./ P(2), 0.5 * n);

%define fitting functions for theory fit
% take uncorrected density as argument for these functions.
% exclude pts
n_fit = n(use_to_fit_index);
nerr_fit = nerr(use_to_fit_index);
nup_fit = nup(use_to_fit_index);
nup_err_fit = zero2one(nup_err(use_to_fit_index));
ndn_fit = ndn(use_to_fit_index);
ndn_err_fit = zero2one(ndn_err(use_to_fit_index));
nup_corr_fit = nup_corr(use_to_fit_index);
nup_corr_err_fit = zero2one(nup_corr_err(use_to_fit_index));
ndn_corr_fit = ndn_corr(use_to_fit_index);
ndn_corr_err_fit = zero2one(ndn_corr_err(use_to_fit_index));
ns_fit = ns(use_to_fit_index);
nserr_fit = zero2one(nserr(use_to_fit_index));
ns_corr_fit = ns_corr(use_to_fit_index);
ns_corr_err_fit = zero2one(ns_corr_err(use_to_fit_index));
d_fit = d(use_to_fit_index);
derr_fit = zero2one(derr(use_to_fit_index));
dd_corr_fit = ddcorr(use_to_fit_index);
dd_corr_err_fit = zero2one(ddcorr_err(use_to_fit_index));

% functions
nup_diff_fn = @(P) (nup_fn(P, n_fit / P(4)) * P(4) - nup_fit) ./ nup_err_fit;
ndn_diff_fn = @(P) (ndn_fn(P, n_fit / P(4)) * P(4) - ndn_fit) ./ ndn_err_fit;
nup_corr_diff_fn = @(P) (nupcorr_fn(P, n_fit / P(4)) * P(4)^2 - nup_corr_fit) ./ nup_corr_err_fit;
ndn_corr_diff_fn = @(P) (ndncorr_fn(P, n_fit / P(4)) * P(4)^2 - ndn_corr_fit) ./ ndn_corr_err_fit;
if ~isempty(singles)
    ns_diff_fn = @(P) (ns_fn(P, n_fit / P(4)) * P(4) - ns_fit) ./ nserr_fit;
    ns_corr_diff_fn = @(P) (ns_corr_fn(P, n_fit / P(4)) * P(4)^2 - ns_corr_fit) ./ ns_corr_err_fit;
else
    ns_diff_fn = @(P) 0;
    ns_corr_diff_fn = @(P) 0;
end
if ~isempty(doubles)
    dd_diff_fn = @(P) (dd_fn(P, n_fit / P(4)) * P(4) * P(5) - d_fit) ./ derr_fit;
    dd_corr_diff_fn = @(P) (dd_corr_fn(P, n_fit / P(4)) * P(4)^2 * P(5)^2 - dd_corr_fit) ./ dd_corr_err_fit;
else
    dd_diff_fn = @(P) 0;
    dd_corr_diff_fn = @(P) 0 ;
end

%Select what to fit. Either singles density, correlator, or both
diff_fn = @(P) [dd_diff_fn(P) * use_quantity_to_fit(1),...
                ns_diff_fn(P) * use_quantity_to_fit(2),...
                nup_corr_diff_fn(P) * use_quantity_to_fit(3),...
                ndn_corr_diff_fn(P) * use_quantity_to_fit(3),...
                dd_corr_diff_fn(P) * use_quantity_to_fit(4),...
                ns_corr_diff_fn(P) * use_quantity_to_fit(5),...
                nup_diff_fn(P) * use_quantity_to_fit(6),...
                ndn_diff_fn(P) * use_quantity_to_fit(6)];

% check for density points that are too large for initial fit parameters
if any(isnan(nup_diff_fn(init_params)))
    nan_index = find( isnan(nup_diff_fn(init_params)) );
    for ii = 1 : length(nan_index)
        warning('Experimental density %0.2f corresponding to corrected density %0.2f appears to be outside density range of theory data',...
                n(nan_index(ii)), n(nan_index(ii)) / init_params(4));
    end
end
            
% do dqmc fitting
[fit_params, std_errs, redchi_sqr] = ...
    lsq_fixedp(diff_fn, init_params, fixed_params, lbs, ubs);

% get paramers at maximum and minimum temp and interaction
temp_upper_fp = fit_params;
temp_upper_fp(2) = temp_upper_fp(2) + std_errs(2);
temp_lower_fp = fit_params;
temp_lower_fp(2) = temp_lower_fp(2) - std_errs(2);
%
u_upper_fp = fit_params;
u_upper_fp(1) = u_upper_fp(1) + std_errs(1);
u_lower_fp = fit_params;
u_lower_fp(1) = u_lower_fp(1) - std_errs(1);

% fit chemical potential to harmonic profile
chem_pot_fn_fidelity_corrected = @(P, n) chem_pot_fn(P, n / P(4));

mus = chem_pot_fn_fidelity_corrected(fit_params, n_fit);
mus_unc = 0.5 * (chem_pot_fn_fidelity_corrected(fit_params, n_fit + nerr_fit)...
               - chem_pot_fn_fidelity_corrected(fit_params, n_fit - nerr_fit));

radius = radius_full(use_to_fit_index);
xx = xx_full(use_to_fit_index);
yy = yy_full(use_to_fit_index);
% mu_struct = fit_harmonic_trap_1d(mus, mus_unc, radius, xx, yy);
if use_az_avg
    mu_struct = fit_harmonic_trap_1d(mus, mus_unc, radius, xx, yy);
else
    % How to do this right?
    mus_full = nan(size(xx_full));
    mus_unc_full = nan(size(xx_full));

    mus_full(use_to_fit_index) = mus;
    mus_unc_full(use_to_fit_index) = mus_unc;
    
    mu_struct = fit_harmonic_trap_2d(mus_full, mus_unc_full, radius_full, xx_full, yy_full);
end
          
% store results in struct
fit_struct = struct();

fit_struct.n_cutoffs = n_cutoffs;

fit_struct.xx = xx_full;
fit_struct.yy = yy_full;
fit_struct.radius = radius_full;

% experimental density and correlators corrected for efficiency
fit_struct.n = n / fit_params(4);
fit_struct.n_err = nerr / fit_params(4);
fit_struct.m = m / fit_params(4);
fit_struct.m_err = merr / fit_params(4);

fit_struct.nup = nup / fit_params(4);
fit_struct.nup_err = nup_err / fit_params(4);
fit_struct.nup_corr = nup_corr / fit_params(4)^2;
fit_struct.nup_corr_err = nup_corr_err / fit_params(4)^2;

fit_struct.ndn = ndn / fit_params(4);
fit_struct.ndn_err = ndn_err / fit_params(4);
fit_struct.ndn_corr = ndn_corr / fit_params(4)^2;
fit_struct.ndn_corr_err = ndn_corr_err / fit_params(4)^2;

fit_struct.ns = ns / fit_params(4);
fit_struct.ns_err = nserr / fit_params(4);
fit_struct.ns_corr = ns_corr / fit_params(4)^2;
fit_struct.ns_corr_err = ns_corr_err / fit_params(4)^2;

fit_struct.nd = d / fit_params(4) / fit_params(5);
fit_struct.nd_err = derr / fit_params(4) / fit_params(5);
fit_struct.nd_corr = ddcorr / fit_params(4)^2 / fit_params(5)^2;
fit_struct.nd_corr_err = ddcorr_err / fit_params(4)^2 / fit_params(5)^2;

% fit results
fit_struct.fixed_params = fixed_params;
fit_struct.use_to_fit_indices = use_to_fit_index;
fit_struct.lbs = lbs;
fit_struct.ubs = ubs;
fit_struct.use_param_to_fit = use_quantity_to_fit;
fit_struct.fit_params = fit_params;
fit_struct.std_errs = std_errs;
fit_struct.redchi_sqr = redchi_sqr;

% store sampled fit functions
fit_struct.n_interp = linspace(0, min([max(fit_struct.n) * 1.2, 2]), 300);

fit_struct.nup_fit = nup_fn(fit_params, fit_struct.n_interp);
fit_struct.ndn_fit = ndn_fn(fit_params, fit_struct.n_interp);

fit_struct.ns_fit = ns_fn(fit_params, fit_struct.n_interp);
fit_struct.ns_fit_maxT = ns_fn(temp_upper_fp, fit_struct.n_interp);
fit_struct.ns_fit_minT = ns_fn(temp_lower_fp, fit_struct.n_interp);
fit_struct.ns_fit_maxU = ns_fn(u_upper_fp, fit_struct.n_interp);
fit_struct.ns_fit_minU = ns_fn(u_lower_fp, fit_struct.n_interp);

fit_struct.nd_fit = dd_fn(fit_params, fit_struct.n_interp);
fit_struct.nd_fit_maxT = dd_fn(temp_upper_fp, fit_struct.n_interp);
fit_struct.nd_fit_minT = dd_fn(temp_lower_fp, fit_struct.n_interp);
fit_struct.nd_fit_maxU = dd_fn(u_upper_fp, fit_struct.n_interp);
fit_struct.nd_fit_minU = dd_fn(u_lower_fp, fit_struct.n_interp);

fit_struct.nup_corr_fit = nupcorr_fn(fit_params, fit_struct.n_interp);
fit_struct.nup_corr_fit_maxT = nupcorr_fn(temp_upper_fp, fit_struct.n_interp);
fit_struct.nup_corr_fit_minT = nupcorr_fn(temp_lower_fp, fit_struct.n_interp);
fit_struct.nup_corr_fit_maxU = nupcorr_fn(u_upper_fp, fit_struct.n_interp);
fit_struct.nup_corr_fit_minU = nupcorr_fn(u_lower_fp, fit_struct.n_interp);

fit_struct.ndn_corr_fit = ndncorr_fn(fit_params, fit_struct.n_interp);
fit_struct.ndn_corr_fit_maxT = ndncorr_fn(temp_upper_fp, fit_struct.n_interp);
fit_struct.ndn_corr_fit_minT = ndncorr_fn(temp_lower_fp, fit_struct.n_interp);
fit_struct.ndn_corr_fit_maxU = ndncorr_fn(u_upper_fp, fit_struct.n_interp);
fit_struct.ndn_corr_fit_minU = ndncorr_fn(u_lower_fp, fit_struct.n_interp);

fit_struct.ns_corr_fit = ns_corr_fn(fit_params, fit_struct.n_interp);
fit_struct.ns_corr_fit_maxT = ns_corr_fn(temp_upper_fp, fit_struct.n_interp);
fit_struct.ns_corr_fit_minT = ns_corr_fn(temp_lower_fp, fit_struct.n_interp);
fit_struct.ns_corr_fit_maxU = ns_corr_fn(u_upper_fp, fit_struct.n_interp);
fit_struct.ns_corr_fit_minU = ns_corr_fn(u_lower_fp, fit_struct.n_interp);

fit_struct.nd_corr_fit = dd_corr_fn(fit_params, fit_struct.n_interp);
fit_struct.nd_corr_fit_maxT = dd_corr_fn(temp_upper_fp, fit_struct.n_interp);
fit_struct.nd_corr_fit_minT = dd_corr_fn(temp_lower_fp, fit_struct.n_interp);
fit_struct.nd_corr_fit_maxU = dd_corr_fn(u_upper_fp, fit_struct.n_interp);
fit_struct.nd_corr_fit_minU = dd_corr_fn(u_lower_fp, fit_struct.n_interp);

% store sampled fg functions
fit_struct.nup_fg = nup_fg_fn(fit_params, fit_struct.n_interp);
fit_struct.ndn_fg = ndn_fg_fn(fit_params, fit_struct.n_interp);
fit_struct.ns_fg = ns_fg_fn(fit_params, fit_struct.n_interp);
fit_struct.nd_fg = d_fg_fn(fit_params, fit_struct.n_interp);

fit_struct.nup_corr_fg = nup_corr_fg_fn(fit_params, fit_struct.n_interp);
fit_struct.ndn_corr_fg = ndn_corr_fg_fn(fit_params, fit_struct.n_interp);
fit_struct.ns_corr_fg = ns_corr_fg_fn(fit_params, fit_struct.n_interp);
fit_struct.nd_corr_fg = dd_corr_fg_fn(fit_params, fit_struct.n_interp);

end

function x = zero2one(x)
    x(x == 0) = 1;
end

function [dens, dens_err, dens_corr, dens_corr_err] = get_azavg_quantities_from_DataFolder(df)
    % extract density, density correlations, and uncertainties from an
    % instance of DataFolder class.
    if isempty(df)
        dens = 0;
        dens_err = 0;
        dens_corr = 0;
        dens_corr_err = 0;
        return;
    end
    
    dens = transpose(df.Occs_AzAvg);
    
    dens_err = transpose(df.Occs_AzAvgUnc);
    
    % dens
    dens_corr = transpose(squeeze(...
                    df.Density_Corr_AzAvg(df.CenterIndex_CorrMatrix,...
                                          df.CenterIndex_CorrMatrix + 1, :) ));
                                      
    dens_corr_err = transpose(squeeze(...
                    df.Density_Corr_AzAvgUnc(df.CenterIndex_CorrMatrix,...
                                             df.CenterIndex_CorrMatrix + 1, :) ));
end

function [dens, dens_err, dens_corr, dens_corr_err] = get_2dbin_quantities_from_DataFolder(df, bin_size)
    if isempty(df)
        dens = 0;
        dens_err = 0;
        dens_corr = 0;
        dens_corr_err = 0;
        return;
    end
%     [dens, dens_err, ~, ~] = binImg(df.Occs_ImgAvg, bin_size, bin_size, 'Normal', df.Occs_ImgAvgUnc);
    [dens, ~, ~, dens_err] = binImg(df.Occs_ImgAvg, bin_size, bin_size, 'Normal', df.Occs_ImgAvgUnc);
    dens = transpose(dens(:));
    dens_err = transpose(dens_err(:));

    [dens_corr, ~, ~, dens_corr_err] = binImg( squeeze( ...
                               df.nnc2D(:, :, df.CenterIndex_CorrMatrix,...
                                              df.CenterIndex_CorrMatrix + 1) ),...
                               bin_size, bin_size, 'Normal',...
                               squeeze(...
                               df.nnc2Dunc(:, :, df.CenterIndex_CorrMatrix,...
                                                  df.CenterIndex_CorrMatrix + 1) ) );
    dens_corr = transpose(dens_corr(:));    
    dens_corr_err = transpose(dens_corr_err(:));
                                         
end