function [fit_params, std_errs, redchi_sqr, mu_struct, fit_struct] = fit_fg(...
    ups, downs, singles, doubles, use_quantity_to_fit, n_cutoffs,...
    init_params, fixed_params, use_az_avg, img_bin_size)
% Fit singles density and single-spin component correlators versus density
% against temperature.
%
% UseToFit = [1, 1, 1, 1, 1, 1] selects which of the doubles density, singles
% density, density correlator, doubles correlator, singles correlator, and 
% individual component densities to use in fitting. In that order
%
% n_cutoffs = [nmin, nmax] limit the densities that are included in the fit
%
% InitParams = [T, h, imaging fidelity, doublon fidelity]. Note that the
% doublon fidelity is only the transfer efficiency of the doublon process.
% The imaging fidelity is also taken into account for the doublons.
% h = 0.5 * (mu_up - mu_dn)
%
% fixed_params specifies which of the initial parameters are to be fixed. It
% should be the same size as InitParams, with a zero for any parameter that
% is allowed to vary, and a 1 for any parameter that is not allowed to vary.
% If a one is given, the parameter will be fixed at the value provided in
% InitParams.
%
% theory_fn_structure 

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

if length(use_quantity_to_fit) ~= 6
    error('use_quantity_to_fit should be 6 elements long');
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

if ~exist('fixed_params', 'var') || isempty(fixed_params)
    fixed_params = zeros(size(init_params));
end

if ~isequal( size(init_params), size(fixed_params) )
    error('init_params and fixed_params must be the same size');
end

% account for possibility we may not have all of ones, singles, doubles,
% etc.
if isempty(downs) && ~isempty(ups)
    downs = ups;
elseif isempty(ups) && ~isempty(downs)
    ups = downs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract densities and correlators from DataFolder instances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if use_az_avg
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

    % radial position
    radius_full = ups.RadialPos;
    
    [xx_full, yy_full] = meshgrid( 1 : ups.ImgCropSize, 1 : ups.ImgCropSize ); 
    xx_full = xx_full - ups.Cx_AzAvg;
    yy_full = yy_full - ups.Cy_AzAvg;

else
    % Get information from 2D cloud
    % this mode is a work in progress.
%     img_bin_size = 6;
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%use density cutoff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmin_cutoff = n_cutoffs(1);
nmax_cutoff = n_cutoffs(2);
use_to_fit_index = (nmax_cutoff > n) & ( n > nmin_cutoff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%select numerical method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_avg_fn = @(P,n) squeeze( fg_mu_avg( 1 / P(1), P(2), n));

mu_up_fn = @(P,n) mu_avg_fn(P,n) + P(2);
mu_dn_fn = @(P,n) mu_avg_fn(P,n) - P(2);

nup_fn = @(P,n) squeeze( fg_density( 1 / P(1), mu_up_fn(P, n) ) );
ndn_fn = @(P,n) squeeze( fg_density( 1 / P(1), mu_dn_fn(P, n) ) );

ns_fn = @(P, n) squeeze( fg_singles( 1 / P(1),...
                         mu_up_fn(P, n),...
                         mu_dn_fn(P, n)) );
                     
ns_corr_fn = @(P, n) squeeze( fg_singles_corr( 1 / P(1),...
                              mu_up_fn(P, n),...
                              mu_dn_fn(P, n),...
                              [0, 1]) );
                          
nupcorr_fn = @(P, n) squeeze( fg_corr( 1 / P(1),...
                              mu_up_fn(P, n),...
                              [0, 1]) );
                          
ndncorr_fn = @(P, n) squeeze( fg_corr( 1 / P(1),...
                              mu_dn_fn(P, n),...
                              [0, 1]) );

dd_fn = @(P, n) squeeze( fg_doubles( 1 / P(1),...
                         mu_up_fn(P, n),...
                         mu_dn_fn(P, n)) );

dd_corr_fn = @(P, n) squeeze( fg_doubles_corr( 1 / P(1),...
                              mu_up_fn(P, n),...
                              mu_dn_fn(P, n),...
                              [0, 1]) ); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define fitting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fixed_params(fixed_params ~= 0) = 1;
lbs = [0, -10, 0.9, 0.9];
ubs = [10, 10, 1, 1];

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

% up and dn densities
nup_diff_fn = @(P) (nup_fn(P, n_fit / P(3)) * P(3) - nup_fit) ./ nup_err_fit;
ndn_diff_fn = @(P) (ndn_fn(P, n_fit / P(3)) * P(3) - ndn_fit) ./ ndn_err_fit;
nup_corr_diff_fn = @(P) (nupcorr_fn(P, n_fit / P(3)) * P(3)^2 - nup_corr_fit) ./ nup_corr_err_fit;
ndn_corr_diff_fn = @(P) (ndncorr_fn(P, n_fit / P(3)) * P(3)^2 - ndn_corr_fit) ./ ndn_corr_err_fit;

% singles correlators
if ~isempty(singles)
    ns_diff_fn = @(P) (ns_fn(P, n_fit / P(3)) * P(3) - ns_fit) ./ nserr_fit;
    ns_corr_diff_fn = @(P) (ns_corr_fn(P, n_fit / P(3)) - ns_corr_fit) ./ ns_corr_err_fit;
else
    ns_diff_fn = @(P) 0;
    ns_corr_diff_fn = @(P) 0;
end

% doublon correlators
dd_fit_fn = @(P,n) dd_fn(P, n / P(3)) * P(3) * P(4);
dd_corr_fit_fn = @(P,n) dd_corr_fn(P, n / P(3)) * P(3)^2 * P(4)^2;
if ~isempty(doubles)
    dd_diff_fn = @(P) (dd_fn(P, n_fit / P(3)) * P(3) * P(4) - d_fit) ./ derr_fit;
    dd_corr_diff_fn = @(P) (dd_corr_fit_fn(P, n_fit / P(3)) * P(3)^2 * P(4)^2 - dd_corr_fit) ./ dd_corr_err_fit;
else
    dd_diff_fn = @(P) 0;
    dd_corr_diff_fn = @(P) 0 ;
end

%Select what to fit. Either singles density, correlator, or
%both
diff_fn = @(P) [dd_diff_fn(P) * use_quantity_to_fit(1),...
                ns_diff_fn(P) * use_quantity_to_fit(2),...
                nup_corr_diff_fn(P) * use_quantity_to_fit(3),...
                ndn_corr_diff_fn(P) * use_quantity_to_fit(3),...
                dd_corr_diff_fn(P) * use_quantity_to_fit(4),...
                ns_corr_diff_fn(P) * use_quantity_to_fit(5),...
                nup_diff_fn(P) * use_quantity_to_fit(6),...
                ndn_diff_fn(P) * use_quantity_to_fit(6)];

% do theory fitting
[fit_params, std_errs, redchi_sqr] = ...
    lsq_fixedp(diff_fn, init_params, fixed_params, lbs, ubs);

% fit chemical potential to harmonic profile
chem_pot_fn_fidelity_corrected = @(P, n) mu_avg_fn(P, n / P(3));

mus = chem_pot_fn_fidelity_corrected(fit_params, n_fit);
mus_unc = 0.5 * (chem_pot_fn_fidelity_corrected(fit_params, n_fit + nerr_fit)...
               - chem_pot_fn_fidelity_corrected(fit_params, n_fit - nerr_fit) );

if use_az_avg
    mu_struct = fit_harmonic_trap_1d(mus, mus_unc, radius_full(use_to_fit_index), xx_full(use_to_fit_index), yy_full(use_to_fit_index));
else
    % How to do this right?
    mus_full = nan(size(xx_full));
    mus_unc_full = nan(size(xx_full));

    mus_full(use_to_fit_index) = mus;
    mus_unc_full(use_to_fit_index) = mus_unc;
    
    mu_struct = fit_harmonic_trap_2d(mus_full, mus_unc_full, radius_full, xx_full, yy_full);
end

% get paramers at maximum and minimum temp and interaction
temp_upper_fp = fit_params;
temp_upper_fp(1) = temp_upper_fp(1) + std_errs(1);
temp_lower_fp = fit_params;
temp_lower_fp(1) = temp_lower_fp(1) - std_errs(1);

% store results in struct
%%% store info in struct
fit_struct = struct();

fit_struct.n_cutoffs = n_cutoffs;

fit_struct.xx = xx_full;
fit_struct.yy = yy_full;
fit_struct.radius = radius_full;

fit_struct.n = n / fit_params(3);
fit_struct.n_err = nerr / fit_params(3);
fit_struct.m = m / fit_params(3);
fit_struct.m_err = merr / fit_params(3);

% component densities
fit_struct.nup = nup / fit_params(3);
fit_struct.nup_err = nup_err / fit_params(3);
fit_struct.nup_corr = nup_corr / fit_params(3)^2;
fit_struct.nup_corr_err = nup_corr_err / fit_params(3)^2;

fit_struct.ndn = ndn / fit_params(3);
fit_struct.ndn_err = ndn_err / fit_params(3);
fit_struct.ndn_corr = ndn_corr / fit_params(3)^2;
fit_struct.ndn_corr_err = ndn_corr_err / fit_params(3)^2;

fit_struct.ns = ns / fit_params(3);
fit_struct.ns_err = nserr / fit_params(3);
fit_struct.ns_corr = ns_corr / fit_params(3)^2;
fit_struct.ns_corr_err = ns_corr_err / fit_params(3)^2;

fit_struct.nd = d / fit_params(3) / fit_params(4);
fit_struct.nd_err = derr / fit_params(3) / fit_params(4);
fit_struct.nd_corr = ddcorr / fit_params(3)^2  / fit_params(4)^2;
fit_struct.nd_corr_err = ddcorr_err / fit_params(3)^2 / fit_params(4)^2;

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

fit_struct.nd_fit = dd_fn(fit_params, fit_struct.n_interp);
fit_struct.nd_fit_maxT = dd_fn(temp_upper_fp, fit_struct.n_interp);
fit_struct.nd_fit_minT = dd_fn(temp_lower_fp, fit_struct.n_interp);

fit_struct.nup_corr_fit = nupcorr_fn(fit_params, fit_struct.n_interp);
fit_struct.nup_corr_fit_maxT = nupcorr_fn(temp_upper_fp, fit_struct.n_interp);
fit_struct.nup_corr_fit_minT = nupcorr_fn(temp_lower_fp, fit_struct.n_interp);

fit_struct.ndn_corr_fit = ndncorr_fn(fit_params, fit_struct.n_interp);
fit_struct.ndn_corr_fit_maxT = ndncorr_fn(temp_upper_fp, fit_struct.n_interp);
fit_struct.ndn_corr_fit_minT = ndncorr_fn(temp_lower_fp, fit_struct.n_interp);

fit_struct.ns_corr_fit = ns_corr_fn(fit_params, fit_struct.n_interp);
fit_struct.ns_corr_fit_maxT = ns_corr_fn(temp_upper_fp, fit_struct.n_interp);
fit_struct.ns_corr_fit_minT = ns_corr_fn(temp_lower_fp, fit_struct.n_interp);

fit_struct.nd_corr_fit = dd_corr_fn(fit_params, fit_struct.n_interp);
fit_struct.nd_corr_fit_maxT = dd_corr_fn(temp_upper_fp, fit_struct.n_interp);
fit_struct.nd_corr_fit_minT = dd_corr_fn(temp_lower_fp, fit_struct.n_interp);

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