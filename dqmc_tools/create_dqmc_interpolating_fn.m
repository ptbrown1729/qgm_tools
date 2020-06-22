function [interp_fn_struct] = create_dqmc_interpolating_fn(dqmc_results_grid)
% dqmc_results_grid is produced by...
% Assume you have ordered everything (Ts, Us, Mu_ups, Mu_dns)
% currently a couple of problematic use cases
% (1) when we only have a single point along some dimensions, e.g.
% interaction. Interpolation doesn't like this very much. Think the
% solution is to add extra point extremely close along this dimension. E.g.
% if we have only T = 0.3, then add T = 0.3 +/- 1e-6 ... corrected this in
% convert_dqmc_to_grid_v2
% (2) when combining data sets where one data set has information along all
% dimensions, but another doesn't. For example, the first data set is a
% large set with varying U, T, Mu, and imbalance. The second dataset has
% the same variable U, T, and Mu but is balanced. ... correct this with
% switching in here...

is_balanced = 0;
if ndims(dqmc_results_grid.Mus_Up_Grid) == 3
    is_balanced = 1;
end

interp_fn_struct = struct();

interp_fn_struct.u_limits = [min(dqmc_results_grid.Us_Grid(:)),...
                             max(dqmc_results_grid.Us_Grid(:))];
interp_fn_struct.mu_up_limits = [min(dqmc_results_grid.Mus_Up_Grid(:)),...
                                 max(dqmc_results_grid.Mus_Up_Grid(:))];
interp_fn_struct.mu_dn_limits = [min(dqmc_results_grid.Mus_Dn_Grid(:)),...
                                 max(dqmc_results_grid.Mus_Dn_Grid(:))];
interp_fn_struct.t_limits = [min(dqmc_results_grid.Ts_Grid(:)),...
                             max(dqmc_results_grid.Ts_Grid(:))];

%ns and n as function of (Mu, T, U)
% here Mus are x-like = (dim 2 like), Ts are y-like = (dim 1 like),
% and Us are z-like = (dim 3 like)
% Order of function arguments is like (x, y, z)
% Order of interpn arguments is like (dim1, dim2, dim3) = (y, x, z)

if ~is_balanced
    % keep function ordered as (x, y, z)
    % but interp is ordered (y, x, z)
    interp_struct_fn = @(U, T, mu_up, mu_dn, quantity) interpn(...
            dqmc_results_grid.Ts_Grid,...
            dqmc_results_grid.Us_Grid,...
            dqmc_results_grid.Mus_Up_Grid,...
            dqmc_results_grid.Mus_Dn_Grid,...
            squeeze(quantity(:, :, :, :, 1)),...
            T, U, mu_up, mu_dn);
else
    interp_struct_fn = @(U, T, mu_up, mu_dn, quantity) interpn(...
        dqmc_results_grid.Ts_Grid,...
        dqmc_results_grid.Us_Grid,...
        dqmc_results_grid.Mus_Up_Grid,...
        squeeze(quantity(:, :, :, 1)),...
        T, U, mu_up);
end
    
% densities
interp_fn_struct.singles_u_t_mu = @(U, T, mu_up, mu_dn) interp_struct_fn(...
    U, T, mu_up, mu_dn, dqmc_results_grid.SingDensity);
    
interp_fn_struct.density_u_t_mu = @(U, T, mu_up, mu_dn) interp_struct_fn(...
    U, T, mu_up, mu_dn, dqmc_results_grid.Density);

interp_fn_struct.doubles_u_t_mu = @(U, T, mu_up, mu_dn) interp_struct_fn(...
    U, T, mu_up, mu_dn, dqmc_results_grid.DoubDensity);

% correlators
interp_fn_struct.nup_nup_corr_u_t_mu = @(U, T, mu_up, mu_dn) interp_struct_fn(...
    U, T, mu_up, mu_dn, dqmc_results_grid.NNupup);

interp_fn_struct.ndn_ndn_corr_u_t_mu = @(U, T, mu_up, mu_dn) interp_struct_fn(...
    U, T, mu_up, mu_dn, dqmc_results_grid.NNupup);

interp_fn_struct.nup_ndn_corr_u_t_mu = @(U, T, mu_up, mu_dn) interp_struct_fn(...
    U, T, mu_up, mu_dn, dqmc_results_grid.NNupdn);

interp_fn_struct.double_double_corr_u_t_mu = @(U, T, mu_up, mu_dn) interp_struct_fn(...
    U, T, mu_up, mu_dn, dqmc_results_grid.NNdd);

interp_fn_struct.sz_sz_corr_u_t_mu = @(U, T, mu_up, mu_dn) interp_struct_fn(...
    U, T, mu_up, mu_dn, dqmc_results_grid.NNSzSz);

interp_fn_struct.ns_ns_corr_u_t_mu = @(U, T, mu_up, mu_dn) interp_struct_fn(...
    U, T, mu_up, mu_dn, dqmc_results_grid.NNss);

% susceptibilities
interp_fn_struct.magnetic_susc = @(U, T, mu_up, mu_dn) interp_struct_fn(...
    U, T, mu_up, mu_dn, dqmc_results_grid.magnetic_susc);

interp_fn_struct.nup_susc = @(U, T, mu_up, mu_dn) interp_struct_fn(...
    U, T, mu_up, mu_dn, dqmc_results_grid.nup_susc);

% sign
interp_fn_struct.sign_u_t_mu = @(U, T, mu_up, mu_dn) interp_struct_fn(...
    U, T, mu_up, mu_dn, dqmc_results_grid.AvgSign);

% interp_fn_struct.singles_mu_t_u = @(U, T, mu_up, mu_dn) interpn(...
%     dqmc_results_grid.Ts_Grid,...
%     dqmc_results_grid.Us_Grid,...
%     dqmc_results_grid.Mus_Up_Grid,...
%     dqmc_results_grid.Mus_Dn_Grid,...
%     dqmc_results_grid.SingDensity(:, :, :, :, 1),...
%     T, U, mu_up, mu_dn);

% we want to use the variables (T, U, n, h) instead of (T, U, mu_up, mu_dn)
% here h = 0.5*(mu_up - mu_dn)
% produce a regularly grid in this new space
% n_limits = [min(dqmc_results_grid.Density(:)),...
%             max(dqmc_results_grid.Density(:))];
if is_balanced
    n_limits = [min(min(min(dqmc_results_grid.Density(:, :, :, 1)))),...
              max(max(max(dqmc_results_grid.Density(:, :, :, 1))))];
else 
    n_limits = [min(min(min(min(dqmc_results_grid.Density(:, :, :, :, 1))))),...
                max(max(max(max(dqmc_results_grid.Density(:, :, :, :, 1)))))];
end

interp_fn_struct.n_limits = n_limits;       
h_limits = 0.5 * [min(dqmc_results_grid.Mus_Up_Grid(:) - dqmc_results_grid.Mus_Dn_Grid(:)),...
                  max(dqmc_results_grid.Mus_Up_Grid(:) - dqmc_results_grid.Mus_Dn_Grid(:))];
    
npts = linspace(n_limits(1), n_limits(2), 30);
hpts = 0.5 * sort(unique(squeeze(dqmc_results_grid.Mus_Up_Grid(1, 1, :, :) ...
                               - dqmc_results_grid.Mus_Dn_Grid(1, 1, :, :))));
% hpts = linspace(h_limits(1), h_limits(2), 31);
% ensure that temperatures and interactions are the same as before. This
% way we only need to do interpolation along the last two dimensions.
if ~is_balanced
    [ts_nhgrid, us_nhgrid, n_nhgrid, h_nhgrid] = ndgrid(...
                                            sort(unique(dqmc_results_grid.Ts_Grid)),...
                                            sort(unique(dqmc_results_grid.Us_Grid)),...
                                            npts, hpts);
else
    [ts_nhgrid, us_nhgrid, n_nhgrid] = ndgrid(...
                                        sort(unique(dqmc_results_grid.Ts_Grid)),...
                                        sort(unique(dqmc_results_grid.Us_Grid)),...
                                        npts);
end

%get nsingles on at each point on our grid by interpolating
% unlike ndgrid and interpn, griddata expects order to be
% (x, y, z) = (dim 1, dim 2, dim 3)
% interpolate in each dimension separately. This seems like a better
% strategy than trying to do a 4D interpolation. And I don't see the
% advantage to that, as the Ts and Us are the same as before the coordinate
% transformation
muup_t_u_n_h = zeros(size(ts_nhgrid));
mudn_t_u_n_h = zeros(size(ts_nhgrid));
nup_t_u_n_h = zeros(size(ts_nhgrid));
ndn_t_u_n_h = zeros(size(ts_nhgrid));
singles_t_u_n_h = zeros(size(ts_nhgrid));
doubles_t_u_n_h = zeros(size(ts_nhgrid));
% correlators
nup_nup_corr_t_u_n_h = zeros(size(ts_nhgrid));
ndn_ndn_corr_t_u_n_h = zeros(size(ts_nhgrid));
nup_ndn_corr_t_u_n_h = zeros(size(ts_nhgrid));
double_double_corr_t_u_n_h = zeros(size(ts_nhgrid));
sz_sz_corr_t_u_n_h = zeros(size(ts_nhgrid));
ns_ns_corr_t_u_n_h = zeros(size(ts_nhgrid));
% susceptibilities
magnetic_susc_t_u_n_h = zeros(size(ts_nhgrid));
nup_susc_t_u_n_h =  zeros(size(ts_nhgrid));
% sign
sign_t_u_n_h = zeros(size(ts_nhgrid));

if ~is_balanced
    griddata_fn = @(ii, jj, quantity) griddata(...
                0.5 * squeeze(dqmc_results_grid.Mus_Up_Grid(ii, jj, :, :) ...
                            - dqmc_results_grid.Mus_Dn_Grid(ii, jj, :, :)),...
                squeeze(dqmc_results_grid.Density(ii, jj, :, :, 1)),...
                squeeze(quantity(ii, jj, :, :, 1)),...
                squeeze(h_nhgrid(ii, jj, :, :)),...
                squeeze(n_nhgrid(ii, jj, :, :)));
            
    for ii = 1:size(ts_nhgrid, 1)
        for jj = 1:size(ts_nhgrid, 2)
            % unlike ndgrid and interpn, griddata expects order to be
            % (x, y, z) = (dim 1, dim 2, dim 3)
            % so when squeezing out the first two dimensions, that means it
            % will want (dim 4, dim 3)
            % mus
            muup_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.Mus_Up_Grid);
            mudn_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.Mus_Dn_Grid);
            % densities
            singles_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.SingDensity);
            nup_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.SpinUpOcc);
            ndn_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.SpinDownOcc);
            doubles_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.DoubDensity);
            % correlators
            nup_nup_corr_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.NNupup);
            ndn_ndn_corr_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.NNdndn);
            nup_ndn_corr_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.NNupdn);
            double_double_corr_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.NNdd);
            sz_sz_corr_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.NNSzSz);
            ns_ns_corr_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.NNss);
            % susceptibilities
            magnetic_susc_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.magnetic_susc);
            nup_susc_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.nup_susc);
            % sign
            sign_t_u_n_h(ii, jj, :, :) = griddata_fn(ii, jj, dqmc_results_grid.avg_sign);
        end
    end
    
   interp_u_t_n_h_fn = @(U, T, n, h, quantity)interpn(...
                               ts_nhgrid, us_nhgrid, n_nhgrid, h_nhgrid,...
                               quantity, T, U, n, h);
else
    %get nsingles on at each point on our grid by interpolating
    % unlike ndgrid and interpn, griddata expects order to be
    % (x, y, z) = (dim 1, dim 2, dime 3)
    griddata_fn = @(quantity) griddata(...
        dqmc_results_grid.Us_Grid,...
        dqmc_results_grid.Ts_Grid,...
        dqmc_results_grid.Density(:, :, :, 1),...
        quantity(:, :, :, 1),...
        us_nhgrid, ts_nhgrid, n_nhgrid);
    
    % mus
    muup_t_u_n_h = griddata_fn(dqmc_results_grid.Mus_Up_Grid);
    mudn_t_u_n_h = griddata_fn(dqmc_results_grid.Mus_Dn_Grid);
    % densities
    singles_t_u_n_h = griddata_fn(dqmc_results_grid.SingDensity);
    nup_t_u_n_h = griddata_fn(dqmc_results_grid.SpinUpOcc);
    ndn_t_u_n_h = griddata_fn(dqmc_results_grid.SpinDownOcc);
    doubles_t_u_n_h = griddata_fn(dqmc_results_grid.DoubDensity);
    % correlators
    nup_nup_corr_t_u_n_h = griddata_fn(dqmc_results_grid.NNupup);
    ndn_ndn_corr_t_u_n_h = griddata_fn(dqmc_results_grid.NNdndn);
    nup_ndn_corr_t_u_n_h = griddata_fn(dqmc_results_grid.NNupdn);
    double_double_corr_t_u_n_h = griddata_fn(dqmc_results_grid.NNdd);
    sz_sz_corr_t_u_n_h = griddata_fn(dqmc_results_grid.NNSzSz);
    ns_ns_corr_t_u_n_h = griddata_fn(dqmc_results_grid.NNss);
    % susceptibilities
    magnetic_susc_t_u_n_h = griddata_fn(dqmc_results_grid.magnetic_susc);
    nup_susc_t_u_n_h = griddata_fn(dqmc_results_grid.nup_susc);
    % sign
    sign_t_u_n_h = griddata_fn(dqmc_results_grid.avg_sign);
    
    interp_u_t_n_h_fn = @(U, T, n, h, quantity)interpn(...
                                ts_nhgrid, us_nhgrid, n_nhgrid,...
                                quantity, T, U, n);
end


%interpolate on the regular space grid of (U, T, n, h) points
% interp_fn_struct.mups_fn_u_t_n_h = @(U, T, n, h) interpn(ts_nhgrid, us_nhgrid, n_nhgrid, h_nhgrid,...
%                                                        muup_t_u_n_h,...
%                                                        T, U, n, h);

% mus
interp_fn_struct.mups_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, muup_t_u_n_h);
interp_fn_struct.mdns_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, mudn_t_u_n_h);
% densities
interp_fn_struct.nup_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, nup_t_u_n_h);
interp_fn_struct.ndn_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, ndn_t_u_n_h);
interp_fn_struct.singles_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, singles_t_u_n_h);
interp_fn_struct.doubles_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, doubles_t_u_n_h);
% correlators
interp_fn_struct.nup_nup_corr_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, nup_nup_corr_t_u_n_h);
interp_fn_struct.ndn_ndn_corr_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, ndn_ndn_corr_t_u_n_h);
interp_fn_struct.nup_ndn_corr_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, nup_ndn_corr_t_u_n_h);
interp_fn_struct.dd_corr_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, double_double_corr_t_u_n_h);
interp_fn_struct.sz_sz_corr_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, sz_sz_corr_t_u_n_h);
interp_fn_struct.ns_ns_corr_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, ns_ns_corr_t_u_n_h);
% susceptibilities
interp_fn_struct.magnetic_susc_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, magnetic_susc_t_u_n_h);
interp_fn_struct.nup_susc_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, nup_susc_t_u_n_h);
% sign
interp_fn_struct.avgsign_fn_u_t_n_h = @(U, T, n, h) interp_u_t_n_h_fn(U, T, n, h, sign_t_u_n_h);

end