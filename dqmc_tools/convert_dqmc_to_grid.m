function dqmc_results_grid = convert_dqmc_to_grid(dqmc_results_structs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%organize dqmc results into regular grids of parameters. Order of
%coordinates is (T, U, Mu_Up, Mu_Dn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%If you vary several parameters, such as interaction (U), chemical potential (Mu)
%and temperature (T), there are different ways of organizing the outputs.
% One convenient way is to store these in regular grids, because these can
% then be directly interpolated using interp2 or interp3

% if ~exist('T_lims', 'var') || isempty(T_lims)
%     T_lims = [-inf, inf];
% end
% 
% if ~exist('U_lims', 'var') || isempty(U_lims)
%     U_lims = [-inf, inf];
% end
% 
% if ~exist('Mu_up_lims', 'var') || isempty(Mu_up_lims)
%     Mu_up_lims = [-inf, inf];
% end
% 
% if ~exist('Mu_dn_lims', 'var') || isempty(Mu_dn_lims)
%     Mu_dn_lims = [-inf, inf];
% end
% 
% if ~exist('Delta_mu_lims', 'var') || isempty(Delta_mu_lims)
%     Delta_mu_lims = [-inf, inf];
% end

% combine DQMC results
fields = fieldnames(dqmc_results_structs(1));
dqmc_results = struct();
for ii = 1:length(fields)
    if ischar(dqmc_results_structs(1).(fields{ii}))
        dqmc_results.(fields{ii}) = {dqmc_results_structs.(fields{ii})}; 
    elseif iscell(dqmc_results_structs(1).(fields{ii}))
        dqmc_results.(fields{ii}) = horzcat(dqmc_results_structs.(fields{ii}));
    else
        dqmc_results.(fields{ii}) = vertcat(dqmc_results_structs.(fields{ii}));
    end
end

% restrict these to certain set of parameters we care about
% can add to this filtering so only get values within certain ranges.
% might also think about doing this already when reading the file. Maybe
% this would speed up the process more.

% full_size = length(dqmc_results.T);
% check_limits = 1;
% if check_limits
%     % limits
%     T_min = T_lims(1);
%     T_max = T_lims(2);
%     U_min = U_lims(1);
%     U_max = U_lims(2);
%     mu_ups_min = Mu_up_lims(1);
%     mu_ups_max = Mu_up_lims(2);
%     mu_dns_min = Mu_dn_lims(1);
%     mu_dns_max = Mu_dn_lims(2);
%     mu_deltas_min = Delta_mu_lims(1);
%     mu_deltas_max = Delta_mu_lims(2);
%     
%     % initial values
%     mu_ups_init = dqmc_results.Mu_Up;
%     mu_dns_init = dqmc_results.Mu_Dn;
%     mu_deltas_init = mu_ups_init - mu_dns_init;
%     T_init = dqmc_results.T;
%     U_init = dqmc_results.U;
%     
%     include_logical_mat = T_init >= T_min & T_init <= T_max ...
%                 & U_init >= U_min & U_init <= U_max ...
%                 & mu_ups_init >= mu_ups_min & mu_ups_init <= mu_ups_max ...
%                 & mu_dns_init >= mu_dns_min & mu_dns_init <= mu_dns_max ...
%                 & mu_deltas_init >= mu_deltas_min & mu_deltas_init <= mu_deltas_max;
%     
%     for ii = 1:length(fields)
%         if size(dqmc_results.(fields{ii}), 1) == full_size
%             dim_sizes = size(dqmc_results.(fields{ii}));
%             extra_dim_sizes = dim_sizes(2:end);
% %             dqmc_results.(fields{ii}) = dqmc_results.(fields{ii})(repmat(include_logical_mat, [1, extra_dim_sizes]));
%              % this only keeps the first dimension of whatever field we are
%              % considering. Want to keep all dimensions. One option is to
%              % move this checking to aggregate_dqmc_output_files
%              dqmc_results.(fields{ii}) = dqmc_results.(fields{ii})(include_logical_mat);
%         end
%     end
% end

temps_list = dqmc_results.T;
mus_up_list = dqmc_results.Mu_Up; %0.5 * (dqmc_results.Mu_Up + dqmc_results.Mu_Dn);
mus_dn_list = dqmc_results.Mu_Dn; %0.5*(dqmc_results.Mu_Up - dqmc_results.Mu_Dn);
us_list = dqmc_results.U;
    

% ensure these are column vectors
temps_list = reshape(temps_list, [length(temps_list), 1]);
mus_up_list = reshape(mus_up_list, [length(mus_up_list), 1]);
mus_dn_list = reshape(mus_dn_list, [length(mus_dn_list), 1]);
us_list = reshape(us_list, [length(us_list), 1]);

mus_up_unique = sort(unique(mus_up_list));
mus_dn_unique = sort(unique(mus_dn_list));
us_unique = sort(unique(us_list));
ts_unique = sort(unique(temps_list));

% create grid of all the points we have...
% ndgrid is just like meshgrid, except returns everything in matrix order,
% so first return is Y-like and second return is X-like (opposite of
% meshgrid)
% [us, mus, ts] = ndgrid(us_unique, mus_unique, ts_unique);
if isequal(mus_up_list, mus_dn_list)
    % if there is no imbalance, mus_dn is redundant
    [ts, us, mus_up] = ndgrid(ts_unique, us_unique, mus_up_unique);
    mus_dn = mus_up;
else
    % if there is imbalance, we need an extra index for mus_dn
    [ts, us, mus_up, mus_dn] = ndgrid(ts_unique, us_unique, mus_up_unique, mus_dn_unique);
end
    
% In the case where parameters were not already generated on a meshgrid type
% grid, I don't know how to nicely vectorize this. Is currently slow.
% The parameters in simulation_param_table should be in the opposite order 
% that they appear in ndgrid above.
% simulation_param_table = horzcat(temps_list, mus_list, us_list);
simulation_param_table = horzcat(mus_dn_list, mus_up_list, us_list, temps_list);
indices = nan(size(mus_up));

%For each point in our grid, find what linear index that is associated
%with.
for ii = 1:numel(mus_up)
%     CurrentParams = [ts(ii), mus(ii), us(ii)];
    CurrentParams = [mus_dn(ii), mus_up(ii), us(ii), ts(ii)];
    index = find(ismember(simulation_param_table, CurrentParams, 'rows'));
    if ~isempty(index)
        if length(index)>1
            fprintf('Warning, found multiple simulations at T = %0.2f, U = %0.2f, Mu mean = %0.2f, Mu delta = %0.2f. Only the first one will be used.\n',...
                ts(ii), us(ii), mus_up(ii), mus_dn(ii));
        end
        indices(ii) = index(1);
    end
    fprintf('%d/%d\n', ii,numel(mus_up));
end

num_nans = sum(isnan(indices(:)));
if num_nans > 0 
    fprintf('Found %d nans in Indices. This means there were (T,U,Mu,DeltaMu) points in our grids that were not in the input DQMC file \n',...
        num_nans);
end
indices_no_nans = indices;
indices_no_nans(isnan(indices_no_nans)) = 1;

% if we have only a single U or T, this will cause us problems later when
% we try to interpolate.
if size(ts, 1) == 1
    t_small = ts(1) - 1e-6;
    t_large = ts(1) + 1e-6;
    ts = cat(1, t_small * ones(size(ts)), ts, t_large * ones(size(ts)));
    us = cat(1, us, us, us);
    mus_up = cat(1, mus_up, mus_up, mus_up);
    mus_dn = cat(1, mus_dn, mus_dn, mus_dn);
    indices = cat(1, indices, indices, indices);
    indices_no_nans = cat(1, indices_no_nans, indices_no_nans, indices_no_nans);
end

if size(us, 2) == 1
    u_small = us(1) - 1e-6;
    u_large = us(1) + 1e-6;
    ts = cat(2, ts, ts, ts);
    us = cat(2, u_small * ones(size(us)), us, u_large * ones(size(us)));
    mus_up = cat(2, mus_up, mus_up, mus_up);
    mus_dn = cat(2, mus_dn, mus_dn, mus_dn);
    indices = cat(2, indices, indices, indices);
    indices_no_nans = cat(2, indices_no_nans, indices_no_nans, indices_no_nans);
end

%Loop over all structures in the DQMC file and reshape them to match the
%size of our grids over their final dimensions.
fields = fieldnames(dqmc_results);
dqmc_results_grid = struct;
dqmc_results_grid.Us_Grid = us;
dqmc_results_grid.Mus_Up_Grid = mus_up;
dqmc_results_grid.Mus_Dn_Grid = mus_dn;
dqmc_results_grid.Ts_Grid = ts;

% assumption: the first dimension of each field corresponds to changing the
% input parameters.
n_input_params = length(mus_up_list);
for ii = 1:numel(fields)    
    old_field = dqmc_results.(fields{ii});
    old_size = size(old_field);
    
    if ischar(fields{ii}) && (strcmp(fields{ii}, 'U') || strcmp(fields{ii}, 'T')...
                       || strcmp(fields{ii}, 'Mu_Up') || strcmp(fields{ii}, 'Mu_Dn'))
        continue;
    end
    
    if (size(old_field, 1) ~= n_input_params) || ~isnumeric(old_field)
        % need to exclude general data fields, which don't correspond to
        % individual qmc runs. These are fields like the directory of the
        % qmc data, or the number of qmc output files.
        new_field = old_field;
    else
        
        new_size = [size(mus_up), old_size(2:end)];
        new_field = nan(new_size);

        if ndims(old_field) == 2
            new_field(:) = old_field(indices_no_nans, :);
            inds = repmat(isnan(indices), [ones(1, ndims(mus_up)), size(old_field, 2)]);
            new_field(inds) = 0;
        elseif ndims(old_field) == 3
            new_field(:) = old_field(indices_no_nans, :, :);
            inds = repmat(isnan(indices), [ones(1, ndims(mus_up)), size(old_field, 2), size(old_field, 3)]);
            new_field(inds) = 0;
        elseif ndims(old_field) == 4
            new_field(:) = old_field(indices_no_nans, :,:,:);
            inds = repmat(isnan(indices), [ones(1, ndims(mus_up)), size(old_field, 2), size(old_field, 3), size(old_field, 4)]);
            new_field(inds) = 0;
        else
            warning('field %s had %d dimensions. This case is not currently handled, so field was ignored.',...
                fields{ii}, ndims(old_field));
        end
    end
    
    dqmc_results_grid.(fields{ii}) = new_field;
end

dqmc_results_grid.identifier = sprintf('dqmc_nsites=%dx%d_T=%0.1f-%0.1f_U=%0.1f-%0.1f_muavg=%0.1f-%0.1f_mudelta=%0.1f-%0.1f_npass=%d-%d',...
    dqmc_results_grid.nx_sites(1), dqmc_results_grid.ny_sites(1),...
    min(dqmc_results_grid.Ts_Grid(:)), max(dqmc_results_grid.Ts_Grid(:)),...
    min(dqmc_results_grid.Us_Grid(:)), max(dqmc_results_grid.Us_Grid(:)),...
    0.5 * min(dqmc_results_grid.Mus_Up_Grid(:) + dqmc_results_grid.Mus_Dn_Grid(:)),...
    0.5 * max(dqmc_results_grid.Mus_Up_Grid(:) + dqmc_results_grid.Mus_Dn_Grid(:)),...
    min(dqmc_results_grid.Mus_Up_Grid(:) - dqmc_results_grid.Mus_Dn_Grid(:)),...
    max(dqmc_results_grid.Mus_Up_Grid(:) - dqmc_results_grid.Mus_Dn_Grid(:)),...
    min(dqmc_results_grid.n_measurement_sweeps(:)), max(dqmc_results_grid.n_measurement_sweeps(:)));

end