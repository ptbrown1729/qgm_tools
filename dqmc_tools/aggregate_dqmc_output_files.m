function [dqmc_results] = aggregate_dqmc_output_files(dir_paths, file_pattern,...
                                                    T_lims, U_lims, Mu_up_lims, Mu_dn_lims, Delta_mu_lims)
% read dqmc output files, aggregate them and store in a new sturcture. Save
% this structure for further processing.
%
% dir_paths is a cell array of paths to search for output files
%
% file_pattern is a wildcard expression describing what file names are
% output files. If you do not provide an argument, '*.out' is assumed.
% 
% dir_paths: cell array of paths to read dqmc output files from
%
% file_pattern: file pattern to match output files. Typically "*.out"
%
% T_lims: length two array giving smallest and largest T values to include
%
% U_limts:
%
% Mu_up_lims:
%
% Mu_dn_lims:
%
% Delta_mu_lims:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% argument checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(dir_paths)
    dir_paths = {dir_paths};
end

if ~exist('file_pattern', 'var') || isempty(file_pattern)
    file_pattern = '*.out';
end

if ~exist('T_lims', 'var') || isempty(T_lims)
    T_lims = [-inf, inf];
end

if ~exist('U_lims', 'var') || isempty(U_lims)
    U_lims = [-inf, inf];
end

if ~exist('Mu_up_lims', 'var') || isempty(Mu_up_lims)
    Mu_up_lims = [-inf, inf];
end

if ~exist('Mu_dn_lims', 'var') || isempty(Mu_dn_lims)
    Mu_dn_lims = [-inf, inf];
end

if ~exist('Delta_mu_lims', 'var') || isempty(Delta_mu_lims)
    Delta_mu_lims = [-inf, inf];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% list files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Files = [];
for ii = 1 : length(dir_paths)
    Files = vertcat(Files, dir(fullfile(dir_paths{ii}, file_pattern)));
end

num_files = length(Files);

if num_files == 0
    dqmc_results = [];
    warning('no files matched the file pattern %s in the supplied directory.', file_pattern);
    return
end

% create stack of structures
dqmc_stack = [];
failed_indices = [];
tic
for ii = 1:num_files
    data = convert_qmc_file(fullfile(Files(ii).folder, Files(ii).name));
    if isempty(data)
        fprintf('File %s was not succesfully converted\n', Files(ii).name);
        failed_indices = horzcat(failed_indices, ii);
    end
    dqmc_stack = vertcat(dqmc_stack, data);
    tt = toc;
    fprintf('processing %d/%d . About %.1f min left.\n', ii, num_files, tt / 60 * (num_files / ii - 1));
end

% remove data that is outside the range we are interested in.
check_limits = 1;
if check_limits
    T_min = T_lims(1);
    T_max = T_lims(2);
    U_min = U_lims(1);
    U_max = U_lims(2);
    mu_ups_min = Mu_up_lims(1);
    mu_ups_max = Mu_up_lims(2);
    mu_dns_min = Mu_dn_lims(1);
    mu_dns_max = Mu_dn_lims(2);
    mu_deltas_min = Delta_mu_lims(1);
    mu_deltas_max = Delta_mu_lims(2);

    U_init = vertcat(dqmc_stack.U);
    T_init = 1./ (vertcat(dqmc_stack.L) .* vertcat(dqmc_stack.dtau));
    mu_ups_init = vertcat(dqmc_stack.mu_up);
    mu_dns_init = vertcat(dqmc_stack.mu_dn);
    mu_deltas_init = 0.5 * (mu_ups_init - mu_dns_init);
    include_logical_mat = T_init >= T_min & T_init <= T_max ...
                    & U_init >= U_min & U_init <= U_max ...
                    & mu_ups_init >= mu_ups_min & mu_ups_init <= mu_ups_max ...
                    & mu_dns_init >= mu_dns_min & mu_dns_init <= mu_dns_max ...
                    & mu_deltas_init >= mu_deltas_min & mu_deltas_init <= mu_deltas_max;

    dqmc_stack = dqmc_stack(include_logical_mat);
    
    if isempty(dqmc_stack)
        warning('Parameter limit arguments excluded all loaded datasets.');
        return;
    end
end

% add substructures using symmetry between spin up and down in certain
% cases
% dqmc_stack = swap_up_dn(dqmc_stack);

% extract substructures
latt = vertcat(dqmc_stack.latt);
equal_time_measurements = vertcat(dqmc_stack.equal_time_measurements);
equal_time_sign = vertcat(dqmc_stack.sign_of_equal_time_measurements);
double_double_correlation_function = vertcat(dqmc_stack.double_double_correlation_function);
single_single_correlation_function = vertcat(dqmc_stack.single_single_correlation_function);
density_density_correlation_fn_up_dn = vertcat(dqmc_stack.density_density_correlation_fn_up_dn);
density_density_correlation_fn = vertcat(dqmc_stack.density_density_correlation_fn_up_up);
zz_spin_correlation_function = vertcat(dqmc_stack.zz_spin_correlation_function);
pairing_correlation_function = vertcat(dqmc_stack.pairing_correlation_function);

% correlation matrices from sub-structures
calculate_corr_matrices = 1;
if calculate_corr_matrices
    double_double_correlation_mat = convert_corr_stack_to_mat(double_double_correlation_function);
    single_single_correlation_mat = convert_corr_stack_to_mat(single_single_correlation_function);
    density_density_up_dn_correlation_mat = convert_corr_stack_to_mat(density_density_correlation_fn_up_dn);
    
    % this is defined as 0.5 * ( <nup nup> + <ndn ndn> )
    density_density_correlation_mat = convert_corr_stack_to_mat(density_density_correlation_fn);
    zz_spin_correlation_mat = convert_corr_stack_to_mat(zz_spin_correlation_function);
    pairing_correlation_mat = convert_corr_stack_to_mat(pairing_correlation_function);
end

% extract data and store it in a new structure
dqmc_results = struct();
dqmc_results.dir_path = dir_paths;
dqmc_results.file_pattern = file_pattern;
dqmc_results.num_files = num_files;
% input parameters
dqmc_results.U = vertcat(dqmc_stack.U);
dqmc_results.T = 1./ (vertcat(dqmc_stack.L) .* vertcat(dqmc_stack.dtau));
dqmc_results.Mu_Up = vertcat(dqmc_stack.mu_up);
dqmc_results.Mu_Dn = vertcat(dqmc_stack.mu_dn);
dqmc_results.t_up = vertcat(dqmc_stack.t_up);
dqmc_results.t_dn = vertcat(dqmc_stack.t_dn);
% system info, size, sign, etc.
dqmc_results.seed = vertcat(dqmc_stack.seed);
dqmc_results.nx_sites = vertcat(latt.Nx);
dqmc_results.ny_sites = vertcat(latt.Ny);
dqmc_results.n_warmup_sweeps = vertcat(dqmc_stack.number_of_warmup_sweep);
dqmc_results.n_measurement_sweeps = vertcat(dqmc_stack.number_of_measurement_sweep);
dqmc_results.avg_sign = vertcat(equal_time_sign.avg_sign);
% n x 2 pairs of values and uncertainties
dqmc_results.Energy = vertcat(equal_time_measurements.potential_energy);
dqmc_results.SpinUpOcc = vertcat(equal_time_measurements.up_spin_occupancy);
dqmc_results.SpinDownOcc = vertcat(equal_time_measurements.down_spin_occupancy);
dqmc_results.DoubDensity = vertcat(equal_time_measurements.double_occupancy);
dqmc_results.Density = vertcat(equal_time_measurements.density);
% derived quantities with uncertainties
dqmc_results.SingDensity = zeros(size(dqmc_results.Density));
dqmc_results.SingDensity(:, 1) = dqmc_results.Density(:, 1) - 2 * dqmc_results.DoubDensity(:, 1);
dqmc_results.SingDensity(:, 2) = sqrt(dqmc_results.Density(:, 2).^2 + (2*dqmc_results.DoubDensity(:, 2)).^2);
dqmc_results.SingFrac = zeros(size(dqmc_results.Density));
dqmc_results.SingFrac(:, 1) = dqmc_results.SingDensity(:, 1) ./ dqmc_results.Density(:, 1);
dqmc_results.SingFrac(:, 2) = dqmc_results.SingFrac(:, 1) .* sqrt( (dqmc_results.SingDensity(:, 2) ./ dqmc_results.SingDensity(:, 1)).^ 2 + (dqmc_results.Density(:, 2) ./ dqmc_results.Density(:, 1)).^ 2);

% correlators
% appropriately normalize correlators
factor = 1;

if calculate_corr_matrices
% correlation matrices
    % doublon correlator
    dqmc_results.dd_corr_mat = factor * corr_mat_subtract_avg(double_double_correlation_mat, dqmc_results.DoubDensity);
    % singles-density correlator
    dqmc_results.ss_corr_mat = factor * corr_mat_subtract_avg(single_single_correlation_mat, dqmc_results.SingDensity);
    % nup nup correlators.
    % TODO: note that this definition is not correct in the case of a
    % spin-imbalanced system.
    dqmc_results.nupnup_corr_mat = factor * corr_mat_subtract_avg(density_density_correlation_mat, dqmc_results.SpinUpOcc);
    % charge susceptibility, spin ups, using fluctuation dissipation
    % theorem
    dqmc_results.nup_susc = (4 * squeeze( sum(sum( density_density_correlation_mat(:, 2 : end, :, 1), 2), 3) ) ...
                                 + squeeze( density_density_correlation_mat(:, 1, 1, 1)) )...
                                 ./ dqmc_results.T;
    % <n_up n_dn>
    dqmc_results.nupndn_corr_mat = factor * corr_mat_subtract_avg(density_density_up_dn_correlation_mat, dqmc_results.SpinUpOcc, dqmc_results.SpinDownOcc);
    % <s_z s_z> where s_z = (n_up - n_dn)
    % TODO: infer there is no factor of 0.5 because <Sz(0,0) Sz(0,1)> -
    % <Sz>^2 = -0.3 for U/t = 8, T/t = 0.25
    dqmc_results.SzSz_corr_mat = factor * corr_mat_subtract_avg(zz_spin_correlation_mat, dqmc_results.SpinUpOcc - dqmc_results.SpinDownOcc);

    % get static spin structure factor
    [qx, qy, sf, ~] = get_structure_fact(permute(dqmc_results.SzSz_corr_mat, [2, 3, 1, 4]), [], 'quarter');
    dqmc_results.SzSz_sfact = permute(sf, [3, 1, 2, 4]);
    dqmc_results.qx = qx;
    dqmc_results.qy = qy;
    
    % we can calculate <n_dn n_dn> from the other correlators
    % TODO: this definition is not correct
    % TODO: calculate uncertainty appropriately
    dqmc_results.ndnndn_corr_mat = dqmc_results.SzSz_corr_mat - dqmc_results.nupnup_corr_mat ...
            + dqmc_results.nupndn_corr_mat + permute(dqmc_results.nupndn_corr_mat, [1, 3, 2, 4]);
    
    % magnetic susceptibility from fluctuation dissipation theorem
    dqmc_results.magnetic_susc = (4 * squeeze(sum(sum(zz_spin_correlation_mat(:, 2:end, :, 1), 2), 3)) ...
                                 + squeeze(zz_spin_correlation_mat(:, 1, 1, 1)))...
                                 ./ dqmc_results.T;
    % pairing correlations
    dqmc_results.pairing_corr_mat = pairing_correlation_mat;
end

% NN and NNN correlators
% these are redundant since we have the correlation matrices
% retain for backwards compatibility
dqmc_results.NNdd = squeeze(dqmc_results.dd_corr_mat(:, 1, 2, :));
dqmc_results.NNNdd = squeeze(dqmc_results.dd_corr_mat(:, 2, 2, :));
dqmc_results.NNss = squeeze(dqmc_results.ss_corr_mat(:, 1, 2, :));
dqmc_results.NNupup = squeeze(dqmc_results.nupnup_corr_mat(:, 1, 2, :));
dqmc_results.NNdndn = squeeze(dqmc_results.ndnndn_corr_mat(:, 1, 2, :));
dqmc_results.NNupdn = squeeze(dqmc_results.nupndn_corr_mat(:, 1, 2, :));
dqmc_results.NNSzSz = squeeze(dqmc_results.SzSz_corr_mat(:, 1, 2, :));
% think this one had a typo, which is why it is the diagonal site instead
% of the nearest neighbor. Leaving it for now.
dqmc_results.NNPairCorr = squeeze(dqmc_results.pairing_corr_mat(:, 2, 2, :));

% create identifier
dqmc_results.identifier = sprintf('dqmc_nsites=%dx%d_T=%0.1f-%0.1f_U=%0.1f-%0.1f_muup=%0.1f-%0.1f_mudn=%0.1f-%0.1f_npass=%d',...
    dqmc_results.nx_sites(1), dqmc_results.ny_sites(1),...
    min(dqmc_results.T(:)), max(dqmc_results.T(:)),...
    min(dqmc_results.U(:)), max(dqmc_results.U(:)),...
    min(dqmc_results.Mu_Up(:, 1)), max(dqmc_results.Mu_Up(:, 1)),...
    min(dqmc_results.Mu_Dn(:, 1)), max(dqmc_results.Mu_Dn(:, 1)),...
    dqmc_results.n_measurement_sweeps(1));

end

function corr_mat = convert_corr_stack_to_mat(corr_struct_stack)
    % convert stack of correlation data to matrix format
    
    inds = [];
    corr_regex = 'point__dx__(\d+)__dy__(\d+)';
    struct_fields = fields(corr_struct_stack);
    for ii = 1:length(struct_fields)
        m = regexp(struct_fields{ii}, corr_regex, 'tokens');
        inds = vertcat(inds, [str2num(m{1}{1}), str2num(m{1}{2})]);
    end

    mat_size = max(inds(:)) + 1;
    corr_mat = zeros(length(corr_struct_stack), mat_size, mat_size, 2);
    for ii = 1:length(inds)
        xx = inds(ii, 1) + 1; % matlab indexing starts at 1 instead of 0
        yy = inds(ii, 2) + 1;
        corr_mat(:, yy, xx, :) = vertcat(corr_struct_stack.(struct_fields{ii}));
        corr_mat(:, xx, yy, :) = corr_mat(:, yy, xx, :);
    end
end

function corr_mat_minus_avg = corr_mat_subtract_avg(corr_mat, avg1, avg2)
    % corr_mat is an N x nsites x nsites x 2 correlation matrix
    % avg is an N x 2 vector

    % expand average values to be the same size as corr_mat
    avg1_expanded = permute(...
                   repmat(avg1, [1, 1, size(corr_mat, 2), size(corr_mat, 3)]),...
                          [1, 3, 4, 2]);
    
    avg1_is_avg2 = 0;
    if ~exist('avg2', 'var') || isempty(avg2)
        avg1_is_avg2 = 1;
    end
        
    if avg1_is_avg2
        avg2_expanded = avg1_expanded;
    else
        avg2_expanded = permute(...
               repmat(avg2, [1, 1, size(corr_mat, 2), size(corr_mat, 3)]),...
                      [1, 3, 4, 2]);
    end
    % value calculation
    corr_mat_minus_avg(:, :, :, 1) = corr_mat(:, :, :, 1) - avg1_expanded(:, :, :, 1) .* avg2_expanded(:, :, :, 1);
    % uncertainty calculation
    if avg1_is_avg2
        corr_mat_minus_avg(:, :, :, 2) = ...
            sqrt( corr_mat(:, :, :, 2) .^ 2 + ...
                  ( 2 * avg1_expanded(:, :, :, 1) .* avg1_expanded(:, :, :, 2) ) .^ 2 ); 
    else
        corr_mat_minus_avg(:, :, :, 2) = ...
            sqrt( corr_mat(:, :, :, 2) .^ 2 + ...
                  ( avg1_expanded(:, :, :, 1) .* avg2_expanded(:, :, :, 2) ) .^ 2 + ...
                  ( avg2_expanded(:, :, :, 1) .* avg1_expanded(:, :, :, 2) ) .^ 2); 
    end
end