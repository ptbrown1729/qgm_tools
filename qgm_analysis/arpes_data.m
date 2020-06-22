classdef arpes_data < handle
    %arpes_data class for storing and manipulating results of arpes
    %experiments
    
    properties
        % data sets
        identifier
        date
        
        % analysis mode
        analysis_mode = 'resample' % 'resample' or 'round'
        nsample_pixels = 17 % only used if 'resample'
        img_bin_size = 4 % only used if 'round'
        n_energy_bins = 20
        use_background = 1
        use_expt_broadening = 1
        broadening_khz = 1
        symmetry = 'd4' % 'd4', 'd2', or 'none'
        
        % physical data
        latt_depth
        rf_time
        rf_power
        heating_holdt
        heating_modt
        
        % arpes data
        init_hfstate_lower_energy = 1
        frqs_hz_sorted
        frqs_offsets_tunits
        data_folder_stack
        bg_folder
        
        % non-interacting gas rf spectroscopy
        rfspec_df
        rfspec_frqs
        rfspec_xfer
        rfspec_xfer_unc
        rfspec_lor_fitp
        rfspec_lor_fitp_unc 
        rfspec_dbl_lor_fitp
        rfspec_dbl_lor_fitp_unc
        frq_nonint_transition_hz
        
        % interacting gas rf spectroscopy
        rfspec_int_df
        rfspec_int_frqs
        
        % values derived from interacting rf spectroscopy
        rfxfer_density_resolved % atom number
        rfxfer_density_resolved_density % density
        rfxfer_density_resolved_density_unc
        rfxfer_density_resolved_lor_fitp
        rfxfer_density_resolved_frqs
        
        % temperature fitting
        static_dqmc_path
        ones
        threes
        singles
        doubles
        temp_fit_struct
        mu_struct
        mu_center        
        
        % constants
        m = 9.9883e-27;
        latt_const = 750e-9;
        hbar = 1.0546e-34;
        
        % parameters
        omega = 2 * pi * 362; % bottom beam trapping frequency
        t_over_h = inf; % hopping rate in Hz
        latt_dispersion = @(kx, ky) 2*(2 - cos(kx) - cos(ky));
        
        k_bz_edge
        x_bz_edge
        y_bz_edge
        cx
        cy
        
        % 2d quantities          
        dx
        dy
        dkx
        dky
        energy_bin_means
        energy_bin_edges
        
        specfn_unsymm % unsymmetrized, unequal energies
        specfn_uneq_es % symmetrized, unequal energies
        specfn_uneq_scaled % symmetrized, unequal energies, scaled
        specfn % symmetrized, equal energies
        specfn_scaled %symmetrized, equal energies, scaled to match theory
        
        % momentum summed quantities
        specfn_sumk
        specfn_sumk_fitp
        specfn_sumk_fitp_err
        
        rfxfer_sumk
        rfxfer_sumk_fitp
        rfxfer_sumk_fitp_err
        
        % lorentzian fits to dispersion
        use_double_lor_dispersion = 0
        
        lor_fit_params_gxmg
        lor_fit_params_uncs_gxmg
        lor_fit_params_gymg
        lor_fit_params_uncs_gymg
        
        qp_dispersion_af_gxmg
        qp_dispersion_af_gxmg_unc
        qp_dispersion_af_gymg
        qp_dispersion_af_gymg_unc
        
        qp_dispersion_a_gxmg
        qp_dispersion_a_gxmg_unc
        qp_dispersion_a_gymg
        qp_dispersion_a_gymg_unc
       
        % bcs dispersion fit
        fit_params_bcs_gxmg
        stderr_bcs_gxmg
        linear_index_bcs_gxmg
        bcs_fn_sample_gxmg
        
        fit_params_bcs_gymg
        stderr_bcs_gymg
        linear_index_bcs_gymg
        bcs_fn_sample_gymg
        
        % dqmc theory parameters
        dqmc_dynamic_path
        dqmc
        theory_temp_index  
        dqmc_density
        density_area_fraction
        theory_scale_fitparams % [scale factor, vertical offset, energy offset]
        theory_scale_stderrs

        % spectralfn structure
        dqmc_aavg % trap avg theory scaled to match normal identities (e.g. integral normalized to 1)
        dqmc_aavg_nobroad
        dqmc_aavg_scaled % comparable to 2D pics
        dqmc_a_dispersion_gxmg_fitp
        dqmc_a_dispersion_gxmg_stderr
        dqmc_af_dispersion_gxmg_fitp
        dqmc_af_dispersion_gxmg_stderr
        max_num_peaks = 4;
        
    end
    
    methods
        
        % load data
        function obj = arpes_data(date_cell, folders, frqs, bin_edges,...
                                  exclude_pics, az_avg_mode, az_avg_grid)
            
            % check arguments
            if ~exist('az_avg_mode', 'var')
                az_avg_mode = '';
            end
            
            if ~exist('az_avg_grid', 'var')
                az_avg_grid = [];
            end
           
            
            % initialize if arguments are supplied. Otherwise, instantiate
            % an empty object. This allows you to change the settings
            % before doing data analysis.
            if ~exist('date_cell', 'var')
                return;
            end
            
            % load folders          
            obj.load_arpes_folders(date_cell, folders, frqs, bin_edges,...
                             exclude_pics, az_avg_mode, az_avg_grid);
            % analyze arpes
            obj.analyze_arpes();
        end
         
        function load_arpes_folders(obj, date_cell, folders, frqs_hz, bin_edges, exclude_pics, roi_start_coords, bg_folder, bg_exclude_pics)
            % load folder data and store results in instance
            % Arguments
            % date_cell: is a cell array of length 3 of the form {yyyy, mm, dd}
            %
            % folders: is an array of indices to the folders to be analyzed.
            % It is assumed all folders are associated with the date in
            % date_cell. 
            %
            % frqs: Frequencies in Hertz for each folder in folders
            %
            % bin_edges: edge distances for azimuthal average bins used for
            % all folders. These should not effect quantities calculated
            % from e.g. the 2D momentum distribution
            %
            % exclude_pics: a cell array of vectors. Each vector contains
            % indices to pictures that should be excluded from the
            % respective folder. e.g. exclude_pics = {[1, 2, 3], [2], []}
            % will exclude picture indices 1, 2, and 3 from folders(1) and
            % picture index 2 from folders(2).
           
            obj.date = date_cell;
            [obj.frqs_hz_sorted, I] = sort(frqs_hz); 
            
            % begin processing
            obj.identifier = sprintf('%4d-%02d-%02d_%03d-%03d_arpes',...
                obj.date{1}, obj.date{2}, obj.date{3}, min(folders), max(folders));
            
            if isempty(folders)
                return;
            end
            
            if ~exist('roi_start_coords', 'var') || isempty(roi_start_coords)
                roi_start_coords = [124, 107];
            end
            
            obj.data_folder_stack = [];
            for ii = 1:length(folders)
                df = DataFolder();
                df.CroppedPicStartCoords = roi_start_coords;
                df.CenterStyle = 'Fixed';
                df.initialize({obj.date{1}, obj.date{2}, obj.date{3}, folders(ii), 1, 1},...
                               exclude_pics{ii}, bin_edges);
                obj.data_folder_stack = horzcat(obj.data_folder_stack, df);
            end
            % sort data folder stack
            obj.data_folder_stack = obj.data_folder_stack(I);
            
            if exist('bg_folder', 'var') && ~isempty(bg_folder)
                if ~exist('bg_exclude_pics', 'var')
                    bg_exclude_pics = [];
                end
                
                df = DataFolder();
                df.CroppedPicStartCoords = roi_start_coords;
                df.CenterStyle = 'Fixed';
                df.initialize({obj.date{1}, obj.date{2}, obj.date{3}, bg_folder, 1, 1},...
                               bg_exclude_pics, bin_edges);
                obj.bg_folder = df;
            else
                if obj.use_background
                    error('no background folder supplied');
                end
            end
        end
        
        function load_rf_nonint(obj, folder, frqs, date, bin_edges, exclude_pics, roi_start_coords)
            % load non-interacting rf spectroscopy data and store in
            % instance
            
            if ~exist('date', 'var')
                if ~isempty(obj.date)
                    date = obj.date;
                else
                    error('must provide date to load_rf_nonint');
                end  
            end
            
            if ~exist('exclude_pics', 'var')
                exclude_pics = [];
            end
            
            if ~exist('roi_start_coords', 'var') || isempty(roi_start_coords)
                roi_start_coords = [123, 111];
            end
            
            obj.rfspec_df = DataFolder();
            obj.rfspec_df.CenterStyle = 'Fixed';
            obj.rfspec_df.CroppedPicStartCoords = roi_start_coords;
            obj.rfspec_df.initialize(...
                {date{1}, date{2}, date{3}, folder, 1, 1},...
                exclude_pics, bin_edges);
            obj.rfspec_df.IndependentVariable = frqs;
        end
        
        function load_rf_int(obj, folder, frqs, date, bin_edges, exclude_pics)
            % load interacting rf folder and store in instance
            
            if ~exist('date', 'var')
                if ~isempty(obj.date)
                    date = obj.date;
                else
                    error('must provide date to load_rf_nonint');
                end  
            end
            
            if ~exist('exclude_pics', 'var')
                exclude_pics = [];
            end
            
            obj.rfspec_int_df = DataFolder();
            obj.rfspec_int_df.CenterStyle = 'Fixed';
            obj.rfspec_int_df.CroppedPicStartCoords = [123, 111];
            obj.rfspec_int_df.initialize(...
                {date{1}, date{2}, date{3}, folder, 1, 1},...
                exclude_pics, bin_edges);
            obj.rfspec_int_df.IndependentVariable = frqs;
        end
        
        function load_temp_folders(obj, folders, date, bin_edges, exclude_pics)
                % load temperature data folders and store in instance.
                
                if ~exist('date', 'var')
                    if ~isempty(obj.date)
                        date = obj.date;
                    else
                        error('must provide date to load_rf_nonint');
                    end  
                end
            
                if ~exist('exclude_pics', 'var')
                    exclude_pics = {[], [], [], []};
                end
            
                obj.ones = DataFolder({date{1}, date{2}, date{3}, folders(1), 1, 1},...
                                       exclude_pics{1}, bin_edges);
                
                obj.threes = DataFolder();
                obj.threes.CenterStyle = 'Fixed';
                obj.threes.CroppedPicStartCoords = obj.ones.CroppedPicStartCoords;
                obj.threes.initialize({date{1}, date{2}, date{3}, folders(2), 1, 1},...
                                       exclude_pics{2}, bin_edges);
                
                obj.singles = DataFolder();
                obj.singles.CenterStyle = 'Fixed';
                obj.singles.CroppedPicStartCoords = obj.ones.CroppedPicStartCoords;
                obj.singles.initialize({date{1}, date{2}, date{3}, folders(3), 1, 1},...
                                        exclude_pics{3}, bin_edges);
                
                obj.doubles = DataFolder();
                obj.doubles.CenterStyle = 'Fixed';
                obj.doubles.CroppedPicStartCoords = obj.ones.CroppedPicStartCoords;
                obj.doubles.initialize({date{1}, date{2}, date{3}, folders(4), 1, 1},...
                                        exclude_pics{4}, bin_edges);
        end
        
        function load_dqmc_arpes_dat(obj, dqmc)
            % load DQMC dynamical data 
            %
            % dqmc: may be either a string giving the path to a matlab file
            % to load, or a matlab structure (i.e. the already loaded
            % file).
            %
            % load DQMC data. Note that dqmc is the one field of this class
            % which should not be saved with the class because it will make
            % the file too big. I made it a field so I can load the dqmc
            % data initially and then perform multiple analyses without
            % reloading (e.g. if I didn't store it I would have to reload
            % it). One solution might be to handle loading outside of the
            % class. But what I'm doing here seems simplest at the moment.
            
            % load data
            if ~isempty(dqmc) && ischar(dqmc)
                % if dqmc is path to dqmc saved data
                obj.dqmc_dynamic_path = dqmc;
                s = load( obj.dqmc_dynamic_path );  
                obj.dqmc =  spectralfn();
                obj.dqmc.load_struct(s);
            else
                % if dqmc is already a struct
                obj.dqmc = dqmc;
            end

        end
        
        function get_trap_avg_thry(obj)           
            % Get theoy quantities which depend on parameters of our
            % dataset. Call this function after load_dqmc_arpes_dat and
            % analyze_arpes have been run.
            
            % get nearest temperature
            if isempty(obj.dqmc)
                return;
            end
            
            % get temperature
            temp = obj.temp_fit_struct.fit_params(2);
            
%             [~, obj.theory_temp_index] = min( abs( temp - 1 ./ obj.dqmc.betas) ); 
            [~, obj.theory_temp_index] = min( abs( temp - obj.dqmc.ts) );

            % bin trap density to get a more accurate measure
            eff = obj.temp_fit_struct.fit_params(4);
            densities = 2 * obj.ones.Occs_ImgAvg / eff;
            densities = binImg(densities, 4, 4, 'Same');
            
            % get trap average quantities
%             a_fn = squeeze( obj.dqmc.a_interp(:, obj.theory_temp_index, :, :, :) );
%             af_fn = squeeze( obj.dqmc.af_interp(:, obj.theory_temp_index, :, :, :) );
%     
%             sfn = spectralfn( obj.dqmc.kxs, obj.dqmc.kys, obj.dqmc.omega_interp, obj.dqmc.ns,...
%                               obj.dqmc.mus, 1 ./ obj.dqmc.betas(obj.theory_temp_index), obj.dqmc.U,...
%                               permute(a_fn, [2, 3, 4, 1]),...
%                               permute(af_fn, [2, 3, 4, 1]), '', '', 'quad');    
    
            sfn = spectralfn( obj.dqmc.kxs, obj.dqmc.kys, obj.dqmc.es, obj.dqmc.ns,...
                              obj.dqmc.mus, obj.dqmc.ts(obj.theory_temp_index), obj.dqmc.us,...
                              squeeze( obj.dqmc.a(:, :, :, :, obj.theory_temp_index) ),...
                              squeeze( obj.dqmc.af(:, :, :, :, obj.theory_temp_index) ),...
                              '', '', 'quad');

            [a, obj.dqmc_density, areas] = get_trapavg_spectralfn(densities, '', permute(sfn.a, [4, 1, 2, 3]), sfn.ns, 'mean');
            [af, ~, ~] = get_trapavg_spectralfn(densities, '', permute(sfn.af, [4, 1, 2, 3]), sfn.ns, 'mean');
            obj.dqmc_aavg = spectralfn( sfn.kxs, sfn.kys, sfn.es, '', '', sfn.ts, sfn.us, a, af, '', '', 'quad');
            obj.dqmc_aavg_nobroad = obj.dqmc_aavg;
            
            obj.density_area_fraction = areas / nansum(areas);

            % include experimental broadening
            if obj.use_expt_broadening && obj.broadening_khz > 0
                sigma = obj.broadening_khz * 1e3 / obj.t_over_h;
%                 dw = obj.dqmc.omega_interp(2) - obj.dqmc.omega_interp(1);
                dw = obj.dqmc.es(2) - obj.dqmc.es(1);

                a = convolve_gauss( permute(obj.dqmc_aavg.a, [3, 1, 2]), dw, sigma);
                a = permute(a, [2, 3, 1]);
                
                af = convolve_gauss( permute(obj.dqmc_aavg.af, [3, 1, 2]), dw, sigma);
                af = permute(af, [2, 3, 1]);
                
                obj.dqmc_aavg = spectralfn( obj.dqmc_aavg.kxs, obj.dqmc_aavg.kys,...
                                            obj.dqmc_aavg.es, '', '',...
                                            obj.dqmc_aavg.ts, obj.dqmc_aavg.us,...
                                            a, af, '', '', 'quad');
            end
            
            obj.dqmc_aavg.quad2full();

           % fit dispersions along high symmetry cuts
%            dist_fn = @(x1, y1, x2, y2) min( mod( x1 - x2, 3), mod(x2 - x1, 3) );
%            order = dist_fn(obj.dqmc_aavg.linear_index, 0, 0, 0);
%            % guesses
%            [~, iii] = max( obj.dqmc_aavg.a_gxmg(1, :)); 
%            init_params = [obj.dqmc_aavg.es(iii), 1, max(obj.dqmc_aavg.a_gxmg(1, :)), 0];
%            fixed_params = [0, 0, 0, 1];
%            lbs = [min(obj.dqmc_aavg.es), 0.5, 0, -inf];
%            ubs = [max(obj.dqmc_aavg.es), inf, inf, inf];
%            % fitting
%            [obj.dqmc_af_dispersion_gxmg_fitp, obj.dqmc_af_dispersion_gxmg_stderr, ~] = ...
%            fit_array_2lor(obj.dqmc_aavg.af_gxmg, '', obj.dqmc_aavg.es, obj.dqmc_aavg.linear_index, '',...
%                           init_params, fixed_params, lbs, ubs, 1, 1, dist_fn, order);  
%            
% 
%            init_params = [init_params, init_params + [2, 0, 0, 0]];
%            fixed_params = [fixed_params, fixed_params];
%            lbs = [lbs, lbs];
%            ubs = [ubs, ubs];
%            [obj.dqmc_a_dispersion_gxmg_fitp, obj.dqmc_a_dispersion_gxmg_stderr, ~] = ...
%            fit_array_2lor(obj.dqmc_aavg.a_gxmg, '', obj.dqmc_aavg.es, obj.dqmc_aavg.linear_index, '',...
%                           init_params, fixed_params, lbs, ubs, 2, 1, dist_fn, order);
       
           % maybe better to do with find peaks in this case
           af_params = zeros( length(obj.dqmc_aavg.linear_index), 4 * obj.max_num_peaks);
           a_params = zeros( length(obj.dqmc_aavg.linear_index), 4 * obj.max_num_peaks);

           min_peak_height = 0.001;
           
           for ii = 1 : length(obj.dqmc_aavg.linear_index)
               [hs, ps, ws, proms] = findpeaks(obj.dqmc_aavg.a_gxmg(ii, :), obj.dqmc_aavg.es);
               params = [ps(1), ws(1), hs(1), proms(1)];
               for jj = 1 : obj.max_num_peaks-1
                   if length(ps) > jj && hs(jj + 1) > min_peak_height
                       params = [params, [ps(jj + 1), ws(jj + 1), hs(jj + 1), proms(jj + 1)] ];
                   else
                       params = [params, [nan, 0, 0, 0]];
                   end
               end
               a_params(ii, :) = params;
               
               
              [hs, ps, ws, proms] = findpeaks(obj.dqmc_aavg.af_gxmg(ii, :), obj.dqmc_aavg.es);
               params = [ps(1), ws(1), hs(1), proms(1)];
               for jj = 1 : obj.max_num_peaks-1
                   if length(ps) > jj && hs(jj + 1) > min_peak_height
                       params = [params, [ps(jj + 1), ws(jj + 1), hs(jj + 1), proms(jj + 1)] ];
                   else
                       params = [params, [nan, 0, 0, 0]];
                   end
               end
               af_params(ii, :) = params;  
           end
          
           obj.dqmc_af_dispersion_gxmg_fitp = af_params;
           obj.dqmc_af_dispersion_gxmg_stderr = zeros( size(af_params) );
           obj.dqmc_a_dispersion_gxmg_fitp = a_params;
           obj.dqmc_a_dispersion_gxmg_stderr = zeros( size(a_params) );
           
                      
            
%            obj.dqmc_a_disp = a_fit_params_gxmg();
                      
%            array_af_interp = zeros( size( obj.dqmc_aavg.af_gxmg) );
%            array_a_interp = zeros( size( obj.dqmc_aavg.af_gxmg) );
%            for ii = 1 : size(array_af_interp, 1)
%                    array_af_interp(ii, :) = lorentzian1D( af_fit_params_gxmg(ii, :), obj.dqmc_aavg.es);
%                    array_a_interp(ii, :) = lorentzian1D( a_fit_params_gxmg(ii, :), obj.dqmc_aavg.es);
%            end    
%            
%             figh = plot_3d_array(obj.dqmc_aavg.af_gxmg, '', obj.dqmc_aavg.linear_index, '', obj.dqmc_aavg.es, '', '', '', 'b.');
%             plot_3d_array(array_af_interp, '', obj.dqmc_aavg.linear_index, '', obj.dqmc_aavg.es, '', figh, '', 'r');
%                 
%             figh2 = plot_3d_array(obj.dqmc_aavg.a_gxmg, '', obj.dqmc_aavg.linear_index, '', obj.dqmc_aavg.es, '', '', '', 'b.');
%             plot_3d_array(array_a_interp, '', obj.dqmc_aavg.linear_index, '', obj.dqmc_aavg.es, '', figh2, '', 'r');
%             
%             figure;
%             errorbar(obj.specfn_uneq_es.linear_index, obj.qp_dispersion_af_gxmg(:, 1), obj.qp_dispersion_af_gxmg_unc(:, 1));
%             hold on;
%             errorbar(obj.dqmc_aavg.linear_index, af_fit_params_gxmg(:, 1), af_fit_se_gxmg(:, 1));
%             errorbar(obj.dqmc_aavg.linear_index, a_fit_params_gxmg(:, 1),a_fit_se_gxmg(:, 1));
%             legend({'expt', 'thry af', 'thry a'});
            
            % get scale factor for theory to expt. 
            % First try: do this by looking at maximum for (kx, ky) = (0, 0) position
%                 obj.theory_scale_trpavg = max(obj.lor_fit_params_gxmg(:, 3)) ./ max( obj.dqmc_aavg.af_gxmg(obj.theory_temp_index, 1, :));
            % Second try: fitting
            [kyky, kxkx, ee] = ndgrid( obj.dqmc_aavg.kys, obj.dqmc_aavg.kxs, obj.dqmc_aavg.es);
            interp_fn = @(kx, ky, e) interp3(kxkx, kyky, ee, obj.dqmc_aavg.af, kx, ky, e);
            
            % get closest energy alignment
            min_inds = zeros( length(obj.energy_bin_means), 1 );
            for ii = 1 : length(min_inds)
                [~, min_inds(ii)] = min( abs(obj.dqmc_aavg.es - obj.energy_bin_means(ii)) );
            end
            
%             af_to_fit = obj.dqmc_aavg.af(:, :, min_inds);
            % ignore any points that have nans
            use_indices = ~isnan( obj.specfn.af );
            
            % using dqmc kxs and kys because expect them to very close to
            % experimental values, and running into trouble with
            % experimental values bc slightly outisde of [-pi, pi]
            [kxkx_act, kyky_act, ee_act] = ndgrid( obj.dqmc_aavg.kys, obj.dqmc_aavg.kxs, obj.specfn.es);
            use_indices( ee_act > max(obj.dqmc_aavg.es) | ee_act < min(obj.dqmc_aavg.es) ) = 0;
            
            kxkx_act = kxkx_act(use_indices);
            kyky_act = kyky_act(use_indices);
            ee_act = ee_act(use_indices);
            
            % define fit function
%             fit_fn = @(P) ( af_to_fit(use_indices) * P(1) + P(2) - obj.specfn.af(use_indices) ) ./ obj.specfn.afunc(use_indices);
            fit_fn = @(P) ( interp_fn(kxkx_act, kyky_act, ee_act - P(3)) * P(1) + P(2) - obj.specfn.af(use_indices) ) ./ obj.specfn.afunc(use_indices);
            % fit settings
            init_params = [max(obj.lor_fit_params_gxmg(:, 3)) ./ max( obj.dqmc_aavg.af_gxmg(1, :)), 0, 0];
            fixed_params = [0, 1, 1];
            lbs = [0, -0.01, -2];
            ubs = [1, 0.01, 2];

            [fit_params, std_errs, chi_sqr, flag] = lsq_fixedp(fit_fn, init_params, fixed_params, lbs, ubs);
            obj.theory_scale_fitparams = fit_params;
            obj.theory_scale_stderrs = std_errs;

            % scale experiment data to match theory height
            af_scaled = ( obj.specfn.af - fit_params(2) ) / fit_params(1);
            afunc_scaled = obj.specfn.afunc / fit_params(1);    
            e_scaled = obj.specfn.es - fit_params(3);
            % divide by fermi function to get approximate a
            e_exp = repmat(permute(e_scaled(:), [3, 2, 1]), [ size(af_scaled, 1), size(af_scaled, 2), 1]);
            temp = obj.temp_fit_struct.fit_params(2);
            fermi = @(e, T) 1 ./ ( 1 + exp(e / T) ); 
            a = af_scaled ./ fermi(e_exp - obj.mu_center, temp);
            aunc = afunc_scaled ./ fermi(e_exp - obj.mu_center, temp);
            
            obj.specfn_scaled = spectralfn( obj.specfn.kxs, obj.specfn.kys, e_scaled,...
                                            obj.specfn.ns, obj.specfn.mus,...
                                            obj.specfn.ts, obj.specfn.us,...
                                            a, af_scaled, aunc, afunc_scaled, 'full');
            
            %
            af_scaled = ( obj.specfn_uneq_es.af - fit_params(2) ) / fit_params(1);
            afunc_scaled = obj.specfn_uneq_es.afunc / fit_params(1);
            e_scaled = obj.specfn_uneq_es.es - fit_params(3);
            %
            temp = obj.temp_fit_struct.fit_params(2);
            fermi = @(e, T) 1 ./ ( 1 + exp(e / T) ); 
            a = af_scaled ./ fermi(e_scaled - obj.mu_center, temp);
            aunc = afunc_scaled ./ fermi(e_scaled - obj.mu_center, temp);
                                        
            obj.specfn_uneq_scaled = spectralfn( obj.specfn_uneq_es.kxs, obj.specfn_uneq_es.kys, e_scaled,...
                                            obj.specfn_uneq_es.ns, obj.specfn_uneq_es.mus,...
                                            obj.specfn_uneq_es.ts, obj.specfn_uneq_es.us,...
                                            a, af_scaled, aunc, afunc_scaled, 'full');
            % fit a and af dispersion
           dist_fn = @(x1, y1, x2, y2) min( mod( x1 - x2, 3), mod(x2 - x1, 3) );
           order = dist_fn(obj.specfn_uneq_scaled.linear_index, 0, 0, 0);
           % initial parameters
           init_params = [obj.lor_fit_params_gxmg(1, 1), 1, max(obj.specfn_uneq_scaled.af_gxmg(1, :)), 0];
           fixed_params = [0, 0, 0, 1];
           lbs = [min(obj.specfn_uneq_scaled.es(:)), 0.5, 0, -inf];
           ubs = [max(obj.specfn_uneq_scaled.es(:)), inf, inf, inf];
                
                     
          if obj.use_double_lor_dispersion
               n_lors = 2;
               init_params = horzcat(init_params, init_params);
               fixed_params = horzcat(fixed_params, fixed_params);
               lbs = horzcat(lbs, lbs);
               ubs = horzcat(ubs, ubs);
           else
               n_lors = 1;
           end
           
           
           % af
          [fit_params_gxmg, fit_params_uncs_gxmg, ~] = ...
           fit_array_2lor(obj.specfn_uneq_scaled.af_gxmg, '',...
                          obj.specfn_uneq_scaled.es_gxmg, obj.specfn_uneq_scaled.linear_index, '',...
                          init_params, fixed_params, lbs, ubs, n_lors, 1, dist_fn, order);
          
          if ~obj.use_double_lor_dispersion
               fit_params_gxmg = horzcat( fit_params_gxmg, zeros( size(fit_params_gxmg) ) );
               fit_params_uncs_gxmg = horzcat( fit_params_uncs_gxmg, zeros( size(fit_params_uncs_gxmg) ) );
          end
                      
          obj.qp_dispersion_af_gxmg = fit_params_gxmg(:, 1:4:end);
          obj.qp_dispersion_af_gxmg_unc = fit_params_uncs_gxmg(:, 1:4:end);          
         
          % gymg
          if strcmp(obj.symmetry, 'd4')
              obj.qp_dispersion_af_gymg = obj.qp_dispersion_af_gxmg;
              obj.qp_dispersion_af_gymg_unc = obj.qp_dispersion_af_gxmg_unc;
          else
            [fit_params_gymg, fit_params_uncs_gymg, ~] = ...
            fit_array_2lor(obj.specfn_uneq_scaled.af_gymg, '',...
                          obj.specfn_uneq_scaled.es_gymg, obj.specfn_uneq_scaled.linear_index, '',...
                          init_params, fixed_params, lbs, ubs, n_lors, 1, dist_fn, order);
                      
            if ~obj.use_double_lor_dispersion
               fit_params_gymg = horzcat( fit_params_gymg, zeros( size(fit_params_gymg) ) );
               fit_params_uncs_gymg = horzcat( fit_params_uncs_gymg, zeros( size(fit_params_uncs_gymg) ) );
            end
            
            obj.qp_dispersion_af_gymg = fit_params_gymg(:, 1:4:end);
            obj.qp_dispersion_af_gymg_unc = fit_params_uncs_gymg(:, 1:4:end);  
          end
          
           % a
           unc = obj.specfn_uneq_scaled.aunc_gxmg;
           unc( unc == 0) = 0.1;
           [fit_params_gxmg, fit_params_uncs_gxmg, ~] = ...
           fit_array_2lor(obj.specfn_uneq_scaled.a_gxmg, unc,...
                          obj.specfn_uneq_scaled.es_gxmg, obj.specfn_uneq_scaled.linear_index, '',...
                          init_params, fixed_params, lbs, ubs, n_lors, 1, dist_fn, order);
                          
              
          if ~obj.use_double_lor_dispersion
               fit_params_gxmg = horzcat( fit_params_gxmg, zeros( size(fit_params_gxmg) ) );
               fit_params_uncs_gxmg = horzcat( fit_params_uncs_gxmg, zeros( size(fit_params_uncs_gxmg) ) );
          end
          
          obj.qp_dispersion_a_gxmg = fit_params_gxmg(:, 1:4:end);
          obj.qp_dispersion_a_gxmg_unc = fit_params_uncs_gxmg(:, 1:4:end); 
            
          % gymg
          if strcmp(obj.symmetry, 'd4')
              obj.qp_dispersion_a_gymg = obj.qp_dispersion_a_gxmg;
              obj.qp_dispersion_a_gymg_unc = obj.qp_dispersion_a_gxmg_unc;
          else
            [fit_params_gymg, fit_params_uncs_gymg, ~] = ...
            fit_array_2lor(obj.specfn_uneq_scaled.a_gymg, '',...
                          obj.specfn_uneq_scaled.es_gymg, obj.specfn_uneq_scaled.linear_index, '',...
                          init_params, fixed_params, lbs, ubs, n_lors, 1, dist_fn, order);
                      
            if ~obj.use_double_lor_dispersion
               fit_params_gymg = horzcat( fit_params_gymg, zeros( size(fit_params_gymg) ) );
               fit_params_uncs_gymg = horzcat( fit_params_uncs_gymg, zeros( size(fit_params_uncs_gymg) ) );
            end
            
            obj.qp_dispersion_a_gymg = fit_params_gymg(:, 1:4:end);
            obj.qp_dispersion_a_gymg_unc = fit_params_uncs_gymg(:, 1:4:end);  
          end
            
            % also compute quantity comparable to our 2D maps
            dqmc_nr = zeros(size(obj.dqmc_aavg.af, 1), size(obj.dqmc_aavg.af, 2), length(obj.frqs_offsets_tunits));
            energies = zeros( size(dqmc_nr) );
            for jj = 1 : length(obj.dqmc_aavg.kxs)
                for ii = 1 : length(obj.dqmc_aavg.kys)
                    if obj.init_hfstate_lower_energy
                        energy_vals = obj.latt_dispersion(obj.dqmc_aavg.kxs(jj), obj.dqmc_aavg.kys(ii)) - obj.frqs_offsets_tunits;
                    else
                        energy_vals = obj.latt_dispersion(obj.dqmc_aavg.kxs(jj), obj.dqmc_aavg.kys(ii)) + obj.frqs_offsets_tunits;
                    end
                    energies(ii, jj, :) = energy_vals;
                    dqmc_nr(ii, jj, :) = interp1(obj.dqmc_aavg.es,...
                        squeeze(obj.dqmc_aavg.af(ii, jj, :)), energy_vals); 
                end
            end
     
            % scale to match density
            dqmc_nr = dqmc_nr * fit_params(1) + fit_params(2);
            energies_scaled = energies + fit_params(3);
            
            
            obj.dqmc_aavg_scaled = spectralfn( obj.dqmc_aavg.kxs, obj.dqmc_aavg.kys,...
                                            energies_scaled, '', '',...
                                            obj.dqmc_aavg.ts, obj.dqmc_aavg.us,...
                                            '', dqmc_nr, '', '', 'full');
        end
        
        % analyze all
        function analyze_all(obj, non_int_trans_frq_MHz)
            % main analysis function which runs all of the subsidiary ones
            % in correct order.
            
            if ~exist('non_int_trans_frq_MHz', 'var')
                non_int_trans_frq_MHz = '';
            end
            
            obj.process_rf_nonint();
            obj.process_rf_int();
            
            % analyze_arpes() must be run after rf_nonint
            obj.analyze_arpes(non_int_trans_frq_MHz);
            
            % get_trap_avg_thry() must be run after analyze_arpes() and
            % fig_temp_dqmc()
            obj.get_trap_avg_thry();
            
        end
        
        % analysis
        function analyze_arpes(obj, non_int_trans_frq_MHz)
            % Analyze raw images to get spectral function. This is the main 
            % experimental analysis function.    
            
            % ensure frequencies are set correctly
            if exist('non_int_trans_frq_MHz', 'var') && ~isempty(non_int_trans_frq_MHz)
                % if the non-interacting transition frequency is provided
                % as an argument, override other determination of it
                obj.frq_nonint_transition_hz =  non_int_trans_frq_MHz * 1e6;
                 if ~isempty(obj.rfspec_lor_fitp)
                     warning('overriding non-interacting spectroscopy fit value with argument value');
                 end
                
            else
                if ~isempty(obj.rfspec_lor_fitp)
                    if ~isempty(obj.frq_nonint_transition_hz) && ...
                       obj.frq_nonint_transition_hz ~= 0 && ...
                       obj.frq_nonint_transition_hz ~= obj.rfspec_lor_fitp(1) * 1e6
                        warning('non-interacting spectroscopy data present. Overriding value of frq_nonint_transition_hz.');
                    end
                   obj.frq_nonint_transition_hz = obj.rfspec_lor_fitp(1) * 1e6;
                else
                    warning('No non-interacting spectroscopy data. Not modifying obj.frq_nonint_transition_hz');
                end
            end
            
            obj.frqs_offsets_tunits = (obj.frqs_hz_sorted - obj.frq_nonint_transition_hz) / obj.t_over_h;
            
            % find mean center of all folders
            if isempty(obj.data_folder_stack)
                error('obj.data_folder_stack is empty. Either load folder data from file or reanalyze');
            end
            [obj.cx, obj.cy] = obj.find_center(obj.data_folder_stack);
            
            % get Brillouin zone edges
            % x(T/4) = hbar/(m*omega) * K(0)
            % We typically do not use the coordinate x, but rather x/a
            % therefore, (x/a) = (hbar / (m*omega*a) ) * k
            % or at the Brillouin zone edge (x_bz/a) = (hbar / (m*omega*a) ) * k/pi
            xtok = (obj.m * obj.omega * obj.latt_const) / obj.hbar * obj.latt_const;
            obj.x_bz_edge = (pi / xtok);
            obj.y_bz_edge = (pi / xtok);

            % get coordinates of pictures
            img_size = obj.data_folder_stack(1).ImgCropSize;
            [xx, yy] = meshgrid( 1 : img_size, 1 : img_size);
            xx = xx - obj.cx;
            yy = yy - obj.cy;
            
            % k-vectors.
            % x(T/4) = hbar/(m*omega) * K(0)
            % We typically do not use the coordinate x, but rather x/a
            % or coordinates k, but rather (k*a)
            % therefore, (k*a) = m*omega*a / hbar * (x/a) * a
            kxx = xtok * xx;
            kyy = xtok * yy;
            
            % get roi nearest to BZ with the stipulation that the number of
            % pixels is divisible by an odd number of bins
            if strcmp(obj.analysis_mode, 'round')
                [x_roi, x_size_binned, y_roi, y_size_binned] = ...
                    get_nearest_roi([obj.cx - obj.x_bz_edge, obj.cx + obj.x_bz_edge],...
                                    [obj.cy - obj.y_bz_edge, obj.cy + obj.y_bz_edge], obj.img_bin_size, 'odd');
                dx = (x_roi(2) - x_roi(1) ) / (x_size_binned - 1);
                dy = (y_roi(2) - y_roi(1) ) / (y_size_binned - 1);
                
            elseif strcmp(obj.analysis_mode, 'resample')
                x_size_binned = obj.nsample_pixels;
                y_size_binned = obj.nsample_pixels;
                
                dx = 2 * obj.x_bz_edge / (x_size_binned - 1);
                dy = 2 * obj.y_bz_edge / (y_size_binned - 1);
                x_roi = [obj.cx - 0.5 * x_size_binned * dx, obj.cx + 0.5 * x_size_binned * dx];
                y_roi = [obj.cy - 0.5 * y_size_binned * dy, obj.cy + 0.5 * y_size_binned * dy];
            else
                error('obj.analysis mode was neither round nor resample. Incorrect setting');
            end
            
            % variables for storing results
            nr_full = zeros(y_size_binned, x_size_binned, length(obj.data_folder_stack));
            nr_full_unc = zeros(y_size_binned, x_size_binned, length(obj.data_folder_stack));
            nr = zeros(y_size_binned, x_size_binned, length(obj.data_folder_stack));
            nr_unc = zeros(size(nr));
            energies_2d = zeros(size(nr));
                      
            % momentum distributions
            for ii = 1 : length(obj.data_folder_stack)
%                 index = obj.sorted_index_frqs(ii);

                % crop, square average, and bin images
                img_full = obj.data_folder_stack(ii).Occs_ImgAvg;
                unc_full = obj.data_folder_stack(ii).Occs_ImgAvgUnc;
                
                if obj.use_background
                    img_full = img_full - obj.bg_folder.Occs_ImgAvg;
                    unc_full = sqrt( unc_full .^2 + obj.bg_folder.Occs_ImgAvgUnc .^ 2);
                end
                
                % choice of two analysis types. (1) crop image as closely
                % to Brillouin zone as possible, or (2) resample in
                % Brillouin zone.
                if strcmp(obj.analysis_mode, 'round')
                    img_crop = img_full(y_roi(1) : y_roi(2), x_roi(1) : x_roi(2));
                    unc_crop = unc_full(y_roi(1) : y_roi(2), x_roi(1) : x_roi(2)); 
                    
                    [img_bin, ~, ~, unc_bin] = binImg(img_crop, obj.img_bin_size, obj.img_bin_size, 'Normal', unc_crop);
                    % to preserve interpetation as a density we need to
                    % average our bins instead of adding them
%                     img_bin = (obj.img_bin_size)^2 * img_bin;
%                     unc_bin = (obj.img_bin_size)^2 * unc_bin;
                    obj.dx = dx * obj.img_bin_size;
                    obj.dy = dy * obj.img_bin_size;
                    
                    kxx_bin = binImg(kxx(y_roi(1) : y_roi(2), x_roi(1) : x_roi(2)),...
                                           obj.img_bin_size, obj.img_bin_size, 'Normal');
                    kyy_bin = binImg(kyy(y_roi(1) : y_roi(2), x_roi(1) : x_roi(2)),...
                                           obj.img_bin_size, obj.img_bin_size, 'Normal');
                elseif strcmp(obj.analysis_mode, 'resample')
%                     [img_bin, obj.dx, obj.dy, ~, unc_bin] = get_resampled_img(img_full, [x_size_binned, y_size_binned], x_roi, y_roi, 'mean', unc_full);
                    [img_bin, obj.dx, obj.dy, ~, ~, unc_bin] = get_resampled_img(img_full, [x_size_binned, y_size_binned], x_roi, y_roi, 'mean', unc_full);

                    kxx_bin = get_resampled_img(kxx, [x_size_binned, y_size_binned], x_roi, y_roi, 'mean');
                    kyy_bin = get_resampled_img(kyy, [x_size_binned, y_size_binned], x_roi, y_roi, 'mean');
                    % in this case my kxx and kyy's don't quite match the kx and
                    % ky values for the DQMC. I think this is a sampling
                    % issue?
                else 
                    error();
                end
                obj.dkx = xtok * obj.dx;
                obj.dky = xtok * obj.dy;
               
                
                % cropped image averaged using symmetry of square
                if strcmp(obj.symmetry, 'd4')
%                     [img_symm_avg, ~, ~, img_symm_unc] = d4_average(img_bin, unc_bin);
                    [img_symm_avg, img_symm_unc, ~, ~] = d4_average(img_bin, unc_bin);
                    [kxx_symm, ~, kxx_symm_std, ~, kyy_symm, ~, kyy_symm_std, ~] = d4_avg_vectorfn(kxx_bin, '', kyy_bin, '');
                elseif strcmp(obj.symmetry, 'd2')
%                     [img_symm_avg, ~, ~, img_symm_unc] = d2_average(img_bin, unc_bin);
                    [img_symm_avg, img_symm_unc, ~, ~] = d2_average(img_bin, unc_bin);
                    [kxx_symm, ~, kxx_symm_std, ~, kyy_symm, ~, kyy_symm_std, ~] = d2_avg_vectorfn(kxx_bin, '', kyy_bin, '');
                elseif strcmp(obj.symmetry, 'none')
                else
                    error('unsupported mode for obj.symmetry');
                end

                % bin image and store results
                 nr_full(:, :, ii) = img_bin;
                 nr_full_unc(:, :, ii) = unc_bin;
                 nr(:, :, ii) = img_symm_avg;
                 nr_unc(:, :, ii) = img_symm_unc;
                  
               if obj.init_hfstate_lower_energy
                    % if our initial hyperfine state is lower in energy
                    % than the final state, this is the correct expression.
                    % In this case, are we doing a |1> -> |2> rf transition
                   sign = -1;
               else
                   sign = 1;
                   % if our initial hyperfine state is higher in energy
                    % than the final state, this is the correct expression.
                    % In this case, are we doing a |3> -> |2> rf transition
               end
                 
               energies_2d(:, :, ii) = obj.latt_dispersion(kxx_symm, kyy_symm) + ...
                                                sign * obj.frqs_offsets_tunits(ii);
            end

            obj.specfn_unsymm = spectralfn(kxx_bin(1, :), kyy_bin(:, 1), energies_2d,...
                                '', '', '', '',...
                                '', nr_full, '', nr_full_unc, 'full');
            obj.specfn_uneq_es = spectralfn(kxx_symm(1, :), kyy_symm(:, 1), energies_2d,...
                                '', '', '', '',...
                                '', nr, '', nr_unc, 'full');

           % fit each Energy Distribution Curves (EDCS), i..e k-slices,  to
           % lorentzians
           % define fitting order
           dist_fn = @(x1, y1, x2, y2) min( mod( x1 - x2, 3), mod(x2 - x1, 3) );
           order = dist_fn(obj.specfn_uneq_es.linear_index, 0, 0, 0);
           % initial parameters
           init_params = [mean(obj.frqs_offsets_tunits), 1, max(obj.specfn_uneq_es.af_gxmg(1, :)), 0];
           fixed_params = [0, 0, 0, 1];
           lbs = [min(obj.frqs_offsets_tunits), 0.5, 0, -inf];
           ubs = [max(obj.frqs_offsets_tunits), inf, inf, inf];
           
           
           if obj.use_double_lor_dispersion
               n_lors = 2;
               init_params = horzcat(init_params, init_params);
               fixed_params = horzcat(fixed_params, fixed_params);
               lbs = horzcat(lbs, lbs);
               ubs = horzcat(ubs, ubs);
           else
               n_lors = 1;
           end
           
           % do fitting
           [obj.lor_fit_params_gxmg, obj.lor_fit_params_uncs_gxmg, ~] = ...
           fit_array_2lor(obj.specfn_uneq_es.af_gxmg, '', obj.frqs_offsets_tunits, obj.specfn_uneq_es.linear_index, '',...
                          init_params, fixed_params, lbs, ubs, n_lors, 1, dist_fn, order);
            
           if ~obj.use_double_lor_dispersion
               obj.lor_fit_params_gxmg = horzcat( obj.lor_fit_params_gxmg, zeros( size(obj.lor_fit_params_gxmg) ) );
               obj.lor_fit_params_uncs_gxmg = horzcat( obj.lor_fit_params_uncs_gxmg, zeros( size(obj.lor_fit_params_uncs_gxmg) ) );
           end
           % alternative way to get qp dispersion
%            [t, te, ~] = ...
%            fit_array_2lor(obj.specfn_uneq_es.af_gxmg, '', obj.specfn_uneq_es.es_gxmg, obj.specfn_uneq_es.linear_index, '',...
%                           init_params, fixed_params, lbs, ubs, n_lors, 1, dist_fn, order);


           % if using d2 symmetry, also look at gymg line. If using d4
           % symmetry, it will be identical to d2.
           if strcmp(obj.symmetry, 'd2')
                [obj.lor_fit_params_gymg, obj.lor_fit_params_uncs_gymg, ~] = ...
                 fit_array_2lor(obj.specfn_uneq_es.af_gymg, '', obj.frqs_offsets_tunits, obj.specfn_uneq_es.linear_index, '',...
                                init_params, fixed_params, lbs, ubs, n_lors, 1, dist_fn, order);
           else
               obj.lor_fit_params_gymg = obj.lor_fit_params_gxmg;
               obj.lor_fit_params_uncs_gymg = obj.lor_fit_params_uncs_gxmg;
           end
                 
            % bin energies in preparation for computing spectral function
            % A(k,w)*f(w)
            obj.energy_bin_edges = linspace( min(obj.specfn_uneq_es.es(:)) - 0.2,...
                                             max(obj.specfn_uneq_es.es(:)) + 0.2,...
                                             obj.n_energy_bins + 1);
            obj.energy_bin_means = 0.5 * (obj.energy_bin_edges(1:end-1) + obj.energy_bin_edges(2:end));
            
            % spectral function over full brillouin zone
            n_pts = zeros([length(obj.specfn_uneq_es.kys), length(obj.specfn_uneq_es.kxs), obj.n_energy_bins]);
            af_ebinned = zeros( size(n_pts) );
            af_ebinned_unc = zeros( size(n_pts) );

            % loop over all k-vectors and accumulate energies
            for ii = 1 : length(obj.specfn_uneq_es.kys)
                for jj = 1 : length(obj.specfn_uneq_es.kxs)
                     [~, ~, af_ebinned(ii, jj, :), ~, af_ebinned_unc(ii, jj, :), n_pts(ii, jj, :), ~]...
                         = azAvg_General(obj.specfn_uneq_es.af(ii, jj, :),...
                                         1 ./ obj.specfn_uneq_es.afunc(ii, jj, :).^2,...
                                         obj.specfn_uneq_es.es(ii, jj, :),...
                                         obj.energy_bin_edges);
%                      get_resampled_img(img, num_pix, x_roi, y_roi, mode, img_unc);
                end
            end
            
            obj.specfn = spectralfn( kxx_symm(1, :), kyy_symm(:, 1),...
                                     obj.energy_bin_means, '', '', '', '',...
                                     '', af_ebinned, '', af_ebinned_unc, 'full');         
            
            % fit \sum_k A(k, w) * f(w) to a lorentzian
            obj.specfn_sumk = squeeze(nansum( nansum( obj.specfn.af, 1), 2));
            
            init_params = [mean(obj.specfn.es), 2, max(obj.specfn_sumk), 0];
            fixed_params = [0, 0, 0, 1];
            lbs = [min(obj.specfn.es), 0, 0, -inf];
            ubs = [max(obj.specfn.es), inf, inf, inf];
       
            [obj.specfn_sumk_fitp, ~, ~, obj.specfn_sumk_fitp_err, ~] = ...
                fit1D(obj.specfn.es, obj.specfn_sumk, [], {'lorentzian1D'},...
                      init_params, fixed_params, lbs, ubs);
            
            % fit \sum_k A(k, ek - rf) * f(ek - rf) to a lorentzian
            obj.rfxfer_sumk = squeeze(sum( sum( obj.specfn_uneq_es.af, 2), 1));
            
            init_params = [mean(obj.frqs_offsets_tunits), 2, max(obj.rfxfer_sumk), 0];
            fixed_params = [0, 0, 0, 1];
            lbs = [min(obj.frqs_offsets_tunits), 0, 0, -inf];
            ubs = [max(obj.frqs_offsets_tunits), inf, inf, inf];
       
            [obj.rfxfer_sumk_fitp, ~, ~, obj.rfxfer_sumk_fitp_err, ~] = ...
                fit1D(obj.frqs_offsets_tunits, obj.rfxfer_sumk, [], {'lorentzian1D'},...
                      init_params, fixed_params, lbs, ubs);                                    
           
           % quasiparticle dispersion
%            obj.qp_dispersion_af_gxmg = horzcat(...
%                obj.latt_dispersion(obj.specfn_uneq_es.kxs_gxmg, obj.specfn_uneq_es.kys_gxmg) + sign * obj.lor_fit_params_gxmg(:, 1),...
%                obj.latt_dispersion(obj.specfn_uneq_es.kxs_gxmg, obj.specfn_uneq_es.kys_gxmg) + sign * obj.lor_fit_params_gxmg(:, 5) );
%            
%            obj.qp_dispersion_af_gxmg_unc = horzcat(...
%                obj.lor_fit_params_uncs_gxmg(:, 1), obj.lor_fit_params_uncs_gxmg(:, 5));
%            
%            if ~strcmp(obj.symmetry, 'd4')
%                obj.qp_dispersion_af_gymg = horzcat(...
%                    obj.latt_dispersion(obj.specfn_uneq_es.kxs_gxmg, obj.specfn_uneq_es.kys_gxmg) + sign * obj.lor_fit_params_gymg(:, 1),...
%                    obj.latt_dispersion(obj.specfn_uneq_es.kxs_gxmg, obj.specfn_uneq_es.kys_gxmg) + sign * obj.lor_fit_params_gymg(:, 5) );
%                obj.qp_dispersion_af_gymg_unc = horzcat(...
%                    obj.lor_fit_params_uncs_gymg(:, 1), obj.lor_fit_params_uncs_gymg(:, 1));  
%            else
%                obj.qp_dispersion_af_gymg = obj.qp_dispersion_af_gxmg;
%                obj.qp_dispersion_af_gymg_unc = obj.qp_dispersion_af_gxmg_unc;  
%            end
               
           % fit qp dispersions to BCS dispersion
%             [obj.fit_params_bcs_gxmg, obj.stderr_bcs_gxmg,...
%                 obj.linear_index_bcs_gxmg, obj.bcs_fn_sample_gxmg] = ...
%                 obj.fit_bcs_dispersion(obj.qp_dispersion_af_gxmg(:, 1), obj.qp_dispersion_af_gxmg_unc(:, 1));
%             
%             if ~strcmp(obj.symmetry, 'd4')
%                 [obj.fit_params_bcs_gymg, obj.stderr_bcs_gymg,...
%                     obj.linear_index_bcs_gymg, obj.bcs_fn_sample_gymg] = ...
%                     obj.fit_bcs_dispersion(obj.qp_dispersion_af_gymg(:, 1), obj.qp_dispersion_af_gymg_unc(:, 1));
%             else
%                 obj.fit_params_bcs_gymg = obj.fit_params_bcs_gxmg;
%                 obj.stderr_bcs_gymg = obj.stderr_bcs_gxmg;
%                 obj.linear_index_bcs_gymg = obj.linear_index_bcs_gxmg;
%                 obj.bcs_fn_sample_gymg = obj.bcs_fn_sample_gxmg;
%             end
                      
        end
       
        function fit_temp_dqmc(obj, init_params, fixed_params, use_to_fit,...
                               density_limits, static_dqmc_path)
            % fit static quantities to DQMC to determine U, T, central
            % chemical potential, etc.
            %
            % init_params: [U, T, h, detection efficiency, doublon transfer efficiency]
            %
            % fixed_params:
            %
            % use_to_fit: []
            %
            % density_limits: vector of length two defining minimum and
            % maximum density values to be included in fit
            %
            % static_dqmc_path: 
                           
            if ~exist('init_params', 'var') || isempty(init_params)
                init_params = [-4.5, 0.5, 0, 0.97, 0.99];
            end
            
            if ~exist('fixed_params', 'var') || isempty(fixed_params)
                fixed_params = [0, 0, 1, 1, 1];
            end

            if ~exist('use_to_fit', 'var') || isempty(use_to_fit)
                use_to_fit = [];
            end
            
            if ~exist('density_limits', 'var') || isempty(density_limits)
                density_limits = [0.1, 2];
            end
            
            if ~exist('static_dqmc_path', 'var') || isempty(static_dqmc_path)
                static_dqmc_path = fullfile('\\128.112.86.75\lithium\Publications',...
                    '06_arpes_and_att_hubb_spin_suc', 'dqmc',...
                'dqmc_nsites=8x8_T=0.3-10.0_U=-14.5--2.0_muavg=-10.0-8.0_mudelta=0.0-0.0_npass=0-20000_interp_fns.mat');
            end
            
            obj.static_dqmc_path = static_dqmc_path;
            
            dqmc = load(obj.static_dqmc_path);
            [~, ~, redchi_sqr, obj.mu_struct, obj.temp_fit_struct] = ...
                fit_dqmc_attractive(obj.ones, obj.threes, obj.singles, obj.doubles,...
                                    use_to_fit, density_limits, init_params, fixed_params, dqmc);

            % offset mu_center. Using convention that zero energy is at the
            % bottom of the tight-binding dispersion (which accounts for
            % the 4). 0.5*U is the chemical potential at half-filling when
            % writing out the Hubbard model as usual. The QuestQMC uses the
            % convention that mu=0 at half-filling
            obj.mu_center = obj.mu_struct.mu_peak + 4 + 0.5 * obj.temp_fit_struct.fit_params(1);

        end
           
        function [cx, cy] = find_center(obj, folder_stack)
            % find average center for all pictures
            
            nimgs = vertcat(folder_stack.nimgs);
            avg_img_stack = cat(3, folder_stack.Occs_ImgAvg);
            avg_img = 0;
            for ii = 1:length(folder_stack)
                avg_img = avg_img + nimgs(ii) * avg_img_stack(:, :, ii);
            end
            avg_img = avg_img / sum(nimgs);
            
            [cx, cy, ~] = get_moment(avg_img, 1);
            
        end
           
        function [fit_params_bcs, stderr_bcs, linear_index_bcs, bcs_fn_sample] = ...
                fit_bcs_dispersion(obj, qp_dispersion, qp_dispersion_err)
            % fit to BCS dispersion. Everything in units of t/h
            %
            
            if any(isnan(qp_dispersion)) || any(isnan(qp_dispersion_err))
                warning('NaNs in qp_dispersion or qp_dispersion_err passed to fit_bcs_dispersion');
            end
           
           indices_to_use = ~isnan(qp_dispersion) & ~isnan(qp_dispersion_err);
           
           kxs_gxmg = obj.specfn_uneq_es.kxs_gxmg(indices_to_use);
           kys_gxmg = obj.specfn_uneq_es.kys_gxmg(indices_to_use);
           qp_dispersion = qp_dispersion(indices_to_use);
           qp_dispersion_err = qp_dispersion_err(indices_to_use); 
           
           bcs_disp_fn = @(kx, ky, P) P(3) - sqrt( ( obj.latt_dispersion(kx, ky) - P(1) ).^2 + P(2)^2 );
           bcs_min_fn = @(P) (qp_dispersion - bcs_disp_fn(kxs_gxmg, kys_gxmg, P) ) ./ qp_dispersion_err;
           
           init_params = [4, 1, 0];
           fixed_params = [0, 0, 0];
           lbs = [0, 0, -inf];
           ubs = [inf, inf, inf];
           [fit_params_bcs, stderr_bcs, ~, ~] = lsq_fixedp(bcs_min_fn, init_params, fixed_params, lbs, ubs);

           linear_index_bcs = linspace(0, 3, 300);
           kxs_interp = interp1(obj.specfn_uneq_es.linear_index, obj.specfn_uneq_es.kxs_gxmg, linear_index_bcs);
           kys_interp = interp1(obj.specfn_uneq_es.linear_index, obj.specfn_uneq_es.kys_gxmg, linear_index_bcs);
           bcs_fn_sample = bcs_disp_fn(kxs_interp, kys_interp, fit_params_bcs);
        end
          
        function [lor_fit_params, lor_fit_uncs] = ...
                fit_dispersion(obj, frqs_tunits, linear_index, nr_vs_frq, nr_vs_frq_err, use_double_lor)
             % Deprecated. TODO: remove
             %
             % fit 'dispersion relation' along the GXMG or GYMG path from rf images vs. k. i.e., fit
             % energy dispersion curves for center coordinates.
             % give results in terms of tunneling
             % Arguments
             % frequencies: are the RF frequencies at which the data were
             % taken. This function should work for any choice of units.
             %
             % nk_vs_frq: a #kvectors x #frequencies array. This function
             % fits each row of nk_vs_frq to a lorentzian.
             %
             % nk_vs_frq_err: are the uncertainties associated with
             %  nk_vs_frq
             %
              
             % TODO: deal with INFs and add uncertainty to fit
             %1 ./ nk_vs_frq_err(ii, obj.sorted_index_frqs) .^2,...
             
             if ~exist('use_double_lor', 'var')
                 use_double_lor = 0;
             end
             
             % When not using double lorentzian, simply fix those
             % parameters to zero and do not use in fitting.    
            init_params = [[mean(frqs_tunits), 1, max(nr_vs_frq(1, :)), 0], [0, 0, 0, 0]];
            fixed_params = [[0, 0, 0, 1], [1, 1, 1, 1]];
            ubs = [[max(frqs_tunits), 0.5 * (max(frqs_tunits) - min(frqs_tunits)), 3 * max(nr_vs_frq(:)),  max(nr_vs_frq(:))], [0, 0, 0, 0]];
            lbs = [[min(frqs_tunits), min(frqs_tunits(2:end) - frqs_tunits(1:end-1)), -0, -1], [0, 0, 0, 0]];

            lor_fit_params = zeros(size(nr_vs_frq, 1), 8);
            lor_fit_uncs = zeros(size(lor_fit_params));

            fn_name = {'lorentzian1D', 'lorentzian1D'};
            
            if use_double_lor
                init_params(5:8) = init_params(1:4);
                fixed_params(5:8) = fixed_params(1:4);
                ubs(5:8) = ubs(1:4);
                lbs(5:8) = lbs(1:4);
            end
            
            % first fit along the GXM path
            indices_GXM = find(linear_index < 2 & linear_index >=0);
            indices_GM = flip( find(linear_index <=3 & linear_index >=2) );
            
            for ii = 1 : length(indices_GXM)
                index = indices_GXM(ii);
                
                % initial guess from previous fit
                if ii == 1
                    guess = init_params;
                else
                    guess = lor_fit_params(indices_GXM(ii - 1) , :);
                end
                
                 [lor_fit_params(index, :), ~, ~, lor_fit_uncs(index, :), ~] = ...
                        fit1D(frqs_tunits, nr_vs_frq(index, :), [], fn_name, guess, fixed_params, lbs, ubs);
                    
                % ensure that the first lorentzian is at lower frequency
                % than the second
                if use_double_lor && lor_fit_params(index, 5) < lor_fit_params(index, 1)

                    helper_params = zeros( size(lor_fit_params(index, :)) );
                    helper_params(1:4) = lor_fit_params(index, 5:8);
                    helper_params(5:8) = lor_fit_params(index, 1:4);
                    lor_fit_params(index, :) = helper_params;

                    helper_err = zeros( size(helper_params) );
                    helper_err(1:4) = lor_fit_uncs(index, 5:8);
                    helper_err(5:8) = lor_fit_uncs(index, 1:4);
                    lor_fit_uncs(index, :) = helper_err;
                end
            end
            
            % Fit along the GM path
            for jj = 1 : length(indices_GM)
                index = indices_GM(jj);
                
                % initial guess from previous fit
                if jj == 1
                    guess = init_params;
                else
                    guess = lor_fit_params(indices_GM(jj - 1) , :);
                end
                
                % do fitting
                 [lor_fit_params(index, :), ~, ~, lor_fit_uncs(index, :), ~] = ...
                        fit1D(frqs_tunits, nr_vs_frq(index, :), [], fn_name, guess, fixed_params, lbs, ubs);
                    
                                        
                % ensure that the first lorentzian is at lower frequency
                % than the second
                if use_double_lor && lor_fit_params(index, 5) < lor_fit_params(index, 1)

                    helper_params = zeros( size(lor_fit_params(index, :)) );
                    helper_params(1:4) = lor_fit_params(index, 5:8);
                    helper_params(5:8) = lor_fit_params(index, 1:4);
                    lor_fit_params(index, :) = helper_params;

                    helper_err = zeros( size(helper_params) );
                    helper_err(1:4) = lor_fit_uncs(index, 5:8);
                    helper_err(5:8) = lor_fit_uncs(index, 1:4);
                    lor_fit_uncs(index, :) = helper_err;
                end
                    
            end

           
        end
       
        % display all and save
        function display_and_save_all(obj, saving)
            % display all plots and save all results. This function uses
            % the individual plot functions below.
            
            if ~exist('saving', 'var') || isempty(saving)
                saving = 1;
            end
            
            fig_handle_dqmc_fit = obj.display_temp_fit();
            fig_handle_nonintrf = obj.plot_rf_nonint();
            fig_handle_kspace = obj.show_momentum_distributions(); 
            fig_handle_lor_fits_gxmg = obj.show_lor_fits_dispersion('gxmg');
            fig_handle_lor_fits_gymg = obj.show_lor_fits_dispersion('gymg');
%             fig_handle_nr_vs_omega_allks = obj.show_full_spectral_fn();
            fig_handle_arpes = obj.show_summary(0);  
            
           fig_handle_spec_fn = obj.show_af_gxmg_expt('', 0);
           suptitle(strrep(obj.identifier, '_', ' '));
            
            % arpes theory comparisons
            fig_handle_thry_compare_gxmg = obj.show_edc_vs_thry('gxmg');
            fig_handle_thry_compare_gymg = obj.show_edc_vs_thry('gymg');
            fig_handle_thry_compare_af2d = obj.compare_af_gxmg_thry_expt();

            if saving

                % create directory for saving
                save_dir = obj.identifier;

                if ~exist(save_dir, 'dir')
                    mkdir(save_dir);
                end

                % non-interacting rf
                if ~isempty(obj.rfspec_df)
                    fname = sprintf('%s_rfspectrum.fig', obj.rfspec_df.DatasetString);
                    fname = fullfile(save_dir, fname);
                    savefig(fig_handle_nonintrf, fname);
                end

                % temperature figure
                % TODO: because of how this function is written, production of this
                % figure is handled by the same function which does the fitting.
                if exist('fig_handle_dqmc_fit', 'var') && isvalid(fig_handle_dqmc_fit)
                    fname = sprintf('%s_dqmc_fit.fig', obj.identifier);
                    fname = fullfile(save_dir, fname);
                    savefig(fig_handle_dqmc_fit, fname);
                end

                % momentum space maps
                fname = sprintf('%s_momentum_distributions.fig', obj.identifier);
                fname = fullfile(save_dir, fname);
                savefig(fig_handle_kspace, fname);

                % lorentzian fits to rf xfer 
                fname = sprintf('%s_lorentzian_fits_gxmg.fig', obj.identifier);
                fname = fullfile(save_dir, fname);
                savefig(fig_handle_lor_fits_gxmg, fname);

                fname = sprintf('%s_lorentzian_fits_gymg.fig', obj.identifier);
                fname = fullfile(save_dir, fname);
                savefig(fig_handle_lor_fits_gymg, fname);

                fname = sprintf('%s_specfn_thry_vs_exp_gxmg.fig', obj.identifier);
                fname = fullfile(save_dir, fname);
                savefig(fig_handle_thry_compare_gxmg, fname);
                
                fname = sprintf('%s_specfn_thry_vs_exp_gymg.fig', obj.identifier);
                fname = fullfile(save_dir, fname);
                savefig(fig_handle_thry_compare_gymg, fname);
                
%                 fname = sprintf('%s_nr_allks.fig', obj.identifier);
%                 fname = fullfile(save_dir, fname);
%                 savefig(fig_handle_nr_vs_omega_allks, fname);

                fname = sprintf('%s_specfn_thry_vs_exp.fig', obj.identifier);
                fname = fullfile(save_dir, fname);
                savefig(fig_handle_thry_compare_af2d, fname);
                
                fname = sprintf('%s_specfn.fig', obj.identifier);
                fname = fullfile(save_dir, fname);
                savefig(fig_handle_spec_fn, fname);
                
                % spectral function             
                fname = sprintf('%s_summary.fig', obj.identifier);
                fname = fullfile(save_dir, fname);
                savefig(fig_handle_arpes, fname); 

                % save text files
                obj.export_expt2txt(save_dir);
                obj.export_thry2txt(save_dir);
                
                % save structure and folders
                obj.saveStruct(save_dir, '', 1);
                obj.saveAllFolders(save_dir);
            end
        end
        
        % display different properties
        function fig_handle = display_temp_fit(obj)
            % show temperature fit figure
            
            fig_handle = display_dqmc_fit(obj.temp_fit_struct, obj.mu_struct);
            
            % get suptitle and add our identifier string to it
            tags = cell(length(fig_handle.Children), 1);
            for ii = 1:length(fig_handle.Children)
                tags{ii} = fig_handle.Children(ii).Tag;
            end
            fn = @(c) strcmp(c, 'suptitle');
            index = find(cellfun(fn, tags));
            fig_handle.Children(index).Children.String = ...
                vertcat(obj.identifier, fig_handle.Children(index).Children.String);
        end
        
        function fig_handle = show_summary(obj, plot_log, fig_size_pixels)
            % Create summary figure for experimental results
            
            if ~exist('plot_log', 'var') || isempty(plot_log)
                plot_log = 1;
            end
            
            if ~exist('fig_size_pixels', 'var')
                a = get(0, 'Screensize');
                fig_size_pixels = a(3:4);
            end       
            
             % 3d plot of A(k, w)
            fig_handle = figure();
            fig_handle.Position = horzcat([0, 0], fig_size_pixels);
            
            nrows = 3;
            ncols = 3;

            % energy distribution curves, i.e. A(k=c, w) * f(w)
            axes_h = subplot(nrows, ncols, [4,1]);
            obj.show_edcs(axes_h);
            
            % density plot of A(k, w)
            axes_h = subplot(nrows, ncols, 2);
            obj.show_af_gxmg_expt(axes_h, plot_log);
            
            % spectral fn, k summed 
            axes_h = subplot(nrows, ncols, 5);
            obj.show_specfn_ksummed(axes_h);

            % density plot of A(k, e_k - rf_frq)
            axes_h = subplot(nrows, ncols, 3);
            obj.show_rfxfer_fn(axes_h, plot_log);

            % momentum unresolved and resolved rf spectra
            axes_h = subplot(nrows, ncols, 6);
            obj.show_rfxfer_ksummed(axes_h);

            % fit centers of rf spectrum at each momentum
            axes_h = subplot(nrows, ncols, 9);
            obj.show_rfcenters(axes_h);  
        
            % show quasiparticle dispersion
            axes_h = subplot(nrows, ncols, 8);
            obj.show_qpdispersions(axes_h);

            % plot HWHM
            axes_h = subplot(nrows, ncols, 7);
            obj.show_rfhwhm(axes_h);        
            
            ident = strrep(obj.identifier, '_', ' ');
            suptitle(sprintf('%s\n U/t = %0.2f(%.0f), T/t = %0.2f(%.0f), t/h = %.0f Hz\n RF = %.0fdBm at %0.1fms latt = %0.2fv',...
                     ident,...
                     obj.temp_fit_struct.fit_params(1), obj.temp_fit_struct.std_errs(1) * 1e2,...
                     obj.temp_fit_struct.fit_params(2), obj.temp_fit_struct.std_errs(2) * 1e2, obj.t_over_h,...
                     obj.rf_power, obj.rf_time, obj.latt_depth));
            
        end
        
        function [fig_handle, axes] = show_af_gxmg_thry(obj, axes, plot_log)
            % show an image in k - E space fo the theory spectral function
            % along the GXMG cut. Also show the non-interacting dispersion
            % and the chemical potential
            
            if ~exist('axes', 'var') || isempty(axes) || ~isvalid(axes)
                fig_handle = figure;
                axes = gca;
            else
                fig_handle = gcf;
            end
            
            if ~exist('plot_log', 'var') || isempty(plot_log)
                plot_log = 0;
            end
            
            if ~isempty(obj.dqmc_aavg.af_gxmg)
                % crop theory to experiment
                omegas = obj.dqmc_aavg.es;
                I = omegas <= max(obj.energy_bin_edges) & omegas >= min(obj.energy_bin_edges); 

                omegas = omegas(I);
                spectral_fn = squeeze( obj.dqmc_aavg.af_gxmg(:, I) );

                % expand theory to experiment
                if max(obj.energy_bin_edges) > max(omegas)
                    dw = omegas(2) - omegas(1);
                    pad_pix = floor( (max(obj.energy_bin_edges) - max(omegas)) / dw );

                    omegas_extra = max(omegas) : dw : max(obj.energy_bin_edges);
                    omegas = horzcat(omegas, omegas_extra(2:end) );
                    spectral_fn = horzcat(spectral_fn, nan(size(spectral_fn, 1), pad_pix) );
                end

                if min(obj.energy_bin_edges) < min(omegas)
                    dw = omegas(2) - omegas(1);
                    pad_pix = floor( ( min(omegas) - min(obj.energy_bin_edges)) / dw );

                    omegas_extra = min(omegas) : -dw : min(obj.energy_bin_edges);
                    omegas_extra = flip(omegas_extra(2:end) );
                    omegas = horzcat(omegas_extra, omegas);
                    spectral_fn = horzcat(nan(size(spectral_fn, 1), pad_pix), spectral_fn);
                end

                if plot_log
                    imagesc(axes, obj.dqmc_aavg.linear_index, omegas, log10(spectral_fn'));
                else
                    imagesc(axes, obj.dqmc_aavg.linear_index, omegas, spectral_fn');
                end
                hold on;
                for jj = 1 : 4
                    plot(axes, obj.dqmc_aavg.linear_index, obj.dqmc_af_dispersion_gxmg_fitp(:, 1 + (jj-1) * 4),...
                         'bx', 'LineWidth', 2);
                    plot(axes, obj.dqmc_aavg.linear_index, obj.dqmc_a_dispersion_gxmg_fitp(:, 1 + (jj-1) * 4),...
                     'gx', 'LineWidth', 2);
                end
                colormap hot;
                hold(axes, 'on');
                axes.YDir = 'normal';

             
                % plot non-interacting dispersion
                interp_indices = linspace(0, 3, 300);
                [kx_interp, ky_interp] = index2kvect(interp_indices);
                 plot(axes, interp_indices, obj.latt_dispersion( kx_interp, ky_interp),...
                            'b', 'LineWidth', 2);   

                % plot chemical potential
                if ~isempty(obj.mu_center)         
                    plot(axes, [0, 3], [obj.mu_center, obj.mu_center],...
                        'b', 'LineWidth', 2);
                end

                % x -axis
                axes.XTick = [0, 1, 2, 3];
                axes.XTickLabel = {'\Gamma', 'X', 'M', '\Gamma'};
                xlabel('Momentum, K');

                % y-axis
                axes.YTick = floor(min(omegas)) : 4 : ceil(max(omegas));
            end
            
            ylabel('Energy (t)');
            
            title(axes, 'A(k, w)*f(w)');
        end
        
        function fig_handle = compare_af_gxmg_thry_expt(obj, plot_log, fig_size_pixels)
            % compare experimental spectral function with DQMC trap
            % averaged spectral function
            
           if ~exist('plot_log', 'var') || isempty(plot_log)
               plot_log = 0;
           end
           
           if ~exist('fig_size_pixels', 'var')
               a = get(0, 'Screensize');
               fig_size_pixels = a(3:4);
           end                
            
           fig_handle = figure;
           fig_handle.Position = horzcat([0, 0], fig_size_pixels);
           
           ax1 = subplot(1, 2, 1);
           ax2 = subplot(1, 2, 2);
           
           obj.show_af_gxmg_expt(ax1, plot_log, 1);
           obj.show_af_gxmg_thry(ax2, plot_log);
           
           ax2.CLim = ax1.CLim;
           
           temp = obj.temp_fit_struct.fit_params(2);
           temp_err = obj.temp_fit_struct.std_errs(2);
           U = obj.temp_fit_struct.fit_params(1);
           U_err = obj.temp_fit_struct.std_errs(1);
           
           ax1.Title.String = sprintf('A(k, w) * f(w), expt\n U/t=%0.2f(%.0f), T/t=%0.2f(%.0f)', U, U_err * 1e2, temp, temp_err * 1e2);
           ax2.Title.String = sprintf('A(k, w) * f(w), thry\n U/t=%0.2f, T/t=%0.2f', obj.dqmc_aavg.us, obj.dqmc_aavg.ts);
           
           suptitle(strrep(obj.identifier, '_', ' '));
            
        end
        
        function [fig_handle] = show_momentum_distributions(obj, fig_size_pixels, img_scale)
                % display 2d momentum distributions: raw pictures obtained
                % from experiment
            
                if ~exist('fig_size_pixels', 'var') || isempty(fig_size_pixels)
                    a = get(0, 'Screensize');
                    fig_size_pixels = a(3:4);
                end                
            
                if ~exist('img_scale', 'var') || isempty(img_scale)
                    img_scale = 0.1;
                end
                
                % plot momentum distributions
                fig_handle = figure;
                fig_handle.Position = horzcat([0, 0], fig_size_pixels);
                
                nrows = floor(sqrt(length(obj.data_folder_stack)));
                ncols = ceil(length(obj.data_folder_stack) / nrows);

                % expand columns and rows
                ncols = 2 * ncols;
                nrows = 2 * nrows;
                
                for ii = 1 : length(obj.data_folder_stack)
                    % raw momentum distributions
                    df = obj.data_folder_stack(ii);
                    subplot(nrows, ncols, (2 * ii - 1) +  ncols * floor((ii-1) / (ncols/2)) );
                    imagesc(df.Occs_ImgAvg, [0, img_scale]);
                    axis equal;
                    axis image;
                    hold on;
                    scatter(obj.cx, obj.cy, 'rx');
                    ax = gca;
                    ax.XTickLabel = '';
                    ax.YTickLabel = '';
                    
                    % brillouin zone edge
                    hold on;
                    scatter(obj.cx, obj.cy, 'rx');
                    scatter(obj.cx + obj.x_bz_edge, obj.cy, 'rx');
                    scatter(obj.cx - obj.x_bz_edge, obj.cy, 'rx');
                    scatter(obj.cx, obj.cy + obj.x_bz_edge, 'rx');
                    scatter(obj.cx, obj.cy - obj.x_bz_edge, 'rx');
                    scatter(obj.cx + obj.x_bz_edge, obj.cy + obj.x_bz_edge, 'rx');
                    scatter(obj.cx + obj.x_bz_edge, obj.cy - obj.x_bz_edge, 'rx');
                    scatter(obj.cx - obj.x_bz_edge, obj.cy + obj.x_bz_edge, 'rx');
                    scatter(obj.cx - obj.x_bz_edge, obj.cy - obj.x_bz_edge, 'rx');
                    title(sprintf('frq = %0.4f MHz, folder %03d', obj.frqs_hz_sorted(ii) / 1e6, df.Dataset{4}));

%                     title(sprintf('frq = %0.4f MHz, folder %03d', obj.frqs_hz(index) / 1e6, obj.folders(index)));

                    % plot momentum distributions after d4 averaging
                    subplot(nrows, ncols, (2 * ii) +  ncols * floor((ii-1) / (ncols/2)) )
                    imagesc(obj.specfn_uneq_es.af(:, :, ii), [0, img_scale]);
                    axis equal;
                    axis image;
                    ax = gca;
                    ax.XTick = [1, size(obj.specfn_uneq_es.af, 2)/2 + 1 , size(obj.specfn_uneq_es.af, 2)];
                    ax.XTickLabel = {'-\pi', '0', '\pi'};
                    ax.YTick = [1, size(obj.specfn_uneq_es.af, 1)/2 + 1 , size(obj.specfn_uneq_es.af, 1)];
                    ax.YTickLabel = {'\pi', '0', '-\pi'};

                    % plot along high symmetry directions gxmg and gymg
                    subplot(nrows, ncols, [(2 * ii - 1), 2*ii] +  ncols * (floor((ii-1) / (ncols/2) + 1)) )
                    errorbar(obj.specfn_uneq_es.linear_index, obj.specfn_uneq_es.af_gxmg(:, ii), obj.specfn_uneq_es.afunc_gxmg(:, ii), '.--');
                    hold on;
                    errorbar(obj.specfn_uneq_es.linear_index, obj.specfn_uneq_es.af_gymg(:, ii), obj.specfn_uneq_es.afunc_gymg(:, ii), '.--');
                    grid on

                    ax = gca;
                    xlim([0, 3]);
                    ax.XTick = [0, 1, 2, 3];
                    ax.XTickLabel = {'\Gamma', 'X', 'M', '\Gamma'};
                    ax.YLim(1) = -0.02;
                end
                
                suptitle(strrep(obj.identifier, '_', ' '));
            
        end
        
        function [fig_handle] = show_edc_vs_thry(obj, mode, fig_size_pixels)
            % show lorentzian fits of atom transfer vs. frequency
            % for each k-vector along the GXMG and GYMG cuts in the
            % Brillouin zone.
            
            if ~exist('mode', 'var')
                mode = 'gxmg';
            end
            
            if ~exist('fig_size_pixels', 'var')
                    a = get(0, 'Screensize');
                    fig_size_pixels = a(3:4);
            end                
                                          
            fig_handle = figure;
            fig_handle.Position = horzcat([0, 0], fig_size_pixels);
              % suptitle
            suptitle( sprintf('%s Lorentzian fits vs. frequency for momenta %s',...
                      strrep(obj.identifier, '_', ' '), mode) );
            
            nplots = length(obj.specfn_uneq_scaled.linear_index);
            nrows = floor( sqrt(nplots) );
            ncols = ceil( nplots / nrows );
                       
            ax_handles = [];
            ax2_handles = [];
            thry_titles = cell(1, nplots);
            for ii = 1 : nplots
                ax = subplot(nrows, ncols, ii);

                grid on;
                xlim( [min(obj.frqs_offsets_tunits), max(obj.frqs_offsets_tunits)] );
                
                % create second y-axis for energy
                ax2 = axes('Position', ax.Position);
                ax2.YLim = ax.YLim;
                
                % plot atom number transfer vs. energy
                if strcmp(mode, 'gxmg')
                    % plot spectral function at unequal energies
%                     errorbar(ax2, obj.specfn_uneq_es.es_gxmg(ii, :), obj.specfn_uneq_es.af_gxmg(ii, :),...
%                              obj.specfn_uneq_es.afunc_gxmg(ii, :), 'bo');
                    
                    errorbar(ax2, obj.specfn_uneq_scaled.es_gxmg(ii, :), obj.specfn_uneq_scaled.af_gxmg(ii, :),...
                             obj.specfn_uneq_scaled.afunc_gxmg(ii, :), 'bo');

                    % also plot spectral function. Should be the same as atom
                    % number transfer except for resampling to get at same
                    % energy points for every k.
                    hold on;
%                     errorbar(ax2, obj.energy_bin_means, obj.specfn.af_gxmg(ii, :),...
%                                   obj.specfn.afunc_gxmg(ii, :), 'ko');
                    errorbar(ax2, obj.specfn_scaled.es, obj.specfn_scaled.af_gxmg(ii, :),...
                                  obj.specfn_scaled.afunc_gxmg(ii, :), 'ko');
                
                    % plot zero
                    plot(ax2, [ min(obj.specfn_uneq_scaled.es_gxmg(ii, :)), max(obj.specfn_uneq_scaled.es_gxmg(ii, :)) ], [0, 0], 'k');
                    ax2.YLim = [min(obj.specfn_uneq_scaled.af_gxmg(ii, :)), max(obj.specfn_uneq_scaled.af_gxmg(ii, :)) + ...
                                                                            0.2 * abs(max(obj.specfn_uneq_scaled.af_gxmg(ii, :))) ];
                    
                elseif strcmp(mode, 'gymg')
%                     errorbar(ax2, obj.specfn_uneq_es.es_gymg(ii, :), obj.specfn_uneq_es.af_gymg(ii, :),...
%                              obj.specfn_uneq_es.afunc_gymg(ii, :), 'bo');
                    errorbar(ax2, obj.specfn_uneq_scaled.es_gymg(ii, :), obj.specfn_uneq_scaled.af_gymg(ii, :),...
                             obj.specfn_uneq_scaled.afunc_gymg(ii, :), 'bo');
                     hold on;
                     % also plot spectral function. Should be the same as atom
                    % number transfer except for resampling to get at same
                    % energy points for every k.
%                     errorbar(ax2, obj.energy_bin_means, obj.specfn.af_gymg(ii, :),...
%                                   obj.specfn.afunc_gymg(ii, :), 'k.');
                    errorbar(ax2, obj.specfn_scaled.es, obj.specfn_scaled.af_gymg(ii, :),...
                                  obj.specfn_scaled.afunc_gymg(ii, :), 'ko');

                    plot(ax2, [ min(obj.specfn_uneq_scaled.es_gymg(ii, :)), max(obj.specfn_uneq_scaled.es_gymg(ii, :)) ], [0, 0], 'k');
                    ax2.YLim = [min(obj.specfn_uneq_scaled.af_gymg(ii, :)),...
                                max(obj.specfn_uneq_scaled.af_gymg(ii, :)) + 0.2 * abs( max(obj.specfn_uneq_scaled.af_gymg(ii, :)) ) ];
                else
                    error();
                end
                hold on;
                       
                % if dqmc data is present, plot it too. Note that due to
                % symmetry no distinction between GXMG and GYMG here.
                thry_ttl = '';
                if ~isempty(obj.dqmc_aavg) 

                    % find nearest k-index
                    [~, min_index] = min( abs( obj.specfn_uneq_scaled.linear_index(ii) - obj.dqmc_aavg.linear_index) );
                    
                    % plot trap average
%                     thry_spec_fn = obj.lor_fit_params_gxmg(ii, 4) + obj.theory_scale_trpavg * squeeze( obj.dqmc_aavg.af_gxmg(min_index, :) ) + obj.theory_offset_trpavg;
%                     plot(ax2, obj.dqmc_aavg.es, thry_spec_fn, 'r' );
                    thry_spec_fn = squeeze(obj.dqmc_aavg.af_gxmg(min_index, :));
                    plot(ax2, obj.dqmc_aavg.es, thry_spec_fn , 'r' );
                    hold on;
                    ax2.YLim(1) = min( ax2.YLim(1), 0);
                    ax2.YLim(2) = max(ax2.YLim(2), 1.2 * max(thry_spec_fn(:)));    
                    
                    % title with k-vectors
                    [kx_curr, ky_curr] = index2kvect( obj.dqmc_aavg.linear_index(min_index) ); 
                    thry_ttl = sprintf('k th = (%0.2f,%0.2f)pi', kx_curr/pi, ky_curr/pi);
                    thry_titles{ii} = thry_ttl;
                end
                
                % plot mu
                if ~isempty(obj.mu_center)
                    plot(ax2, [obj.mu_center, obj.mu_center], [ax2.YLim(1), ax2.YLim(2)], ...
                            'LineWidth', 2);     
                end
                
                % set up second axis correctly relative to first.
                ax2.Color = 'none';
                if strcmp(mode, 'gxmg')
                    ax2.XLim = [ min(obj.specfn_uneq_scaled.es_gxmg(ii, :)), max(obj.specfn_uneq_scaled.es_gxmg(ii, :)) ];
                elseif strcmp(mode, 'gymg')
                    ax2.XLim = [ min(obj.specfn_uneq_scaled.es_gymg(ii, :)), max(obj.specfn_uneq_scaled.es_gymg(ii, :)) ];
                end
                    
                if obj.init_hfstate_lower_energy
                    % if the initial hyperfine state has lower energy than
                    % the final hyperfine state, E = e_k - nu, so the
                    % energy axis will be opposite to the nu axis (i.e.
                    % increasing rf => decreasing energy)
                    ax2.XDir = 'reverse';
                else
                    % if the initial hyperfine state has higher energy than
                    % the final hyperfine state, E = e_k + nu, so the
                    % energy axis will be the same direction as the nu axis
                    ax2.XDir = 'normal';
                end
                ax2.XAxisLocation = 'top';
                ax2.YAxisLocation = 'right';    
                ax2.YTickLabel = '';
                
                if ii == 1
                    legend({'expt', 'expt energy bins', '', 'thry'});
                end
                
                ax_handles = cat(1, ax_handles, ax);
                ax2_handles = cat(1, ax2_handles, ax2);
            end
            
            % add momentum labels as annotation
            for ii = 1 : nplots
                ax = ax_handles(ii);
                height = 0.05;
                sep = 0.0;
                dim = [ax.Position(1), ax.Position(2) + ax.Position(4) - height - sep,...
                       ax.Position(3), height];
                  
                if strcmp(mode, 'gxmg')
                    text = sprintf('k exp = (%0.2f,%0.2f)pi \n %s',...
                        obj.specfn_uneq_scaled.kxs_gxmg(ii) / pi, obj.specfn_uneq_scaled.kys_gxmg(ii) / pi, thry_titles{ii});
                elseif strcmp(mode, 'gymg')
                    text = sprintf('k exp = (%0.2f,%0.2f)pi \n %s',...
                        obj.specfn_uneq_scaled.kxs_gymg(ii) / pi, obj.specfn_uneq_scaled.kys_gymg(ii) / pi, thry_titles{ii});
                end
                
                annotation('textbox', dim, 'String', text, 'FitBoxToText', 'on', 'FontSize', 10);
            end
           
            for ii = 1 : nplots
                ax_handles(ii).YLim = ax2_handles(ii).YLim;
                ax_handles(ii).Position = ax2_handles(ii).Position;
            end
                        
        end
        
        function fig_handle = show_lor_fits_dispersion(obj, mode, fig_size_pixels)
            % show lorentzian fits to experimental EDCs
            
            if ~exist('mode', 'var')
                mode = 'gxmg';
            end
            
            if ~exist('fig_size_pixels', 'var')
                    a = get(0, 'Screensize');
                    fig_size_pixels = a(3:4);
            end                
            
            frqs_interp = linspace(min(obj.frqs_offsets_tunits),...
                                   max(obj.frqs_offsets_tunits), 300);
                               
            lor_fn = @(P, f) lorentzian1D(P(1:4), f) + lorentzian1D(P(5:8), f);
            
            fig_handle = figure;
            fig_handle.Position = horzcat([0, 0], fig_size_pixels);
            
            suptitle( sprintf('%s Lorentzian fits vs. frequency for momenta %s',...
                              strrep(obj.identifier, '_', ' '), mode) );
            
            nplots = size(obj.lor_fit_params_gxmg, 1);
            nrows = floor( sqrt(nplots) );
            ncols = ceil( nplots / nrows );
            
            ax_handles = [];
            ax2_handles = [];
            thry_titles = cell(1, nplots);
            for ii = 1 : nplots
                subplot(nrows, ncols, ii);
                
                if strcmp(mode, 'gxmg')
                    errorbar(obj.frqs_offsets_tunits, obj.specfn_uneq_es.af_gxmg(ii, :),...
                             obj.specfn_uneq_es.afunc_gxmg(ii, :), 'bo');
                    hold on;
                    plot(frqs_interp, lor_fn(obj.lor_fit_params_gymg(ii, :), frqs_interp), 'b' );
                elseif strcmp(mode, 'gymg')
                    errorbar(obj.frqs_offsets_tunits, obj.specfn_uneq_es.af_gymg(ii, :),...
                             obj.specfn_uneq_es.afunc_gymg(ii, :), 'bo');
                    hold on;
                    plot(frqs_interp, lor_fn(obj.lor_fit_params_gymg(ii, :), frqs_interp), 'b' );
                else
                    error();
                end

                grid on;
                xlabel('frq offsets');
                ylabel('rf transfer');
                
                % create second y-axis
                ax = gca;
                ax2 = axes('Position', ax.Position, 'XAxisLocation', 'top',...
                           'Xdir', 'rev',...
                           'YAxisLocation', 'right', 'Color', 'none');
                
                ax2.YLim = ax.YLim;
                ax2.Color = 'none';
                if strcmp(mode, 'gxmg')
                    ax2.XLim = [ min(obj.specfn_uneq_es.es_gxmg(ii, :)), max(obj.specfn_uneq_es.es_gxmg(ii, :)) ];
                                    title( sprintf('(kx,y) = (%0.2f,%0.2f)pi',...
                    obj.specfn_uneq_es.kxs_gxmg(ii) / pi, obj.specfn_uneq_es.kys_gxmg(ii) / pi) );
                elseif strcmp(mode, 'gymg')
                    ax2.XLim = [ min(obj.specfn_uneq_es.es_gymg(ii, :)), max(obj.specfn_uneq_es.es_gymg(ii, :)) ];
                                    title( sprintf('(kx,ky) = (%0.2f,%0.2f)pi',...
                    obj.specfn_uneq_es.kxs_gymg(ii) / pi, obj.specfn_uneq_es.kys_gymg(ii) / pi) );
                end
                ax2.XDir = 'reverse';
                ax2.XAxisLocation = 'top';

                ax_handles = cat(1, ax_handles, ax);
                ax2_handles = cat(1, ax2_handles, ax2);
            end

            % add momentum labels as annotation
            for ii = 1 : nplots
                ax = ax_handles(ii);
                height = 0.05;
                sep = 0.0;
                dim = [ax.Position(1), ax.Position(2) + ax.Position(4) - height - sep,...
                       ax.Position(3), height];
                  
                if strcmp(mode, 'gxmg')
                    text = sprintf('k exp = (%0.2f,%0.2f)pi \n %s',...
                        obj.specfn_uneq_es.kxs_gxmg(ii) / pi, obj.specfn_uneq_es.kys_gxmg(ii) / pi, thry_titles{ii});
                elseif strcmp(mode, 'gymg')
                    text = sprintf('k exp = (%0.2f,%0.2f)pi \n %s',...
                        obj.specfn_uneq_es.kxs_gymg(ii) / pi, obj.specfn_uneq_es.kys_gymg(ii) / pi, thry_titles{ii});
                end
                
                annotation('textbox', dim, 'String', text, 'FitBoxToText', 'on', 'FontSize', 10);
            end
            
        end
            
        function [fig_handle] = show_full_spectral_fn(obj)
            % display atom transfer at each k-vector as a function of
            % frequency. Grid of plots with each plot in the coordinates of
            % the pixel it represents.

%             fig_handle = plot_3d_array(obj.nr, obj.nr_unc, obj.kxx_mean(1, :), obj.kyy_mean(:, 1), obj.energies_2d, '', 'b.');
%             thry_spec_fn = obj.theory_scale_trpavg * squeeze( obj.dqmc_aavg.af(obj.theory_temp_index, :, :, :) ) + obj.theory_offset_trpavg;     
%             fig_handle = plot_3d_array(thry_spec_fn, '', obj.dqmc_kxs, obj.dqmc_aavg.kys, obj.dqmc_aavg.es, fig_handle, 'r-');
            
            nrows = size(obj.specfn_uneq_scaled.af, 1);
            ncols = size(obj.specfn_uneq_scaled.af, 2);
            nplots = nrows * ncols;
            
            y_limit = [0, max(obj.specfn_uneq_scaled.af(:)) * 1.2];
            
             ax_handles = [];
            ax2_handles = [];
            
            fig_handle = figure;
            for ii = 1:nrows
                for jj = 1:ncols
                    % get plot number. Plots are given in 
                    plot_index = jj + (ii - 1) * (ncols);
                    ax = subplot(nrows, ncols, plot_index);
                    ax.XTick = '';
                    ax.YTick = '';
        
                    grid on;
                    xlim( [min(obj.frqs_offsets_tunits), max(obj.frqs_offsets_tunits)] );

                    % create second y-axis for energy
                    ax2 = axes('Position', ax.Position);

                    ax2.YLim = ax.YLim;
                    % plot chemical potential
                    errorbar(ax2, squeeze(obj.specfn_uneq_scaled.es(ii, jj, :)), squeeze(obj.specfn_uneq_scaled.af(ii, jj, :)),...
                             squeeze(obj.specfn_uneq_scaled.afunc(ii, jj, :)), 'b.');
                    hold on;
                    plot(ax2, [ min(obj.specfn_uneq_scaled.es(ii, jj, :)), max(obj.specfn_uneq_scaled.es(ii, jj, :)) ], [0, 0], 'k');
                    ax2.YLim = [min(obj.specfn_uneq_scaled.af(ii, jj, :)), 1.2 * max(obj.specfn_uneq_scaled.af(ii, jj, :))];

                    % if dqmc data is present, plot it too
                    thry_ttl = '';
                    % one idea: fit fractions of cloud at different densities
                    % and see what works best. Verify it is close to the
                    % observed density distribution.
                    if ~isempty(obj.dqmc_aavg) 

                        % find nearest k-index
                        [~, min_index_kx] = min( (abs(obj.specfn_uneq_scaled.kxs(jj)) - obj.dqmc_aavg.kxs).^2 );
                        [~, min_index_ky] = min( (abs(obj.specfn_uneq_scaled.kys(ii)) - obj.dqmc_aavg.kys).^2 );

                        kx_curr = obj.dqmc_aavg.kxs(min_index_kx);
                        ky_curr = obj.dqmc_aavg.kys(min_index_ky);

                        % plot trap average
%                         thry_spec_fn = obj.theory_scale_trpavg * squeeze( obj.dqmc_aavg.af(min_index_ky, min_index_kx, :) ) + obj.theory_offset_trpavg;
                        thry_spec_fn = squeeze( obj.dqmc_aavg.af(min_index_ky, min_index_kx, : ));
                        plot(ax2, obj.dqmc_aavg.es, thry_spec_fn, 'r' );
                        hold on;
                        ax2.YLim(1) = min( ax2.YLim(1), 0);
                        ax2.YLim(2) = max(ax2.YLim(2), 1.2 * max(thry_spec_fn(:)));

                         % title with k-vectors
                        thry_ttl = sprintf('k th = (%0.2f,%0.2f)pi', kx_curr/pi, ky_curr/pi);
                        thry_titles{plot_index} = thry_ttl;
                    end

                    if ~isempty(obj.mu_center)
                        plot(ax2, [obj.mu_center, obj.mu_center], [ax2.YLim(1), ax2.YLim(2)], ...
                                'LineWidth', 2);

                    end

                    % set up second axis correctly relative to first.
                    ax2.XTick = '';
                    ax2.YTick = '';
                    ax2.Color = 'none';
                    ax2.XLim = [ min(obj.specfn_uneq_scaled.es(ii, jj, :)), max(obj.specfn_uneq_scaled.es(ii, jj, :)) ];
                    if obj.init_hfstate_lower_energy
                        % if the initial hyperfine state has lower energy than
                        % the final hyperfine state, E = e_k - nu, so the
                        % energy axis will be opposite to the nu axis (i.e.
                        % increasing rf => decreasing energy)
                        ax2.XDir = 'reverse';
                    else
                        % if the initial hyperfine state has higher energy than
                        % the final hyperfine state, E = e_k + nu, so the
                        % energy axis will be the same direction as the nu axis
                        ax2.XDir = 'normal';
                    end
                    ax2.XAxisLocation = 'top';
                    ax2.YAxisLocation = 'right';    
                    ax2.YTickLabel = '';


                    ax_handles = cat(1, ax_handles, ax);
                    ax2_handles = cat(1, ax2_handles, ax2);
                                        % plot chemical potential
                    if ~isempty(obj.mu_center)
                        plot(ax2, [obj.mu_center, obj.mu_center], [ax2.YLim(1), ax2.YLim(2)], ...
                                'LineWidth', 2);
                        hold on;
                    end
                end
            end
            
                    % add momentum labels as annotation
%             for ii = 1 : nplots
%                 ax = ax_handles(ii);
%                 height = 0.05;
%                 sep = 0.0;
%                 dim = [ax.Position(1), ax.Position(2) + ax.Position(4) - height - sep,...
%                        ax.Position(3), height];
%                 text = sprintf('k exp = (%0.2f,%0.2f)pi \n %s',...
%                     obj.specfn_uneq_es.kxs_gxmg(ii) / pi, obj.specfn_uneq_es.kys_gxmg(ii) / pi, thry_titles{ii});
%                 annotation('textbox', dim, 'String', text, 'FitBoxToText', 'on', 'FontSize', 10);
%             end
           
            for ii = 1 : nplots
                ax_handles(ii).YLim = ax2_handles(ii).YLim;
                ax_handles(ii).Position = ax2_handles(ii).Position;
            end    
        
            suptitle( sprintf('%s, full spectral fn', strrep(obj.identifier, '_', ' ')) );
        
        end
        
        function [fig_handle, axes] = show_edcs(obj, axes)
            % plot energy distribution curves (EDCs), which are A(k,w)*f(w)
            % at constant k values.
            
            if ~exist('axes', 'var') || isempty(axes)
                fig_handle = figure;
                axes = gca;
            else
                fig_handle = gcf;
            end
            
            % offset for each EDC
            offset_size = 0.2 * max( obj.specfn.af_gxmg(:) );
            offsets = offset_size * (1:size(obj.specfn.af_gxmg, 1));
            
            % for plotting peak line
            peak1_xs = obj.qp_dispersion_af_gxmg(:, 1);
            peak1_ys = zeros(size(peak1_xs));
            
            if obj.use_double_lor_dispersion
                peak2_xs = obj.qp_dispersion_af_gxmg(:, 2);
                peak2_ys = zeros(size(peak2_xs));
            end
            
            line_width = 2;
            % plot each curve
            for ii = 1:size(obj.specfn.af_gxmg, 1)

                plot(axes, obj.energy_bin_means, obj.specfn.af_gxmg(ii, :) + offsets(ii),...
                    'b-', 'LineWidth', line_width)
                hold on;
                plot(axes, [obj.energy_bin_means(1), obj.energy_bin_means(end)],...
                    [offsets(ii), offsets(ii)], 'k--', 'LineWidth', line_width/2);
                
                peak1_ys(ii) = interp1(obj.energy_bin_means,...
                           obj.specfn.af_gxmg(ii, :) + offsets(ii),...
                           peak1_xs(ii));
                       
                if obj.use_double_lor_dispersion
                    peak2_ys(ii) = interp1(obj.energy_bin_means,...
                               obj.specfn.af_gxmg(ii, :) + offsets(ii),...
                               peak2_xs(ii));
                end
            end
            grid on;
            
            % plot centers of EDC fits
            plot(axes, peak1_xs, peak1_ys, 'k', 'LineWidth', line_width);
            if obj.use_double_lor_dispersion
                plot(peak2_xs, peak2_ys, 'k', 'LineWidth', line_width);
            end
            
            % x-axis
            xlim( [min(obj.energy_bin_means), max(obj.energy_bin_means)] );
            xlabel('Energy (t)');
            
            % y-axis
            fn_klinearindex2y = @(kindex) interp1(obj.specfn_uneq_scaled.linear_index, offsets, kindex);
            axes.YTick = fn_klinearindex2y([0, 1, 2, 3]);
            axes.YTickLabel = {'\Gamma', 'X', 'M', '\Gamma'};
            ylabel('spectral intensity (arb)');

            title('Energy distribution curves');   
        end
        
        function [fig_handle, axes] = show_af_gxmg_expt(obj, axes, plot_log, use_thry_scaled)
            % Show an image of the spectral function along the GXMG path
            % throubh the Brillouin zone. The y-axis is energy and the
            % x-axis is quasimomentum.
            %
            % axes: an axis object on which to place the plot.
            %
            % plot_log: Boolean. If 1, plot on a log scale. Default is 0.
            
            if ~exist('axes', 'var') || isempty(axes)
                fig_handle = figure;
                axes = gca;
            else
                fig_handle = axes.Parent;
            end
            
            if ~exist('plot_log', 'var') || isempty(plot_log)
                plot_log = 0;
            end
            
            if ~exist('use_thry_scaled', 'var')
                use_thry_scaled = 0;
            end
            
            if use_thry_scaled
                spec_fn = obj.specfn_scaled.af_gxmg;
            else
                spec_fn = obj.specfn.af_gxmg;
            end
            
            if plot_log
                imagesc(axes, obj.specfn_uneq_es.linear_index, obj.energy_bin_means, log10(spec_fn'));
            else
                imagesc(axes, obj.specfn_uneq_es.linear_index, obj.energy_bin_means, spec_fn');
                ax.CLim(1) = 0;
            end
            axes.YDir = 'normal';
            colormap hot;
            hold(axes, 'on');

            % plot qusiparticle dispersion centers
            scatter(axes, obj.specfn_uneq_es.linear_index,...
                        obj.qp_dispersion_af_gxmg(:, 1), 40, 'b', 'x', 'LineWidth', 2);
            scatter(axes, obj.specfn_uneq_es.linear_index,...
                        obj.qp_dispersion_a_gxmg(:, 1), 40, 'g', 'x', 'LineWidth', 2);
            if obj.use_double_lor_dispersion
                scatter(axes, obj.specfn_uneq_es.linear_index,...
                        obj.qp_dispersion_af_gxmg(:, 2), 40, 'b', 'x', 'LineWidth', 2);
                scatter(axes, obj.specfn_uneq_es.linear_index,...
                        obj.qp_dispersion_a_gxmg(:, 2), 40, 'g', 'x', 'LineWidth', 2);
            end

            % plot non-interacting dispersion
            interp_indices = linspace(0, 3, 300);
            [kx_interp, ky_interp] = index2kvect(interp_indices);
             plot(axes, interp_indices, obj.latt_dispersion( kx_interp, ky_interp), 'b',...
                'LineWidth', 2);   
            
            % plot chemical potential
            if ~isempty(obj.mu_center)         
                plot(axes, [0, 3], [obj.mu_center, obj.mu_center],...
                    'b', 'LineWidth', 2);
            end
            
            % x -axis
            axes.XTick = [0, 1, 2, 3];
            axes.XTickLabel = {'\Gamma', 'X', 'M', '\Gamma'};
            xlabel('Momentum, K');

            % y-axis
            ytics = -12 : 4 : 12;
            axes.YTick = ytics;
            ylabel('Energy (t)');
            
            title(axes, 'A(k, w)*f(w)');
        end
             
        function [fig_handle, axes] = show_rfxfer_fn(obj, axes, plot_log)
            % show rf transfer along GXMG
            % A(k, ek - nu) * f(ek - nu)
            %
            % x-axis is quasimomentum and y-axis is rf frequency nu
            % subtracted from the non-interacting zeeman splitting in units
            % of the hopping t.
            
            if ~exist('axes', 'var') || isempty(axes)
                fig_handle = figure;
                axes = gca;
            else
                fig_handle = gcf;
            end
            
            if ~exist('plot_log', 'var') || isempty(plot_log)
                plot_log = 0;
            end
            
            if plot_log
                im = log10(obj.specfn_uneq_es.af_gxmg');
                im( ~isreal(im)) = nan;
                im = real(im);
                imagesc(axes, obj.specfn_uneq_es.linear_index, obj.frqs_offsets_tunits, im);
            else
                imagesc(axes, obj.specfn_uneq_es.linear_index, obj.frqs_offsets_tunits, obj.specfn_uneq_es.af_gxmg');
            end
            colormap hot;
            hold on;
            
            % plot fit centers for each k-vector
            scatter(axes, obj.specfn_uneq_es.linear_index,...
                          obj.lor_fit_params_gxmg(:, 1),...
                          40, 'b', 'x', 'LineWidth', 2);
            if obj.use_double_lor_dispersion
                scatter(axes, obj.specfn_uneq_es.linear_index,...
                              obj.lor_fit_params_gxmg(:, 5),...
                              40, 'b', 'x', 'LineWidth', 2);
            end
            plot([0, 3], [0, 0] , 'b', 'LineWidth', 2);
            
            % x-axis
            axes.XTick = [0, 1, 2, 3];
            axes.XTickLabel = {'\Gamma', 'X', 'M', '\Gamma'};
            xlabel('Momentum, K');
            
            % y-axis
            axes.YTick = obj.frqs_offsets_tunits(1:2:end);
            axes.YTickLabel = obj.frqs_hz_sorted(1:2:end) / 1e6;
            ylabel('rf frq (MHz)');
            title('transfer vs (nu, k) = A(k, e_k - nu) * f(e_k - nu)');     
        end
        
        function [fig_handle, axes] = show_specfn_ksummed(obj, axes)
            
            if ~exist('axes', 'var') || isempty(axes)
                fig_handle = figure;
                axes = gca;
            else
                fig_handle = gcf;
            end
            
            plot(axes, obj.energy_bin_means, obj.specfn_sumk, 'o--');
            hold on;
            interp_es = linspace(min(obj.energy_bin_means), max(obj.energy_bin_means), 300);
            plot(interp_es, lorentzian1D(obj.specfn_sumk_fitp ,interp_es));
            
            ylim([0 1.1 * max(obj.specfn_sumk)]);
            grid on;
            xlabel('Energy (t)');
            ylabel('sum_k A(k,w)f(w)');
            title( sprintf('sum_k A(k,w)f(w)\n HWHM = %0.1f(%.0f) t = %0.1f(%.0f) KHz',...
                obj.specfn_sumk_fitp(2), obj.specfn_sumk_fitp_err(2) * 10,...
                obj.specfn_sumk_fitp(2) * obj.t_over_h / 1e3,...
                obj.specfn_sumk_fitp_err(2) * obj.t_over_h / 1e3 * 10) );

        end
        
        function [fig_handle, axes] = show_rfxfer_ksummed(obj, axes)
            
            if ~exist('axes', 'var') || isempty(axes)
                fig_handle = figure;
                axes = gca;
            else
                fig_handle = gcf;
            end
            
            interp_frqs = linspace(min(obj.frqs_offsets_tunits), max(obj.frqs_offsets_tunits), 300);
            
            phs = [];
            strs = {};
            
            ph1 = plot(axes, obj.frqs_offsets_tunits, obj.dx * obj.dy * obj.rfxfer_sumk, 'bo--');
            hold on;
            plot(axes, interp_frqs, obj.dx * obj.dy * lorentzian1D(obj.rfxfer_sumk_fitp, interp_frqs), 'b');
            %
            phs = horzcat(phs, ph1);
            strs = horzcat(strs, 'momentum integrated'); 
            if ~isempty( obj.rfspec_int_df)
                int_frqs_tunits = (obj.rfxfer_density_resolved_frqs * 1e6 - obj.frq_nonint_transition_hz) / obj.t_over_h; 
                ph2 = plot(axes, int_frqs_tunits, obj.rfxfer_density_resolved, 'ro--');
                
                phs = horzcat(phs, ph2);
                strs = horzcat(strs, 'position resolved');
            end

            grid on;
            xlabel('Frequency (t)');
            ylabel('atom number');
            legend(phs, strs);
            axes.XLim = [ min(obj.frqs_offsets_tunits), max(obj.frqs_offsets_tunits)];
            
            title( sprintf('rf spectra \n HWHM = %0.1f(%.0f) t = %0.1f(%.0f) KHz',...
                   obj.rfxfer_sumk_fitp(2), obj.rfxfer_sumk_fitp_err(2) * 10,...
                   obj.rfxfer_sumk_fitp(2) * obj.t_over_h / 1e3,...
                   obj.rfxfer_sumk_fitp_err(2) * obj.t_over_h / 1e3 * 10) );
        end
        
        function [fig_handle, axes, plot_handles, legend_entries] = show_rfcenters(obj, axes)
            
            if ~exist('axes', 'var') || isempty(axes)
                fig_handle = figure;
                axes = gca;
            else
                fig_handle = gcf;
            end
            
            plot_handles = [];
            legend_entries = {};
            
            hold on;
            
            line_width = 2;
            marker_size = 20;
            
            % GXMG, first peak
            ph_gxmg = errorbar(axes, obj.specfn_uneq_es.linear_index, obj.lor_fit_params_gxmg(:, 1),...
                     obj.lor_fit_params_uncs_gxmg(:, 1), '.--',...
                     'LineWidth', line_width, 'MarkerSize', marker_size);
            plot_handles = horzcat(plot_handles, ph_gxmg);
            leg_entry = sprintf('gxmg, pk1, T = %0.2f(%.0f)', obj.temp_fit_struct.fit_params(2), obj.temp_fit_struct.std_errs(2) * 100);
            legend_entries = horzcat(legend_entries, leg_entry);     
            
            if obj.use_double_lor_dispersion
                ph_pk2_gxmg = errorbar(axes, obj.specfn_uneq_es.linear_index, obj.lor_fit_params_gxmg(:, 5),...
                     obj.lor_fit_params_uncs_gxmg(:, 5), '.--', 'Color', ph_gxmg.Color,...
                     'LineWidth', line_width, 'MarkerSize', marker_size);
                
                plot_handles = horzcat(plot_handles, ph_pk2_gxmg);
                leg_entry = sprintf('gxmg, pk2, T = %0.2f(%.0f)', obj.temp_fit_struct.fit_params(2), obj.temp_fit_struct.std_errs(2) * 100);
                legend_entries = horzcat(legend_entries, leg_entry); 
            end
                  
            ph_gymg = errorbar(axes, obj.specfn_uneq_es.linear_index, obj.lor_fit_params_gymg(:, 1),...
                     obj.lor_fit_params_uncs_gymg(:, 1), '.--',...
                     'LineWidth', line_width, 'MarkerSize', marker_size);
            plot_handles = horzcat(plot_handles, ph_gymg);
            leg_entry = sprintf('gymg, pk1, T = %0.2f(%.0f)', obj.temp_fit_struct.fit_params(2), obj.temp_fit_struct.std_errs(2) * 100);
            legend_entries = horzcat(legend_entries, leg_entry);    
                 
            if obj.use_double_lor_dispersion
                ph_pk2_gymg = errorbar(axes, obj.specfn_uneq_es.linear_index, obj.lor_fit_params_gymg(:, 5),...
                     obj.lor_fit_params_uncs_gymg(:, 5), '.--', 'Color', ph_gymg.Color,...
                     'LineWidth', line_width, 'MarkerSize', marker_size);
               
                 plot_handles = horzcat(plot_handles, ph_pk2_gymg);
                leg_entry = sprintf('gymg, pk2, T = %0.2f(%.0f)', obj.temp_fit_struct.fit_params(2), obj.temp_fit_struct.std_errs(2) * 100);
                legend_entries = horzcat(legend_entries, leg_entry);    
 
            end
                 
            grid on;
            ymin = min(obj.lor_fit_params_gxmg(:, 1)) - 0.2 * abs(min(obj.lor_fit_params_gxmg(:, 1)));
            
            if obj.use_double_lor_dispersion
                ymax = max(max(obj.lor_fit_params_gxmg(:, 1)) + 0.2 * abs(max(obj.lor_fit_params_gxmg(:, 1))),...
                        max(obj.lor_fit_params_gxmg(:, 5)) + 0.2 * abs(max(obj.lor_fit_params_gxmg(:, 5))) );
            else
                ymax = max(obj.lor_fit_params_gxmg(:, 1)) + 0.2 * abs(max(obj.lor_fit_params_gxmg(:, 1)));
            end
            
            ylim([ymin, ymax]);
              
%             axes.YLim(1) = 0;
            axes.XTick = [0, 1, 2, 3];
            axes.XTickLabel = {'\Gamma', 'X', 'M', '\Gamma'};
            xlabel('Momentum, k');
            ylabel('frq offset(t)');
            legend([ph_gxmg, ph_gymg], {'gxmg', 'gymg'});          
            
            title('peak position vs. k');
        end
        
        function [fig_handle, axes, plot_handles, legend_entries] = show_qpdispersions(obj, axes)
            % display expt and theory dispersions for a and af
            
            if ~exist('axes', 'var') || isempty(axes)
                fig_handle = figure;
                axes = gca;
            else
                fig_handle = gcf;
            end
            
            plot_handles = [];
            legend_entries = {};  
            
            line_width = 2;
            marker_size = 20;
            
            % non-interacting dispersion
            plot(axes, obj.specfn_uneq_es.linear_index, obj.latt_dispersion(obj.specfn_uneq_es.kxs_gxmg, obj.specfn_uneq_es.kys_gxmg),...
                'k', 'LineWidth', line_width);
            hold on;
            
            % plot mu
            plot(axes, [0, 3], [obj.mu_center, obj.mu_center], 'color', 'blue', 'LineWidth', line_width);
%             legend_entries = horzcat(legend_entries, 'non-int', 'mu');
            
            % plot qp dispersion
            ph = errorbar(axes, obj.specfn_uneq_es.linear_index, obj.qp_dispersion_af_gxmg(:, 1),...
                     obj.qp_dispersion_af_gxmg_unc(:, 1), '.--',...
                     'LineWidth', line_width, 'MarkerSize', marker_size);
            
            plot_handles = horzcat(plot_handles, ph);
            leg_entry = sprintf('gxmg, T = %0.2f(%.0f)', obj.temp_fit_struct.fit_params(2), obj.temp_fit_struct.std_errs(2) * 100);
            legend_entries = horzcat(legend_entries, leg_entry);
                 
            hold on;
            
            if obj.use_double_lor_dispersion
                errorbar(axes, obj.specfn_uneq_es.linear_index, obj.qp_dispersion_af_gxmg(:, 2),...
                     obj.qp_dispersion_af_gxmg_unc(:, 2), '.--', 'Color', ph.Color,...
                     'LineWidth', line_width, 'MarkerSize', marker_size);
            end
            
            % gymg
            ph2 = errorbar(axes, obj.specfn_uneq_es.linear_index, obj.qp_dispersion_af_gymg(:, 1),...
                     obj.qp_dispersion_af_gymg_unc(:, 1), '.--',...
                     'LineWidth', line_width, 'MarkerSize', marker_size);
            plot_handles = horzcat(plot_handles, ph2);
            leg_entry = sprintf('gymg, T = %0.2f(%.0f)',...
                        obj.temp_fit_struct.fit_params(2), obj.temp_fit_struct.std_errs(2) * 100);
            legend_entries = horzcat(legend_entries, leg_entry);
                
            
            if obj.use_double_lor_dispersion
                errorbar(axes, obj.specfn_uneq_es.linear_index, obj.qp_dispersion_af_gymg(:, 2),...
                     obj.qp_dispersion_af_gymg_unc(:, 2), '.--', 'Color', ph2.Color,...
                     'LineWidth', line_width, 'MarkerSize', marker_size);
            end
               
            % expt, after dividing fermi function
            ph = errorbar(axes, obj.specfn_uneq_scaled.linear_index, obj.qp_dispersion_a_gxmg(:, 1),...
                          obj.qp_dispersion_a_gxmg_unc(:, 1), '.--',...
                          'LineWidth', line_width, 'MarkerSize', marker_size);
            plot_handles = horzcat(plot_handles, ph);
            leg_entry = sprintf('a, gxmg expt');
            legend_entries = horzcat(legend_entries, leg_entry);
            
            
            % theory
            if ~isempty(obj.dqmc_a_dispersion_gxmg_fitp)
                for jj = 1 : obj.max_num_peaks
                        ph_a = errorbar(axes, obj.dqmc_aavg.linear_index, obj.dqmc_af_dispersion_gxmg_fitp(:, 1 + (jj-1) * 4),...
                                   obj.dqmc_af_dispersion_gxmg_stderr(:, 2 + (jj-1) * 4),...
                                   'bx', 'LineWidth', 2);
                        ph_b = errorbar(axes, obj.dqmc_aavg.linear_index, obj.dqmc_a_dispersion_gxmg_fitp(:, 1 + (jj-1) * 4),...
                                   obj.dqmc_af_dispersion_gxmg_stderr(:, 2 + (jj-1) * 4),...
                                   'gx', 'LineWidth', 2);
                end
                
                plot_handles = horzcat(plot_handles, ph_a, ph_b);
                legend_entries = horzcat(legend_entries, 'thry af', 'thry a');    
            end
            
            ylim([-5, 8]);
            ax = gca;
            ax.XTick = [0, 1, 2, 3];
            ax.XTickLabel = {'\Gamma', 'X', 'M', '\Gamma'};
            
            xlabel('Momentum, k');
            ylabel('Energy (t)');
            
            legend(plot_handles, legend_entries);
            title ('dispersions');
        end
          
        function [fig_handle, axes] = show_rfhwhm(obj, axes)
            
            if ~exist('axes', 'var') || isempty(axes)
                fig_handle = figure;
                axes = gca;
            else
                fig_handle = gcf;
            end
            
            line_width = 2;
            marker_size = 20;
            
            ph = errorbar(axes, obj.specfn_uneq_es.linear_index, obj.lor_fit_params_gxmg(:, 2),...
                     obj.lor_fit_params_uncs_gxmg(:, 2), '.--',...
                     'LineWidth', line_width, 'MarkerSize', marker_size);
            if obj.use_double_lor_dispersion
                errorbar(axes, obj.specfn_uneq_es.linear_index, obj.lor_fit_params_gxmg(:, 6),...
                     obj.lor_fit_params_uncs_gxmg(:, 6), '.--', 'Color', ph.Color,...
                     'LineWidth', line_width, 'MarkerSize', marker_size);
            end
            
            hold on;
            
            ph2 = errorbar(axes, obj.specfn_uneq_es.linear_index, obj.lor_fit_params_gymg(:, 2),...
                     obj.lor_fit_params_uncs_gymg(:, 2), '.--',...
                     'LineWidth', line_width, 'MarkerSize', marker_size);
            if obj.use_double_lor_dispersion
                errorbar(axes, obj.specfn_uneq_es.linear_index, obj.lor_fit_params_gymg(:, 6),...
                 obj.lor_fit_params_uncs_gymg(:, 6), '.--', 'Color', ph2.Color,...
                 'LineWidth', line_width, 'MarkerSize', marker_size);
            end
                 
            grid on;
            legend({'gxmg', 'gymg'});
            axes.XTick = [0, 1, 2, 3];
            axes.XTickLabel = {'\Gamma', 'X', 'M', '\Gamma'};
            axes.YLim = [0, 4];
            legend({'gxmg', 'gymg'});
            xlabel('Momentum, k');
            ylabel('HWHM (t)');                      
            
            title('HWHM');
        end
        
        % non-interacting rf
        function process_rf_nonint(obj)
            % analyze non-interacting rf spectroscopy data. Should already
            % be stored in instance. You must call load_rf_nonint before
            % you call this function.
            
            [frqs_unique, dens_unique, dens_unique_unc,...
                num_points, std_dev, fig_handle] = obj.rfspec_df.showAllVsIndVar;
            
            % store some results in instance
            obj.rfspec_frqs = frqs_unique;
            obj.rfspec_xfer = dens_unique;
            obj.rfspec_xfer_unc = dens_unique_unc;
            
            nbins = obj.rfspec_df.NBins;
            min_frq = min(frqs_unique);
            max_frq = max(frqs_unique);
            min_width = 0.05e-3;
            max_width = max_frq - min_frq;
            frq_guess = mean(frqs_unique);

            % parameters for single lorentzian fit
            init_params_lor = [frq_guess, 0.002, 0.1, 0];
            fixed_params_lor = [0, 0, 0, 0];
            lbs_lor = [min_frq, min_width, 0, 0];
            ubs_lor = [max_frq, max_width, inf, inf];
            
            fit_params_lor = zeros(nbins, length(init_params_lor));
            std_errs_lor = zeros(nbins, length(init_params_lor));
            
            % parameters for double lorentzian fit
            init_params_dbllor = [ [frq_guess, 0.002, 0.1, 0], [frq_guess, 0.002, 0.0, 0] ];
            fixed_params_dbllor = [ [0, 0, 0, 0], [0, 0, 0, 1] ];
            lbs_dbllor = [ [min_frq, min_width, 0, 0],...
                           [min_frq, min_width, 0, -1] ];
            ubs_dbllor = [ [max_frq, max_width, 0.2, 0.02],...
                           [max_frq, max_width, 0.2, 1] ];

            fit_params_dbllor = zeros(nbins, length(init_params_dbllor));
            std_errs_dbllor = zeros(nbins, length(init_params_dbllor));
    
            for ii = 1 : nbins
                fit_unc = dens_unique_unc(:, ii);
                fit_unc(fit_unc == 0) = mean(fit_unc);
                
                % single lorentzian fit
                [fit_params_lor(ii, :), ~, ffh_lor, std_errs_lor(ii, :), ~] = ...
                    fit1D(frqs_unique, dens_unique(:, ii), 1 ./ fit_unc.^2, ...
                    {'lorentzian1D'}, init_params_lor, fixed_params_lor,...
                    lbs_lor, ubs_lor);
                
                % double lorentzian fit
                [fit_params_dbllor(ii, :), ~, ffh_dbllor, std_errs_dbllor(ii, :), chi_sqr] = fit1D(...
                    frqs_unique, dens_unique(:, ii), 1 ./ fit_unc.^2,...
                    {'lorentzian1D', 'lorentzian1D'}, init_params_dbllor,...
                    fixed_params_dbllor, lbs_dbllor, ubs_dbllor);
                
                
                % sort double lorentzian fit parameters by lower peak
                if fit_params_dbllor(ii, 5) < fit_params_dbllor(ii, 1)
                    fitp_helper = zeros(size(fit_params_dbllor(ii, :)));
                    fitp_helper(1:3) = fit_params_dbllor(ii, 5:7);
                    fitp_helper(5:7) = fit_params_dbllor(ii, 1:3);
                    fit_params_dbllor(ii, :) = fitp_helper;

                    std_err_helper = zeros(size(std_errs_dbllor(ii, :)));
                    std_err_helper(1:3) = std_errs_dbllor(ii, 5:7);
                    std_err_helper(5:7) = std_errs_dbllor(ii, 1:3);
                    std_errs_dbllor(ii, :) = std_err_helper;
                end
            end

            % store fit results in instance
            obj.rfspec_lor_fitp = fit_params_lor;
            obj.rfspec_lor_fitp_unc = std_errs_lor; 
            obj.rfspec_dbl_lor_fitp = fit_params_dbllor;
            obj.rfspec_dbl_lor_fitp_unc = std_errs_dbllor;
           
        end
        
        function [fig_handle] = plot_rf_nonint(obj, fig_size_pixels)
            
            if ~exist('fig_size_pixels', 'var')
                a = get(0, 'Screensize');
                fig_size_pixels = a(3:4);
            end
            
            fig_handle = figure;
            fig_handle.Position = horzcat([0, 0], fig_size_pixels);
            
            if isempty(obj.rfspec_frqs) || isempty(obj.rfspec_xfer) || ...
               isempty(obj.rfspec_xfer_unc) || isempty(obj.rfspec_lor_fitp) || ...
               isempty(obj.rfspec_lor_fitp_unc) || isempty(obj.rfspec_dbl_lor_fitp) || ...
               isempty(obj.rfspec_lor_fitp_unc) || isempty(obj.rfspec_df)
               warning('Prematurely returning from plot_rf_nonint because necessary data was not present');
               return;
            end
            
              % interpolation frequencies
            interp_frqs = linspace(min(obj.rfspec_frqs), max(obj.rfspec_frqs), 300);
            
            nbins = obj.rfspec_df.NBins;
           % figure
            nrows = floor(sqrt(nbins));
            ncols = ceil(nbins / nrows);
            
            for ii = 1 : nbins
                subplot(nrows, ncols, ii);
                errorbar(obj.rfspec_frqs, obj.rfspec_xfer(:, ii),...
                                      obj.rfspec_xfer_unc(:, ii), '-o');
                hold on;
                plot(interp_frqs, lorentzian1D(obj.rfspec_lor_fitp, interp_frqs));
                plot(interp_frqs, lorentzian1D(obj.rfspec_dbl_lor_fitp(1:4), interp_frqs) + ...
                                  lorentzian1D(obj.rfspec_dbl_lor_fitp(5:8), interp_frqs) );
                xlabel('Frequency (MHz)');
                ylabel('Filling');
                xlim([min(obj.rfspec_frqs), max(obj.rfspec_frqs)])
                ax = gca;
                ax.YLim(1) = 0;
                ax.YLim(2) = 1.2 * max(obj.rfspec_xfer(:));
                grid on;               
                
                ttl = sprintf('Frq = %0.6f(%0.0f) MHz, HWHM = %0.3f(%0.0f) KHz\n %.0f - %0.f sites', ...
                    obj.rfspec_lor_fitp(ii, 1),...
                    round(obj.rfspec_lor_fitp_unc(ii, 1), 6) * 1e6,...
                    obj.rfspec_lor_fitp(ii, 2) * 1e3,...
                    round(obj.rfspec_lor_fitp_unc(ii, 2) * 1e3, 3) * 1e3,...
                    obj.rfspec_df.BinEdges(ii), obj.rfspec_df.BinEdges(ii + 1) );
                
                title(ttl, 'FontSize', 8);
            end
            
            suptitle(sprintf('%s', obj.rfspec_df.DatasetString));
            
        end
        
        % interacting rf
        function process_rf_int(obj)
            
            if ~isempty(obj.rfspec_int_df.Occs_ImgAvg)
                [frqs_unique, atom_nums, atom_nums_unc, num_points, std_dev, rep_num]...
                 = obj.rfspec_int_df.get_versus_ind_var(obj.rfspec_int_df.AtomNumbers);

                
                [~, ns, ns_unc, ~, ~, ~]...
                 = obj.rfspec_int_df.get_versus_ind_var(obj.rfspec_int_df.Occs_AzAvgStack');
                
                % fit lorentzian to determine peak transfer
                [~, max_index] = max(atom_nums);
                span = max(frqs_unique) - min(frqs_unique);

                init_params = [ frqs_unique(max_index), span/10, max(atom_nums), 0];
                fixed_params = [0, 0, 0, 0];
                lbs = [min(frqs_unique), min( frqs_unique(2:end) - frqs_unique(1:end-1) ), 0, -inf];
                ubs = [max(frqs_unique), span, 2 * max(atom_nums), inf];

%                 [fit_params, ~, ffh, std_err, chi_sqr] = fit1D(frqs_unique, atom_nums,...
%                     [], {'lorentzian1D'}, init_params, fixed_params, lbs, ubs);
                [fit_params, ~, ffh, std_err, chi_sqr] = fit1D(frqs_unique, ns,...
                    [], {'lorentzian1D'}, init_params, fixed_params, lbs, ubs);

                obj.rfxfer_density_resolved_density = ns;
                obj.rfxfer_density_resolved_density_unc = ns_unc;
                obj.rfxfer_density_resolved_lor_fitp = fit_params;
                obj.rfxfer_density_resolved = atom_nums;
                obj.rfxfer_density_resolved_frqs = frqs_unique;
            end
        end
        
        % comparison functions
        function [fig_handle] = compare_sets(instance_stack)
            
            fig_handle = figure;
            nrows = 2;
            ncols = 4;
            
            leg = cell(length(instance_stack), 1);
            for ii = 1:length(instance_stack)
                inst = instance_stack(ii);
                leg{ii} = strrep(inst.identifier, '_', ' ');
                
                % qp dispersion gxmg
                p1 = subplot(nrows, ncols, 1);
                errorbar(inst.linear_index, inst.qp_dispersion_af_gxmg(:, 1),...
                         inst.qp_dispersion_af_gxmg_unc, '--o');
                grid on;
                hold on;
                
                % qp dispersion gymg
                p2 = subplot(nrows, ncols, 2);
                errorbar(inst.linear_index, inst.qp_dispersion_af_gymg(:, 1),...
                         inst.qp_dispersion_af_gymg_unc, '--o');
                grid on;
                hold on;
                
                p3 = subplot(nrows, ncols, 3);
                errorbar(inst.linear_index, inst.lor_fit_params_gxmg(:, 2),...
                         inst.lor_fit_params_uncs_gxmg(:, 2), '--o');
                grid on;
                hold on;
                
                p4 = subplot(nrows, ncols, 4);
                errorbar(inst.linear_index, inst.lor_fit_params_gymg(:, 2),...
                         inst.lor_fit_params_uncs_gymg(:, 2), '--o');
                grid on;
                hold on;
                
                p5 = subplot(nrows, ncols, 5);
                errorbar(inst.linear_index, inst.lor_fit_params_gxmg(:, 1),...
                         inst.lor_fit_params_uncs_gxmg(:, 1), '--o');
                grid on;
                hold on;
                
                p6 = subplot(nrows, ncols, 6);
                errorbar(inst.linear_index, inst.lor_fit_params_gymg(:, 1),...
                         inst.lor_fit_params_uncs_gymg(:, 1), '--o');
                grid on;
                hold on;
                
                p7 = subplot(nrows, ncols, 7);
                plot(inst.energy_bin_means, inst.specfn_sumk, '--o');
                grid on;
                hold on;
                
                p8 = subplot(nrows, ncols, 8);
                plot(inst.frqs_offsets_tunits, inst.rfxfer_sumk, '--o');
                grid on;
                hold on;
           
                
            end
            legend(leg);

            p1.Title.String = 'qp dispersion GXMG';
            p1.YLabel.String = 'Energy (t)';
            p1.XTick = [0, 1, 2, 3];
            p1.XTickLabel = {'\Gamma', 'X', 'M', '\Gamma'};
            p1.YLim = [-0.5, 8.5];
            
            p2.Title.String = 'qp dispersion GYMG';
            p2.YLabel.String = 'Energy (t)';
            p2.XTick = [0, 1, 2, 3];
            p2.XTickLabel = {'\Gamma', 'Y', 'M', '\Gamma'};
            p2.YLim = [-0.5, 8.5];
            
            p3.Title.String = 'dispersion HWHM GXMG';
            p3.YLabel.String = 'HWHM (KHz)';
            p3.XTick = [0, 1, 2, 3];
            p3.XTickLabel = {'\Gamma', 'Y', 'M', '\Gamma'};
            p3.YLim = [-0.2, 5];
            
            p4.Title.String = 'dispersion HWHM GXMG';
            p4.YLabel.String = 'HWHM (KHz)';
            p4.XTick = [0, 1, 2, 3];
            p4.XTickLabel = {'\Gamma', 'Y', 'M', '\Gamma'};
            p4.YLim = [-0.2, 5];
            
            p5.Title.String = 'energy shift GXMG';
            p5.YLabel.String = 'energy shift(KHz)';
            p5.XTick = [0, 1, 2, 3];
            p5.XTickLabel = {'\Gamma', 'Y', 'M', '\Gamma'};
            
            p6.Title.String = 'energy shift GYMG';
            p6.YLabel.String = 'energy shift(KHz)';
            p6.XTick = [0, 1, 2, 3];
            p6.XTickLabel = {'\Gamma', 'Y', 'M', '\Gamma'};
            
            p7.Title.String = '\Sigma_k A(k,\omega) f(\omega)';
            p7.YLabel.String = '';
            
            p8.Title.String = '\Sigma_K A(k, \epsilon_k - \nu) f(\epsilon_k - \nu)';
            p8.YLabel.String = '';

            suptitle('comparison, arpes sets');
        end
        
        function [fig_handle] = compare_specfns(instance_stack, fig_size_pix, log_scale)
            
            if ~exist('fig_size_pix', 'var') || isempty(fig_size_pix)
                fig_size_pix = [800, 800];
            end
            
            if ~exist('log_scale', 'var') || isempty(log_scale)
                log_scale = 1;
            end
            
            fig_pos = horzcat([100, 100], fig_size_pix);

            % plot EDCs
            fig_handle = figure();
            fig_handle.Position = fig_pos;

            nrows = ceil(sqrt(length(instance_stack)));
            ncols = ceil(length(instance_stack) / nrows);
            
            for ii = 1:length(instance_stack)
                inst = instance_stack(ii);

                axes_h = subplot(nrows, ncols, ii);
                inst.show_af_gxmg_expt(axes_h, log_scale);
                axes_h.FontSize = 16;
                id = strrep(inst.identifier, '_', ' ');
                title(axes_h, sprintf('T/t = %0.2f(%.0f), U/t = %0.2f(%.f)',...
                    inst.temp_fit_struct.fit_params(2), inst.temp_fit_struct.std_errs(2) * 100,...
                    inst.temp_fit_struct.fit_params(1), inst.temp_fit_struct.std_errs(1) * 100));
            end
            suptitle('spectral functions');
            
        end
        
        function [fig_handle] = compare_edcs(instance_stack, fig_size_pix)
            
            if ~exist('fig_size_pix', 'var') || isempty(fig_size_pix)
                fig_size_pix = [800, 800];
            end
            fig_pos = horzcat([100, 100], fig_size_pix);

            % plot EDCs
            fig_handle = figure();
            fig_handle.Position = fig_pos;

            nrows = ceil(sqrt(length(instance_stack)));
            ncols = ceil(length(instance_stack) / nrows);
            
            for ii = 1:length(instance_stack)   
                inst = instance_stack(ii);
                
                axes_h = subplot(nrows, ncols, ii);
                inst.show_edcs(axes_h);
                axes_h.FontSize = 16;
                id = strrep(inst.identifier, '_', ' ');
                title(axes_h, sprintf('T/t = %0.2f(%.0f), U/t = %0.2f(%.f)',...
                    inst.temp_fit_struct.fit_params(2), inst.temp_fit_struct.std_errs(2) * 100,...
                    inst.temp_fit_struct.fit_params(1), inst.temp_fit_struct.std_errs(1) * 100));
            end
            suptitle('Energy distribution curves (EDCs)');
            
        end
        
        function [fig_handle] = compare_rfxfer_fns(instance_stack, fig_size_pix)
            
            if ~exist('fig_size_pix', 'var') || isempty(fig_size_pix)
                fig_size_pix = [800, 800];
            end
            fig_pos = horzcat([100, 100], fig_size_pix);

            % plot EDCs
            fig_handle = figure();
            fig_handle.Position = fig_pos;

            nrows = ceil(sqrt(length(instance_stack)));
            ncols = ceil(length(instance_stack) / nrows);
            
            for ii = 1:length(instance_stack)
                inst = instance_stack(ii);
                
                axes_h = subplot(nrows, ncols, ii);
                inst.show_rfxfer_fn(axes_h);
                axes_h.FontSize = 16;
                id = strrep(inst.identifier, '_', ' ');
                title(axes_h, sprintf('T/t = %0.2f(%.0f), U/t = %0.2f(%.f)',...
                    inst.temp_fit_struct.fit_params(2), inst.temp_fit_struct.std_errs(2) * 100,...
                    inst.temp_fit_struct.fit_params(1), inst.temp_fit_struct.std_errs(1) * 100));
            end
            suptitle('RF Transfers Fns');
            
        end
        
        function [fig_handle] = compare_energyshift_k(instance_stack, fig_size_pix)
            
            if ~exist('fig_size_pix', 'var') || isempty(fig_size_pix)
                fig_size_pix = [800, 800];
            end
            fig_pos = horzcat([100, 100], fig_size_pix);

            % plot EDCs
            fig_handle = figure();
            fig_handle.Position = fig_pos;
            axes_h = gca;
            
            plot_handles = [];
            leg_entries = {};
            
            for ii = 1:length(instance_stack)
                inst = instance_stack(ii);
                [~, ~, ph, leg_e] = inst.show_rfcenters(axes_h);
                hold on;     
                plot_handles = horzcat(plot_handles, ph(1));
                leg_entries = horzcat(leg_entries, leg_e{1});
            end
            title('dispersion');
            
            % TODO: fix legend
            legend(plot_handles, leg_entries);
        end
        
        function [fig_handle] = compare_dispersions(instance_stack, fig_size_pix, plot_relative_to_mu)
            
            if ~exist('fig_size_pix', 'var') || isempty(fig_size_pix)
                fig_size_pix = [800, 800];
            end
            fig_pos = horzcat([100, 100], fig_size_pix);

            fig_handle = figure();
            fig_handle.Position = fig_pos;
            axes_h = gca;
            
            plot_handles = [];
            leg_entries = {};
            
            for ii = 1:length(instance_stack)
                inst = instance_stack(ii);
                [~, ~, ph, leg_e] = inst.show_qpdispersions(axes_h);
                hold on;     
                plot_handles = horzcat(plot_handles, ph(1));
                leg_entries = horzcat(leg_entries, leg_e{1});
            end
            
            title('dispersion');
            legend(plot_handles, leg_entries);
            
        end
        
        % exporting
        function export_expt2txt(obj, save_dir)
            if ~exist('save_dir', 'var')
                save_dir = '';
            end
            
            % parameters
            fname = sprintf('%s_dqmc_fit_params.txt', obj.identifier);
            fname = fullfile(save_dir, fname);
            names_cell = {'mu (t)', 'U (t)', 'U unc', 'T (t)', 'T unc',...
                          'h (t)', 'h unc', 'efficiency', 'efficiency unc',...
                          'doubles efficiency', 'doubles efficiency unc'};
            dqmc_params = vertcat(obj.temp_fit_struct.fit_params, obj.temp_fit_struct.std_errs);
            data = horzcat(obj.mu_center, dqmc_params(:)');
            save_data_file(data, names_cell, '\t', '', '', fname);
            
            % export density
            fname = sprintf('%s_expt_density.txt', obj.identifier);
            fname = fullfile(save_dir, fname);
            export_gnuplot_mat(fname, obj.ones.Occs_ImgAvg + obj.threes.Occs_ImgAvg,...
                               1 : size(obj.dqmc_density, 2), 1 : size(obj.dqmc_density, 1) );
            
            % export rf transfer density images
            for ii = 1 : size(obj.specfn_uneq_es.af, 3)
                fname = sprintf('%s_nrsymm_frq=%0.3fMHz=%0.1ftoverh.txt',...
                        obj.identifier, obj.frqs_hz_sorted(ii) / 1e6, obj.frqs_offsets_tunits(ii));
                fname = fullfile(save_dir, fname);
                dlmwrite(fname, obj.specfn_uneq_es.af(:, :, ii),'\t');
                
                fname = sprintf('%s_nrsymm_frq=%0.3fMHz=%0.1ftoverh_splot.txt',...
                        obj.identifier, obj.frqs_hz_sorted(ii) / 1e6, obj.frqs_offsets_tunits(ii));
                fname = fullfile(save_dir, fname);
                [xx, yy] = meshgrid( 1 : size(obj.specfn_uneq_es.af, 2), 1 : size(obj.specfn_uneq_es.af, 2));
                export_gnuplot_splot_mat(fname, obj.specfn_uneq_es.af(:, :, ii), xx, yy);
            end
            
            
            if isempty( obj.specfn_scaled)
                return;
            end
            
            % A*f GXMG
            fname = sprintf('%s_af_gxmg.txt', obj.identifier);
            fname = fullfile(save_dir, fname);
            export_gnuplot_mat(fname, obj.specfn_scaled.af_gxmg', obj.specfn_scaled.linear_index, obj.specfn_scaled.es);

            
           % A*f for each k-vector
            for ii = 1 : length(obj.specfn_uneq_scaled.kys)
                for jj = 1 : length(obj.specfn_uneq_scaled.kxs)  
                    kx = round(obj.specfn_uneq_scaled.kxs(jj), 12);
                    ky = round(obj.specfn_uneq_scaled.kys(ii), 12);
                    
                    % sometimes rounds to zero but sprintf will give
                    % '-0.00'. This seems to fix it.
                    if kx == 0
                        kx = 0;
                    end
                    
                    if ky == 0
                        ky = 0;
                    end
                    
                    if kx >= 0 && ky >= 0
                        fname = sprintf('%s_af_kx=%0.2f_ky=%0.2f_vs_frq.txt', obj.identifier, kx, ky);
                        fname = fullfile(save_dir, fname);
                        names_cell = {'frqs (MHz)', 'frqs - nonint (KHz)', 'frqs - nonint (t/h)', 'energy (t)',...
                                      'Af(kx, ky) (1/a^2)', 'Af unc (1/a^2)' };
                                  
                        data = horzcat(obj.frqs_hz_sorted'/1e6,...
                                       obj.frqs_offsets_tunits' * obj.t_over_h / 1e3,...
                                       obj.frqs_offsets_tunits',...
                                       squeeze(obj.specfn_uneq_scaled.es(ii, jj, :)),...
                                       squeeze(obj.specfn_uneq_scaled.af(ii, jj, :)),...
                                       squeeze(obj.specfn_uneq_scaled.afunc(ii, jj, :)) );

                                   
                        save_data_file(data, names_cell, '\t', '', '', fname);
                    end
                end
            end
            
            % export af dispersion
            fname = sprintf('%s_af_dispersion.txt', obj.identifier);
            fname = fullfile(save_dir, fname);
            names_cell = {'linear index', 'center (t)', 'center unc', 'hwhm (t)', 'hwhm unc',...
                                          'center (t)', 'center unc', 'hwhm (t)', 'hwhm unc'};
            data = [obj.specfn.linear_index', obj.qp_dispersion_af_gxmg(:, 1), obj.qp_dispersion_af_gxmg_unc(:, 1),...
                                              obj.lor_fit_params_gxmg(:, 2), obj.lor_fit_params_uncs_gxmg(:, 2),...
                                              obj.qp_dispersion_af_gxmg(:, 2), obj.qp_dispersion_af_gxmg_unc(:, 2),...
                                              obj.lor_fit_params_gxmg(:, 6), obj.lor_fit_params_uncs_gxmg(:, 6)];
            save_data_file(data, names_cell, '\t', '', '', fname);
            
            % export a dispersion
            fname = sprintf('%s_a_dispersion.txt', obj.identifier);
            fname = fullfile(save_dir, fname);
            names_cell = {'linear index', 'center (t)', 'center unc', 'hwhm (t)', 'hwhm unc',...
                                          'center (t)', 'center unc', 'hwhm (t)', 'hwhm unc'};
            data = [obj.specfn.linear_index', obj.qp_dispersion_a_gxmg(:, 1), obj.qp_dispersion_a_gxmg_unc(:, 1),...
                                              obj.lor_fit_params_gxmg(:, 2), obj.lor_fit_params_uncs_gxmg(:, 2),...
                                              obj.qp_dispersion_a_gxmg(:, 2), obj.qp_dispersion_a_gxmg_unc(:, 2),...
                                              obj.lor_fit_params_gxmg(:, 6), obj.lor_fit_params_uncs_gxmg(:, 6)];
            save_data_file(data, names_cell, '\t', '', '', fname);
             
            % export temperature fit experiment data
            fname = sprintf('%s_temp_fit_exp.txt', obj.identifier);
            fname = fullfile(save_dir, fname);
            names_cell = {'n', 'nerr', 'nup', 'nup err', 'upup', 'upup err',...
                          'ndn', 'ndn err', 'dndn', 'dndn err', 'ns', 'ns err',...
                          'nsns', 'nsns err', 'd', 'd err', 'dd', 'dd err'};
            data = [obj.temp_fit_struct.n', obj.temp_fit_struct.n_err',...
                    obj.temp_fit_struct.nup', obj.temp_fit_struct.nup_err',...
                    obj.temp_fit_struct.nup_corr', obj.temp_fit_struct.nup_corr_err',...
                    obj.temp_fit_struct.ndn', obj.temp_fit_struct.ndn_err',...
                    obj.temp_fit_struct.ndn_corr', obj.temp_fit_struct.ndn_corr_err',...
                    obj.temp_fit_struct.ns', obj.temp_fit_struct.ns_err',...
                    obj.temp_fit_struct.ns_corr', obj.temp_fit_struct.ns_corr_err',...
                    obj.temp_fit_struct.nd', obj.temp_fit_struct.nd_err',...
                    obj.temp_fit_struct.nd_corr', obj.temp_fit_struct.nd_corr_err'];
            
            save_data_file(data, names_cell, '\t', '', '', fname);
            
            % temperature theory fit
            fname = sprintf('%s_temp_fit_thry_U=%0.2f_T=%0.2f_Fid=%0.2f_DFid=%0.2f.txt', obj.identifier,...
                            obj.temp_fit_struct.fit_params(1), obj.temp_fit_struct.fit_params(2),...
                            obj.temp_fit_struct.fit_params(4), obj.temp_fit_struct.fit_params(5));
            fname = fullfile(save_dir, fname);
            names_cell = {'n', 'nup', 'upup', 'ndn', 'dndn', 'ns', 'nsns', 'd', 'dd'};
            data = [obj.temp_fit_struct.n_interp', obj.temp_fit_struct.nup_fit',...
                    obj.temp_fit_struct.nup_corr_fit', obj.temp_fit_struct.ndn_fit',...
                    obj.temp_fit_struct.ndn_corr_fit', obj.temp_fit_struct.ns_fit',...
                    obj.temp_fit_struct.ns_corr_fit', obj.temp_fit_struct.nd_fit',...
                    obj.temp_fit_struct.nd_corr_fit'];
            save_data_file(data, names_cell, '\t', '', '', fname);
            
            if ~isempty(obj.rfspec_int_df)
                % export rf transfer
                fname = sprintf('%s_interacting_rf_density_resolved.txt', obj.identifier);
                fname = fullfile(save_dir, fname);
                names_cell = {'frqs (t)', '<n>', '<n> unc'};

                int_frqs_tunits = (obj.rfxfer_density_resolved_frqs * 1e6 - obj.frq_nonint_transition_hz) / obj.t_over_h; 
                data = horzcat(int_frqs_tunits', obj.rfxfer_density_resolved_density, obj.rfxfer_density_resolved_density_unc);
                save_data_file(data, names_cell, '\t', '', '', fname);
                
                fname = sprintf('%s_interacting_rf_density_resolved_fit.txt', obj.identifier);
                fname = fullfile(save_dir, fname);
                names_cell = {'frqs (t)', '<n>'};
                
                fit_params = obj.rfxfer_density_resolved_lor_fitp;
                fit_params(1) = (fit_params(1) * 1e6 - obj.frq_nonint_transition_hz) / obj.t_over_h; 
                fit_params(2) = fit_params(2) * 1e6 / obj.t_over_h;
                
                frq_interp = linspace(min(int_frqs_tunits), max(int_frqs_tunits), 300);
                fit = lorentzian1D(fit_params, frq_interp);
                data = horzcat(frq_interp', fit');
                save_data_file(data, names_cell, '\t', '', '', fname);
                
            end
            
        end
            
        function export_thry2txt(obj, save_dir)
                if ~exist('save_dir', 'var')
                    save_dir = '';
                end
            
                if obj.use_expt_broadening
                    broadening = obj.broadening_khz;
                else
                    broadening = 0;
                end
                dqmc_desc = sprintf('broad=%0.1fKHz_U=%0.2f_T=%0.2f',...
                                    broadening, obj.dqmc_aavg.us, obj.dqmc_aavg.ts);
                
                % export images comparable to 2D functions, i.e. these are
                % A(k, ek \pm nu) * f(ek \pm nu)
                for ii = 1 : size(obj.dqmc_aavg_scaled.af, 3)
                    fname = sprintf('%s_dqmc_nr_frq=%0.1ftoverh_%s.txt',...
                            obj.identifier, obj.frqs_offsets_tunits(ii), dqmc_desc);
                    fname = fullfile(save_dir, fname);
                    dlmwrite(fname, obj.dqmc_aavg_scaled.af(:, :, ii),'\t');
                end
                
                % export inferred density
                fname = sprintf('%s_dqmc_density_%s.txt', obj.identifier, dqmc_desc);
                fname = fullfile(save_dir, fname);
                export_gnuplot_mat(fname, obj.dqmc_density, 1 : size(obj.dqmc_density, 2), 1 : size(obj.dqmc_density, 1) );
                
                % GXMG a*f function
                fname = sprintf('%s_dqmc_af_trpavg_gxmg_%s.txt',...
                                obj.identifier, dqmc_desc);
                fname = fullfile(save_dir, fname);
                export_gnuplot_mat(fname, obj.dqmc_aavg.af_gxmg', obj.dqmc_aavg.linear_index, obj.dqmc_aavg.es);
                
                % GXMG a function
                fname = sprintf('%s_dqmc_a_trpavg_gxmg_%s.txt',...
                                obj.identifier, dqmc_desc);
                fname = fullfile(save_dir, fname);
                export_gnuplot_mat(fname, obj.dqmc_aavg.a_gxmg', obj.dqmc_aavg.linear_index, obj.dqmc_aavg.es);
                
                % a and af for each k-vector
                for ii = 1 : length(obj.dqmc_aavg.kys)
                    for jj = 1 : length(obj.dqmc_aavg.kxs)  
                        kx = round(obj.dqmc_aavg.kxs(jj), 12);
                        ky = round(obj.dqmc_aavg.kys(ii), 12);
                        
                    
                        % sometimes rounds to zero but sprintf will give
                        % '-0.00'. This seems to fix it.
                        if kx == 0
                            kx = 0;
                        end

                        if ky == 0
                            ky = 0;
                        end
                        
                        if kx < 0 || ky < 0
                            continue;
                        end
                        
                        fname = sprintf('%s_dqmc_a_af_trapavg_%s_kx=%0.2f_ky=%0.2f.txt',...
                                obj.identifier, dqmc_desc, kx, ky);
                        fname = fullfile(save_dir, fname);
                        names_cell = {'Energy (t)', 'nu-nu_0 (t/h)', 'a*f', 'a'};
                        
                        if obj.init_hfstate_lower_energy
                            nu_offset = obj.latt_dispersion(kx, ky) - obj.dqmc_aavg.es';
                        else
                            nu_offset = obj.dqmc_aavg.es' - obj.latt_dispersion(kx, ky);
                        end
                        
                        data = horzcat(obj.dqmc_aavg.es', nu_offset,...
                            squeeze(obj.dqmc_aavg.af(ii, jj, :)),...
                            squeeze(obj.dqmc_aavg.a(ii, jj, :)));
                        save_data_file(data, names_cell, '\t', '', '', fname);
                    end
                end
                
                
                 % af full
                fname = sprintf('%s_dqmc_af_trapavg_%s.txt', obj.identifier, dqmc_desc);
                fname = fullfile(save_dir, fname);
                names_cell = {'Energy (t)'};
                data = obj.dqmc_aavg.es';
                for ii = 1 : size(obj.dqmc_aavg.af, 1)
                    for jj = 1 : size(obj.dqmc_aavg.af, 2)
                        if obj.dqmc_aavg.kxs(jj) >=0 && obj.dqmc_aavg.kys(jj) >= 0
                            data = horzcat(data, squeeze(obj.dqmc_aavg.af(ii, jj, :)) );
                            names_cell = horzcat(names_cell, sprintf('kx/pi=%0.2f ky/pi=%0.2f', obj.dqmc_aavg.kxs(jj), obj.dqmc_aavg.kys(ii)));
                        end
                    end
                end
                save_data_file(data, names_cell, '\t', '', '', fname);
                
                % aa dispersion
                fname = sprintf('%s_dqmc_a_gxmg_dispersion_%s.txt', obj.identifier, dqmc_desc);
                fname = fullfile(save_dir, fname);
                names_cell = {'linear index'};
                data = [obj.dqmc_aavg.linear_index'];
                for ii = 1 : obj.max_num_peaks
%                     names_cell = horzcat(names_cell, {'center (t)', 'width (t)'});
%                     data = horzcat(data, obj.dqmc_a_dispersion_gxmg_fitp(:, 1 + (ii - 1) * 4),...
%                                          obj.dqmc_a_dispersion_gxmg_fitp(:, 2 + (ii - 1) * 4));
                    names_cell = horzcat(names_cell, {'center (t)', 'width (t)', 'height', 'prominence'});
                    data = horzcat(data, obj.dqmc_a_dispersion_gxmg_fitp(:, 1 + (ii - 1) * 4 : 4 + (ii - 1) * 4));
                end
                save_data_file(data, names_cell, '\t', '', '', fname);
                
                % af dispersion
                fname = sprintf('%s_dqmc_af_gxmg_dispersion_%s.txt', obj.identifier, dqmc_desc);
                fname = fullfile(save_dir, fname);
                names_cell = {'linear index'};
                data = [obj.dqmc_aavg.linear_index'];
                for ii = 1 : obj.max_num_peaks
%                     names_cell = horzcat(names_cell, {'center (t)', 'width (t)'});
%                     data = horzcat(data, obj.dqmc_af_dispersion_gxmg_fitp(:, 1 + (ii - 1) * 4),...
%                                          obj.dqmc_af_dispersion_gxmg_fitp(:, 2 + (ii - 1) * 4));
                    names_cell = horzcat(names_cell, {'center (t)', 'width (t)', 'height', 'prominence'});
                    data = horzcat(data, obj.dqmc_af_dispersion_gxmg_fitp(:, 1 + (ii - 1) * 4 : 4 + (ii - 1) * 4));
                end
                save_data_file(data, names_cell, '\t', '', '', fname);
             
        end
        
        % saving and loading
        function saveStruct(obj, save_dir, fname, save_all_folders)
            % saveStruct(obj, save_dir, fname, save_all_folders)
            %
            % save as struct. Better than saving as class, because the
            % struct can be parsed even if you don't have DiffuseSet.m, or
            % if the saved file comes from a different version of
            % DiffuseSet.m
            %
            % save_dir is the directory to save structure in
            %
            % fname is file name of saved structure
            %
            % save_all_folders is bool. Specifies whether or not to save
            % DataFolder sets also. DataFolder sets are not saved in this
            % structure, but in their own separate structures to reduce
            % overhead for working with DiffuseSet data. Typically they are
            % not needed after initial processing.
            %
            
            % check arguments
            if ~exist('save_dir','var') || isempty(save_dir)
                save_dir = pwd;
            end
            
            if ~exist('fname','var') || isempty(fname)
                fname = fullfile([sprintf('%s_struct', obj.identifier), '.mat']);
            end
            
            if ~exist('save_all_folders', 'var') || isempty(save_all_folders)
                save_all_folders = 1;
            end
            
            fpath = fullfile(save_dir, fname);
            
            AsStruct = struct();
            fields = fieldnames(obj);
            for ii = 1:length(fields)
                % member of this class that aren't DataFolder classes are
                % stored together. DataFolder objects are stored on their
                % own, one structure per object.
                if isa(obj.(fields{ii}), 'DataFolder')
                    if save_all_folders
                        struct_array = [];
                        for jj = 1:length(obj.(fields{ii}))
                            struct_current = obj.(fields{ii})(jj).getStruct();
                            struct_array = vertcat(struct_array, struct_current);
%                             obj.(fields{ii})(jj).saveStruct(save_dir);
                        end
                        AsStruct.(fields{ii}) = struct_array;
                    end
                elseif strcmp(fields{ii}, 'dqmc')
                    % don't save this field
                else
                     AsStruct.(fields{ii}) = obj.(fields{ii});
                end
            end
            
            fprintf('Saved to %s\n',fpath);
            save(fpath, '-struct', 'AsStruct');
        end
        
        function saveAllFolders(obj, save_dir)
            % save data from all DataFolder objects as structs in save_dir
            if ~exist('save_dir','var') || isempty(save_dir)
                save_dir = pwd;
            end
            
            fields = fieldnames(obj);
            for ii = 1:length(fields)
                % member of this class that aren't DataFolder classes are
                % stored together. DataFolder objects are stored on their
                % own, one structure per object.
                if isa(obj.(fields{ii}), 'DataFolder')
                    for jj = 1:length(obj.(fields{ii}))
                        obj.(fields{ii})(jj).saveStruct(save_dir);
                    end
                end
            end
        end
        
        function loadStruct(obj, StructToLoad, load_folders)
            % loadStruct(obj, StructToLoad, load_folders)
            %
            % load class from structure.
            % 
            % StructToLoad: may either be the actual structure, or it may be
            % the file path to the structure
            %
            % load_folders: Boolean, if 1 load DataFolder() fields of
            % class. If not, do not load these.
            
            % if StructToLoad is a string, load the structure
            if ischar(StructToLoad)
                StructToLoad = load(StructToLoad);
            end
            
            % now StructToLoad should be our desired structure
            if ~isstruct(StructToLoad)
                error('StructToLoad was not a structure');
            end
            
            if ~exist('load_folders', 'var')
                load_folders = 1;
            end
            
            % list of DataFolder() fields, which will need to be handled
            % differently than other fields
            df_fields = {'data_folder_stack', 'rfspec_df', 'rfspec_int_df', 'ones', 'threes', 'singles', 'doubles'};
            fn = @(s1, s2) strcmp(s1, s2);
            
            % clear object prior to loading to avoid accidentally keeping
            % any data.
            ObjFields = fieldnames(obj);
            for ii = 1:length(ObjFields)
                obj.(ObjFields{ii}) = [];
            end
            
            % load struct and get field names
            StructFields = fieldnames(StructToLoad);
            % loop over struct field names and assign any shared to our
            % class
            for ii = 1:length(StructFields)
                cell_fn = @(c) fn(StructFields{ii}, c);
                if ~any( cell_fn(df_fields) )
%                 if ~any(strcmp(df_fields, StructFields{ii}))
                   
                    try   
                        obj.(StructFields{ii}) = StructToLoad.(StructFields{ii});
                    catch
                        fprintf('field %s was present in loaded struct, but it is not a field of arpes_data class. Skipped this field. \n',StructFields{ii});
                    end
                end
            end
            
            if ~load_folders        
                return;
            end
            
            for ii = 1 : length(df_fields)
                % get field in question
                fld = StructToLoad.(df_fields{ii});

                % check if it is empty. If not, load it as a DataFolder
                % instance and store it in our class
                if isempty(fld)  
                    continue;
                end
                
                if length(fld) == 1
                    % fld is a single DataFolder
                    df = DataFolder();
                    df.loadStruct(fld);
                    obj.(df_fields{ii}) = df;
                else
                    % fld is a stack of DataFolder's
                    df_stack = [];
                    for jj = 1 : length(fld)
                        df_temp = DataFolder();
                        df_temp.loadStruct(fld(jj));
                        % add to stack
                        df_stack = vertcat(df_stack, df_temp);
                    end
                    obj.(df_fields{ii}) = df_stack; 
                end
            end
            
        end
                
    end
    
end

