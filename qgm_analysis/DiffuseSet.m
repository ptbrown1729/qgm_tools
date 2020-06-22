classdef DiffuseSet < handle
    % class for organizing DataFolder objects and analyzing them for
    % diffusion/transport experiment.
    % Do full analysis of dataset using DiffuseSet and variable defining
    % script files. Assume all sub-directories of
    % this main directory house all the data for a given temperature.
    %
    % DiffuseSetSummary will enter each directory. It will then run DiffuseSet
    % for each k-vector. All the excluded pictures etc for a given k-vector are
    % defined in a matlab script file of the form set*.m. As the analysis is
    % done, most of the results will be saved in a subdirectory set*.
    %
    % Before performing a new analysis, run "clean_old_analysis" to move the
    % older analysis to the 'hide' subfolders.
    
    properties
        % Define datasets to use
        DescriptStr
        DateCell
        dfs = []
        dfs_folder_numbers = []
        dfs_index_folders = []
        dfs_index_normfolders = []
        dfs_index_ones = []
        dfs_index_threes = []
        dfs_index_singles = []
        
        Ones = []
        Threes = []
        Singles = []
        Folders = []
        NormFolders = []
        Times = []
        ExcludePics = {};
        ExcludePics_Norm = {};
        
        % analysis settings
        FixPhase = 1
        Angle = 90.5%+90
        Saving = 1
        save_path
        
        % Import params for full set
        Period
        PeriodUnc
        Temp
        TempUnc  
        % density information
        DensAvg
        DensAvgUnc
        DensGrads
        DensGradsUnc
        % also store filling/occupation info for ALL folders
        AllDensAvg
        AllDensAvgUnc
        AllDensGrads
        AllDensGradsUnc
        
        % radial data
        Density % real density azimuthal average, from blowing shots
        DensityUnc
        occ_vstime % actually the density for the time folders
        occ_vstime_unc
        
        % Leaning towards shouldn't be storing these quantities here...
        BinAvgDist % average distance of points in bin. Should be the same as RadPos
        BinAvgDistUnc
        
        % density profiles summed along one direction
        AllPos
        AllCuts
        AllUnc
        AllCuts_NoBgPos
        AllCuts_NoBg
        AllCuts_NoBgUnc
        AllCutsCorr
        AllCutsCorrUnc
        % also store cuts at right angles to see how flat things are
        AllOrthPos
        AllOrthCuts
        AllOrthCutsUnc
        AllOrthCutsCorr
        AllOrthCutsCorrUnc
        
        % fits to 1D sinusoids/modulation profiles
        fp_modulation
        stderr_modulation
        fp_modulation_nobgsub
        stderr_modulation_nobgsub
        fp_mod_names = {'Frequency','Amplitude','Phase','Background'};
        
        % fits to amplitude of modulation profiles versus times
        Fpdo % Gamma and D
        SEdo
        Fpdo_bs
        SEdo_bs
        Fpdo_offset % Gamma and D, allowing for offset
        SEdo_offset
        offset
        offset_unc
       
        DQMCPath = '\\128.112.86.75\lithium\Publications\05_Charge-Diffusion\DQMC\DQMC_n=8_T=0.25to14.97_U=8to8_mu=-107.92to0.00_GridParams_InterpolatingFn.mat' 
    end
    
    methods
        function obj = DiffuseSet()
            %instantiate class.
        end
        
        function analyzeDataFolders(obj, DateCell, Folders, NormFolders, Times,...
                ExcludePics, ExcludePics_Norm, ...
                Cx, Cy, OnesFolder, ThreesFolder, ExcludeOnes,...
                ExcludeThrees, BinWidth, ComputeCorrelator,...
                CutoffRad, SinglesFolder, ExcludeSingles,...
                save_dir, BinEndPts)
            %analyze all datafolders for a single k-vector for a single
            %Temperature     
            
            if ~exist('NormFolders','var') || isempty(NormFolders)
                NormFolders = Folders(end)*ones(1,length(Folders));
            end
            
            if length(NormFolders) == 1
                NormFolders = NormFolders*ones(1,length(Folders));
            end
            
            if length(NormFolders) ~= length(Folders)
                error('NormFolders was a different length than Folders');
            end
            
            if length(ExcludePics_Norm) == 1
                ExcludePics_Norm = repmat(ExcludePics_Norm,[1,length(NormFolders)]);
            end
            
            if length(ExcludePics_Norm) ~= length(NormFolders)
                ExcludePics_Norm = cell([1,length(NormFolders)]);
            end
            
            if length(ExcludePics) < length(Folders)
                %ensure ExcludePics are the right length
                ExPics = cell([1,length(Folders)]);
                for ii = 1:length(ExcludePics)
                    ExPics{ii} = ExcludePics{ii};
                end
                ExcludePics = ExPics;
            end
            
            if ~exist('ExcludePics','var') || isempty(ExcludePics)
                ExcludePics = cell(length(Folders),1);
            end
            
            if length(ExcludePics) ~= length(Folders)
                ExcludePics = cell([1,length(Folders)]);
            end
            
            if ~exist('Cx','var')
                Cx = [];
            end
            
            if ~exist('Cy','var')
                Cy = [];
            end
            
            if ~exist('OnesFolder','var')
                OnesFolder = [];
            end
            
            if ~exist('ThreesFolder','var')
                ThreesFolder = [];
            end
            
            if ~exist('ExcludeOnes','var')
                ExcludeOnes = [];
            end
            if ~exist('ExcludeThrees','var')
                ExcludeThrees = [];
            end
            
            if ~exist('SinglesFolder','var')
                SinglesFolder = [];
            end
            
            if ~exist('ExcludeSingles','var')
                ExcludeSingles = [];
            end
            
            if ~exist('ComputeCorrelator','var') || isempty(ComputeCorrelator)
                ComputeCorrelator = 0;
            end
            
            if ~exist('BinWidth', 'var') || isempty(BinWidth)
                BinWidth = 1;
            end
            
            if ~exist('CutoffRad', 'var') || isempty(CutoffRad)
                CutoffRad = 14;
            end       
            
            if ~exist('save_dir','var') || isempty(save_dir)
                save_dir = pwd;
            end
            
            if ~exist('BinEndPts', 'var') || isempty(BinEndPts)
                BinEndPts = sqrt(linspace(0, 40^2, 40));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % assign arguments to class members
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.save_path = save_dir;
            obj.DateCell = DateCell;
            obj.Folders = Folders;
            obj.Times = Times;
            obj.ExcludePics = ExcludePics;
            
            % what is the point of this?
            obj.dfs = [];
            obj.AllPos = [];
            obj.AllCuts = [];
            obj.AllUnc = [];
            obj.AllCutsCorr = [];
            obj.AllCutsCorrUnc = [];
            
            obj.AllOrthPos = [];
            obj.AllOrthCuts = [];
            obj.AllOrthCutsUnc = [];
            obj.AllOrthCutsCorr = [];
            obj.AllOrthCutsCorrUnc = [];
            
            obj.AllDensAvg = [];
            obj.AllDensAvgUnc = [];
            obj.AllDensGrads = [];
            obj.AllDensGradsUnc = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % combine all folders we want to analyze, keeping only unique
            % entries to avoid replicating analysis
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [AllFolders, ia, ic] = unique([Folders, NormFolders,...
                OnesFolder, ThreesFolder, SinglesFolder], 'stable');
            AllExcludePics = {ExcludePics{:}, ExcludePics_Norm{:}, ...
                ExcludeOnes, ExcludeThrees, ExcludeSingles};
            AllExcludePicsUnique = cell([1,length(AllFolders)]);
            
            % get all exclude pics for a given folder, no matter which
            % exclusion field is used to define them. e.g. if you have a
            % "ones" folder which is also a normalization folder, and you
            % defined different pictures to be excluded in "ExcludesOnes"
            % and "ExcludePics_Norm", then we want to exclude all of those
            % pictures.
            for ii = 1:length(AllFolders)
                combined_excludes = AllExcludePics(ic == ii);
                AllExcludePicsUnique{ii} = horzcat(combined_excludes{:});
            end
            
            % get indices of different kinds of folders. The idea is that,
            % if you want to access only that specific kind of folder, you
            % should index obj.dfs(obj.dfs_index_kind-of-folder)
            obj.dfs_folder_numbers = AllFolders;
            func = @(n) find(obj.dfs_folder_numbers == n);
            obj.dfs_index_folders = arrayfun(func, Folders);
            obj.dfs_index_normfolders = arrayfun(func, NormFolders);
            if ~isempty(OnesFolder)
                obj.dfs_index_ones = find(obj.dfs_folder_numbers == OnesFolder);
            else
                obj.dfs_index_ones = [];
            end
            
            if ~isempty(ThreesFolder)
                obj.dfs_index_threes = find(obj.dfs_folder_numbers == ThreesFolder);
            else
                obj.dfs_index_threes = [];
            end
            
            if ~isempty(SinglesFolder)
                obj.dfs_index_singles = find(obj.dfs_folder_numbers == SinglesFolder);
            else
                obj.dfs_index_singles = [];
            end
            
            obj.DescriptStr = sprintf('%04d-%02d-%02d_Folders=%03d-%03d',...
                obj.DateCell{1}, obj.DateCell{2}, obj.DateCell{3},...
                min(obj.dfs_folder_numbers), max(obj.dfs_folder_numbers));
            
            %check if there is already a saved version of this dataset...if
            %so, load it and test each folder to see if it needs to be
            %reanalyzed.
            fname = fullfile(save_dir, ['*_struct', '.mat']);
            files = dir(fname);
            %             if exist(fname,'file')
            if length(files) == 1
                fprintf('Found and loaded previous analysis file.\n');
                loaded_set = load(fullfile(files.folder, files.name));
                loaded_set = loaded_set.AsStruct;
            end
            
            % Find the folder we want to process in previously saved
            % data. If nothing has changed in our analysis settings, we
            % don't need to analyze the data again from scratch.
            for ii = 1:length(AllFolders)
                df = DataFolder();
                DSet = obj.DateCell;
                DSet{4} = AllFolders(ii);
                DSet{5} = 1; DSet{6} = 1;
                
                % assign properties to df
                df.ExcludedPictures = AllExcludePicsUnique{ii};
                if ii == 1
                    df.AzAvgType = 'spatial';
                    DistGrid = [];
                else
                    df.AzAvgType = 'external';
                    df.CenterStyle = 'Fixed';
                    df.CroppedPicStartCoords = obj.dfs(1).CroppedPicStartCoords;
                    df.Cx_AzAvg = obj.dfs(1).Cx_AzAvg;
                    df.Cy_AzAvg = obj.dfs(1).Cy_AzAvg;
                    df.AzAvgCentering = 'external';
                end
                
                df_saved = DataFolder();
                if exist('loaded_set', 'var')
                    try 
                        df_saved.loadStruct(loaded_set.dfs(AllFolders(ii) ==...
                            loaded_set.dfs_folder_numbers));
                    catch err
                        disp(err.message);
                        df_saved.CenterStyle = 'not a real center style...';
                    end
                else
                    %so won't accidentally be equal to the other set
                    df_saved.CenterStyle = 'not a real center style...';
                end
                
                % if sets are the same, we don't need to do the initial
                % processing again.
                if df.compare_dsets(df_saved)
                    df = df_saved;
                    fprintf('Sets equal, using saved set for %s\n', df.DatasetString);
                else
                    df.initialize(DSet, df.ExcludedPictures,...
                        BinEndPts, df.AzAvgType, DistGrid);
                end
                DistGrid = df.DistGrid;
                
                % main function for analyzing the modulation. This deals
                % with defining the cropping region, computing dens and
                % corr and averaging along the direction orthogonal to the
                % modulation. Can select an arbitrary angle.
                [fig_handle, XPos, CutSum, Unc, CutsCorr, CutsCorrUnc, LineFitPX, LineFitPXErr,...
                    YPos, OrthCut, OrthCutUnc, OrthCutsCorr, OrthCutsCorrUnc,...
                    LineFitPOrth, LineFitPOrthErr,...
                    CropRegDens, CropRegDensUnc, CropRegDensSD] ...
                    = df.sumDirection(obj.Angle * pi / 180 , Cx, Cy, CutoffRad,...
                    BinWidth, ComputeCorrelator);
                
                % save diagnostic picture for info about azimuthal
                % averaging, reconstruction errors, atom number stability,
                % lattice geometry, etc.
                if obj.Saving
                    fname = sprintf('%s_sum_direction.fig', df.DatasetString);
                    fpath = fullfile(obj.save_path, fname);
                    savefig(fig_handle, fpath);
                    
                    df.showAzAvg(obj.Saving, obj.save_path);
                end
                % save flattening analysis for folders that should not have
                % a modulation.
                if any( ii == [obj.dfs_index_normfolders, obj.dfs_index_ones, obj.dfs_index_threes] )
                    [flat_fig, ~, ~] = df.showFlatness([], [], obj.Saving, obj.save_path, [], 0);
                end
                
                obj.AllPos = cat(2, obj.AllPos, XPos);
                obj.AllCuts = cat(2, obj.AllCuts, CutSum);
                obj.AllUnc = cat(2, obj.AllUnc, Unc);
                obj.AllCutsCorr = cat(2, obj.AllCutsCorr, CutsCorr);
                obj.AllCutsCorrUnc = cat(2, obj.AllCutsCorrUnc, CutsCorrUnc);
                
                obj.AllOrthPos = cat(2, obj.AllOrthPos, YPos);
                obj.AllOrthCuts = cat(2, obj.AllOrthCuts, OrthCut);
                obj.AllOrthCutsUnc = cat(2, obj.AllOrthCutsUnc, OrthCutUnc);
                obj.AllOrthCutsCorr = cat(2, obj.AllOrthCutsCorr, OrthCutsCorr);
                obj.AllOrthCutsCorrUnc = cat(2, obj.AllOrthCutsCorrUnc, OrthCutsCorrUnc);
                
                obj.AllDensAvg = cat(2, obj.AllDensAvg, CropRegDens);
                obj.AllDensAvgUnc = cat(2, obj.AllDensAvgUnc, CropRegDensUnc);
                obj.AllDensGrads = cat(2, obj.AllDensGrads, [LineFitPX(2); LineFitPOrth(2)]);
                obj.AllDensGradsUnc = cat(2, obj.AllDensGradsUnc, [LineFitPXErr(2); LineFitPOrthErr(2)]);
                
                obj.dfs = cat(1, obj.dfs,df);
            end
            
            % extract summary data from all folders
            obj.occ_vstime = {obj.dfs(obj.dfs_index_folders).Occs_AzAvg};
            obj.occ_vstime_unc = {obj.dfs(obj.dfs_index_folders).Occs_AzAvgUnc};
            obj.BinAvgDist = {obj.dfs(obj.dfs_index_folders).BinAvg};
            obj.BinAvgDistUnc = {obj.dfs(obj.dfs_index_folders).BinUnc};
            
            % keeping these for backwards compatibility and ease of access.
            obj.Ones = obj.dfs(obj.dfs_index_ones);
            obj.Threes = obj.dfs(obj.dfs_index_threes);
            obj.Singles = obj.dfs(obj.dfs_index_singles);
            % get full density
            if ~isempty(obj.Ones) && ~isempty(obj.Threes)
                obj.Density = obj.Ones.Occs_AzAvg + obj.Threes.Occs_AzAvg;
                obj.DensityUnc = sqrt(obj.Ones.Occs_AzAvgUnc.^2 + obj.Threes.Occs_AzAvgUnc.^2);
                
                ones_avg_dens = obj.AllDensAvg(obj.dfs_index_ones);
                ones_avg_dens_unc = obj.AllDensAvgUnc(obj.dfs_index_ones);
                threes_avg_dens = obj.AllDensAvg(obj.dfs_index_threes);
                threes_avg_dens_unc = obj.AllDensAvgUnc(obj.dfs_index_threes);
                ones_grads = obj.AllDensGrads(:,obj.dfs_index_ones);
                ones_grads_unc = obj.AllDensGradsUnc(:,obj.dfs_index_ones);
                threes_grads = obj.AllDensGrads(:,obj.dfs_index_threes);
                threes_grads_unc = obj.AllDensGradsUnc(:,obj.dfs_index_threes);
                
                obj.DensAvg = ones_avg_dens + threes_avg_dens;
                obj.DensAvgUnc = sqrt(ones_avg_dens_unc.^2 + ...
                    threes_avg_dens_unc.^2);
                
                obj.DensGrads = ...
                    (ones_avg_dens*ones_grads + threes_avg_dens*threes_grads) / ...
                    (ones_avg_dens + threes_avg_dens);
                %TODO write correct expression for uncertainty
                obj.DensGradsUnc = sqrt(ones_grads_unc.^2 + ...
                    threes_grads_unc.^2);
            end
            
            if obj.Saving
                for ii = 1:length(obj.dfs)
                    obj.dfs(ii).saveStruct(save_dir);
                end
            end
        end
        
        function analyzeSingleWavelength(obj, PeriodGuess, PeriodFixed,...
                n_low_cutoff, use_quantity_to_fit, temp_fit_init_p, save_dir,...
                do_bootstrap, use_offsets)
            %analyzeSingleWavelength(obj, PeriodGuess, FixedPeriod,...
            %  n_low_cutoff, use_quantity_to_fit, temp_fit_init_p, save_dir)
            %
            % fit results of analyzeDataFolders to diffusion models, etc.
            % for a single wavevector
            

            if ~exist('do_bootstrap', 'var') || isempty(do_bootstrap)
                do_bootstrap = 1;
            end
            
            if ~exist('use_offsets', 'var') || isempty(use_offsets)
                use_offsets = 0;
            end
            
            
            Times = obj.Times;
            XPositions = obj.AllPos(:, obj.dfs_index_folders);
            FitCuts = obj.AllCuts(:, obj.dfs_index_folders);
            FitCutsUnc = obj.AllUnc(:, obj.dfs_index_folders);
            NormCuts = obj.AllCuts(:, obj.dfs_index_normfolders);
            NormCutsUnc = obj.AllUnc(:, obj.dfs_index_normfolders);
            
            export_fname_stem = sprintf('%s_1d_averages', obj.DescriptStr);
            [fig_handle, obj.fp_modulation, obj.stderr_modulation,...
                obj.AllCuts_NoBgPos, obj.AllCuts_NoBg, obj.AllCuts_NoBgUnc]...
                = obj.fitCuts(PeriodGuess, PeriodFixed,...
                Times, XPositions,...
                FitCuts, FitCutsUnc,...
                NormCuts, NormCutsUnc,...
                save_dir, obj.Saving, export_fname_stem);                     
            
            if obj.Saving
                % figure
                fig_fname = sprintf('%s_FitCuts.fig', obj.DescriptStr);
                fpath = fullfile(save_dir, fig_fname);
                savefig(fig_handle, fpath);
                
                % png
                [f, n, ~] = fileparts(fpath);
                png_fname = fullfile(f, [n, '.png']);
                print(png_fname, '-dpng');
            end
            
            obj.Period = 1 / obj.fp_modulation(1, 1);
            obj.PeriodUnc = obj.stderr_modulation(1, 1) / obj.fp_modulation(1, 1)^2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % fit one dimesional averages, no background subtraction
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            export_fname_stem = sprintf('%s_1d_averages_nosub', obj.DescriptStr);
            [fig_handle, obj.fp_modulation_nobgsub, obj.stderr_modulation_nobgsub,...
                ~, ~, ~]...
                = obj.fitCuts(PeriodGuess, PeriodFixed,...
                Times, XPositions,...
                FitCuts, FitCutsUnc,...
                [], [],...
                save_dir, obj.Saving, export_fname_stem);
            
                        
           if obj.Saving
               % figure
                fig_fname = sprintf('%s_FitCuts_nobg.fig', obj.DescriptStr);
                fpath = fullfile(save_dir, fig_fname);
                savefig(fig_handle, fpath);
                
                % png
                [f, n, ~] = fileparts(fpath);
                png_fname = fullfile(f, [n, '.png']);
                print(png_fname, '-dpng');
           end
            
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % Compare modulation decay with and without background
           % subtraction
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           fig_handle = obj.compareDiffuseAnalysis();
           
           if obj.Saving
                % figure
                fig_fname = sprintf('%s_1davgs_bgsubtraction_vs_not.fig',...
                    obj.DescriptStr);
                fpath = fullfile(save_dir, fig_fname);
                savefig(fig_handle, fpath);
                
                % png
                [f, n, ~] = fileparts(fpath);
                png_fname = fullfile(f, [n, '.png']);
                print(png_fname, '-dpng');
           end
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % fit amplitude decays of 1D modulation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(obj.dfs_index_folders) > 3
                
                % new style
                period = obj.Period;
                times = obj.Times;
                amplitude = obj.fp_modulation(2,:);
                amplitude_unc = obj.stderr_modulation(2,:);
                
                % exponential fit
                % order of params needs to be [ModelParams, Amp, Offset, Period]
                % re-order arguments here and ignore period
                model_fn = @(P, M, T) exp1D([0, M(1), 1./P(1), M(2)], T);
                init_model_params = [0.01];
                fixed_model_params = [0];
                model_lbs = [0];
                model_ubs = [inf];
                num_curve_params = 3;
                init_curve_params = [0.1, 0, period];
                fixed_curve_params = [0, 1, 1];
                curve_lbs = [0, 0, 0];
                curve_ubs = [1, 1, inf];
                
               [model_fit_params, model_std_errs, curve_fit_params, curve_fit_std_errs, chi_sqr, exitflag] ...
                = fit_model_simultaneous(times, amplitude, amplitude_unc,...
                model_fn,...
                init_model_params, fixed_model_params, model_lbs, model_ubs,...
                init_curve_params, fixed_curve_params, curve_lbs, curve_ubs, num_curve_params);    
                
                fighandle = obj.plotDiffuseModel(times,...
                    amplitude, amplitude_unc, period,...
                    model_fn, model_fit_params, model_std_errs,...
                    curve_fit_params, curve_fit_std_errs, chi_sqr);
                
                FpExp = [model_fit_params, curve_fit_params];
                FpSE = [model_std_errs, curve_fit_std_errs];
                
                if obj.Saving
                    % figure
                    fname = sprintf('%s_fit_decay_exponential.fig', obj.DescriptStr);
                    fpath = fullfile(save_dir, fname);
                    savefig(fighandle, fpath);
                    
                    % png
                    [f, n, ~] = fileparts(fpath);
                    png_fname = fullfile(f, [n, '.png']);
                    print(png_fname, '-dpng');
                end
                    
                % D/Gamma variant
                % Damped Osc 1, Gamma and D
                % can't get any k dependence from a single
                % wavevector, so may as well fit only the combinations
                % P = [Gamma, D, Amplitude]
                % dampedOsc1D([InitPos, InitVel, Omega0^2, DecayLen, Cx, Bg]
                model_fn_2 = @(P, M, T) dampedOsc1D([M(1), 0,...
                    sqrt( P(2) * P(1) ) * (2 * pi / M(3)),...
                    2 / P(1), 0, M(2)], T);
                num_model_params_2 = 2;
                init_model_params_2 = [6e-3, 10e-3];
                fixed_model_params_2 = [0, 0];
                model_lbs_2 = [0, 0];
                model_ubs_2 = [1, 1];
                
                num_curve_params_2 = 3;
                init_curve_params_2 = [0.1, 0, period];
                fixed_curve_params_2 = [0, 1, 1];
                curve_lbs_2 = [0, -inf, 0]; 
                curve_ubs_2 = [1, inf, inf];
                
                % initial fit
                [model_fit_params_2, model_std_errs_2, curve_fit_params_2, curve_fit_std_errs_2, chi_sqr_2, exitflag] ...
                    = fit_model_simultaneous(times, amplitude, amplitude_unc,...
                    model_fn_2,...
                    init_model_params_2, fixed_model_params_2, model_lbs_2, model_ubs_2,...
                    init_curve_params_2, fixed_curve_params_2, curve_lbs_2, curve_ubs_2, num_curve_params_2);    

                if use_offsets
                    Gamma = model_fit_params_2(1);
                    D = model_fit_params_2(2);
                    [offset, offset_unc] = obj.estimate_offset(times, amplitude, amplitude_unc, period, Gamma, D);
                    % refit with offset
                    init_curve_params_2(2) = offset;
                    fixed_curve_params_2(2) = 1;
                    [model_fit_params_2, model_std_errs_2, curve_fit_params_2, curve_fit_std_errs_2, chi_sqr_2, exitflag] ...
                        = fit_model_simultaneous(times, amplitude, amplitude_unc,...
                        model_fn_2,...
                        init_model_params_2, fixed_model_params_2, model_lbs_2, model_ubs_2,...
                        init_curve_params_2, fixed_curve_params_2, curve_lbs_2, curve_ubs_2, num_curve_params_2);    
                end
                Gamma = model_fit_params_2(1);
                D = model_fit_params_2(2);
                
                fighandle = obj.plotDiffuseModel(times,...
                    amplitude, amplitude_unc, period,...
                    model_fn_2, model_fit_params_2, model_std_errs_2,...
                    curve_fit_params_2, curve_fit_std_errs_2, chi_sqr_2);
                % add -Gamma/2 decaying exponential to curve
                ax = gca;
                hold on;
                t_interp = linspace(0, max(times), 300);
                exp_amp = curve_fit_params_2(1);
                plot(t_interp, exp1D([0, 1, 2/Gamma, 0], t_interp), 'r');
                
                obj.Fpdo = [model_fit_params_2, curve_fit_params_2];
                obj.SEdo = [model_std_errs_2, curve_fit_std_errs_2];

                
                if obj.Saving
                    % figure
                    fname = sprintf('%s_fit_decay_gamma_D.fig', obj.DescriptStr);
                    fpath = fullfile(save_dir, fname);
                    savefig(fighandle, fpath);
                    
                    % png
                    [f, n, ~] = fileparts(fpath);
                    png_fname = fullfile(f, [n, '.png']);
                    print(png_fname, '-dpng');
                    
                    % export results
                    fname = sprintf('%s_period=%0.1f_decay_versus_time.txt', obj.DescriptStr, period);
                    
                    data = [times', amplitude', amplitude_unc',...
                        (amplitude' - curve_fit_params_2(2)) / curve_fit_params_2(1), amplitude_unc' / curve_fit_params_2(1)];
                    names_cell = {'time(us)', 'amplitude(delta n_up)', 'amplitude_unc',...
                        'amplitude_normalized', 'amplitude_normalized_unc'};
                    
                    save_data_file(data, names_cell, '\t', '', save_dir, fname, 0);
                end
                
                % do boot-strap
                ntrials = 1000;
                if do_bootstrap
                    % two possible approaches to this error analysis: take
                    % dependent variables with errorbars, and select normally
                    % distributed points around the actual datapoint with
                    % standard deviation given by the errorbar. Fit many of
                    % these and look at the distribution of results.
                    % Other approach: first fit the points, then use the fitted
                    % values to generate the normally distributed points.
                    
                    [model_fit_params_bs, model_std_errs_bs, model_fit_params_distribution,...
                        curve_fit_params_bs, curve_std_errs_bs, curve_fit_params_distribution,...
                        chi_sqr_bs, num_failed_fits] ...
                        = bootstrap_model_simultaneous(times, amplitude, amplitude_unc,...
                        model_fn_2,...
                        model_fit_params_2, fixed_model_params_2, model_lbs_2, model_ubs_2,...
                        curve_fit_params_2, fixed_curve_params_2, curve_lbs_2, curve_ubs_2, num_curve_params_2,...
                        ntrials, 'fit');
                    
                    obj.Fpdo_bs = [model_fit_params_bs, curve_fit_params_bs];
                    obj.SEdo_bs = [model_std_errs_bs, curve_std_errs_bs];
                    
                    fighandle = obj.plotDiffuseModel(times,...
                        amplitude, amplitude_unc, period,...
                        model_fn, model_fit_params_bs, model_std_errs_bs,...
                        curve_fit_params_bs, curve_std_errs_bs, chi_sqr_bs);
                    
                    % look at distribution of model parameters
                    num_sets = 1;
                    
                    fighandle_model_dist = figure;
                    ncols = ceil(sqrt(num_model_params_2));
                    nrows = ceil(num_model_params_2 / ncols);
                    for kk = 1 : num_model_params_2
                        ax = subplot(nrows, ncols, kk);
                        h = histogram(model_fit_params_distribution(:, kk));
                        title(sprintf('Model Param %d', kk));
                    end
                    suptitle(sprintf('%s\n model parameters', obj.DescriptStr));
                    
                    % llok at distribution of curve parameters
                    fighandle_curve_dist = figure;
                    ncols = ceil(sqrt(num_curve_params_2 * num_sets));
                    nrows = ceil(num_curve_params_2 * num_sets / ncols);
                    for kk = 1 : num_curve_params_2 * num_sets
                        ax = subplot(nrows, ncols, kk);
                        h = histogram(curve_fit_params_distribution(:, kk));
                        title(sprintf('Curve Param %d', kk));
                    end
                    suptitle(sprintf('%s\n curve parameters', obj.DescriptStr));
                    
                    if obj.Saving 
                        % save figure
                        fig_fname = sprintf('%s_fit_decay_gamma_D_bootstrap.fig', obj.DescriptStr);
                        fpath = fullfile(save_dir, fig_fname);
                        savefig(fighandle, fpath);
                        % save png file
                        [f, n, ~] = fileparts(fpath);
                        png_fpath = fullfile(f, [n, '.png']);
                        print(png_fpath, '-dpng');
                        
                        % save histogram, model fit params
                        fig_fname = sprintf('%s_fit_decay_gamma_D_no_offset_bootstrap_model_param_histogram.fig', obj.DescriptStr);
                        fpath = fullfile(save_dir, fig_fname);
                        savefig(fighandle_model_dist, fpath);
                        % png
                        [f, n, ~] = fileparts(fpath);
                        png_fpath = fullfile(f, [n, '.png']);
                        print(png_fpath, '-dpng');
                        
                        % save histogram, curve fit params
                        fig_fname = sprintf('%s_fit_decay_gamma_D_no_offset_bootstrap_curve_param_histogram.fig', obj.DescriptStr);
                        fpath = fullfile(save_dir, fig_fname);
                        savefig(fighandle_curve_dist, fpath);
                        % png
                        [f, n, ~] = fileparts(fpath);
                        png_fpath = fullfile(f, [n, '.png']);
                        print(png_fpath, '-dpng');
                    end  
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % fit temperature to DQMC
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             export_data = 1;
            if ~isempty(obj.Density) && ~isempty(obj.Singles)
                [fighandle, FitParams,SErrs] = obj.fitSinglesVsDens(...
                    n_low_cutoff, 'DQMC',...
                    use_quantity_to_fit, temp_fit_init_p,[0,1],...
                    obj.Saving, save_dir);
                obj.Temp = FitParams(1);
                obj.TempUnc = SErrs(1);
                
                if obj.Saving
                    % figure
                    fname = sprintf('%s_fit_temp_DQMC.fig', obj.DescriptStr);
                    fpath = fullfile(save_dir, fname);
                    savefig(fighandle, fpath);
                    
                    %png
                    [f, n, ~] = fileparts(fpath);
                    png_fname = fullfile(f, [n, '.png']);
                    print(png_fname, '-dpng');
                end
                
            else
                fprintf('Skipped temperature fitting because not all necessary quantities were calculated\n');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % display azimuthal average
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isempty(obj.dfs)
                fighandle = obj.showAzAvg();
                if obj.Saving
                    % figure
                    fname = sprintf('%s_showAzAvg.fig', obj.DescriptStr);
                    fpath = fullfile(save_dir, fname);
                    savefig(fighandle, fpath);

                    % png
                    [f, n, ~] = fileparts(fpath);
                    png_fname = fullfile(f, [n, '.png']);
                    print(png_fname, '-dpng');
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % display all one-d averages
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fighandle = obj.showCuts();
            if obj.Saving
                % figure
                fname = sprintf('%s_showCuts.fig', obj.DescriptStr);
                fpath = fullfile(save_dir, fname);
                savefig(fighandle, fpath);
                
                % png
                [f, n, ~] = fileparts(fpath);
                png_fname = fullfile(f, [n, '.png']);
                print(png_fname, '-dpng');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % save structure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.Saving
                % folders should have been saved in the previous function.
                % This function should be independent of those as much as
                % possible.
                obj.saveStruct(save_dir, '', 0);
            end
            
        end
        
        function analyzeMultiWavelengths(obj, analysis_dir, script_file_pattern, analyze_data_folders)
            % analyzeMultiWavelengths(obj, analysis_dir, script_file_pattern)
            % fit results of analyzeDataFolders and analyzeSingleWavelength
            % to diffusion models etc using all wavevectors
            % For each analysis script of the form set*.m in a given folder
            % it, run it and run DiffuseSet with the variables derived from
            % the script. This performs the analysis for a single
            % k-vector at a given temperature. This function loops over all
            % such scripts in the folder, and therefore does the analysis
            % for the full temperature.
            % 
            % This function doesn't load any data into this class object.
            % It only saves the DiffuseSet objects. This data can be loaded
            % and processed using the obj.intialize function
            if ~exist('analysis_dir','var') || isempty(analysis_dir)
                analysis_dir = pwd;
            end
            
            if ~exist('script_file_pattern', 'var') || isempty(script_file_pattern)
                script_file_pattern = 'set*.m';
            end
            
            if ~exist('analyze_data_folders', 'var') || isempty(analyze_data_folders)
                analyze_data_folders = 1;
            end
            
            %load files
            AnalysisFiles = dir(fullfile(analysis_dir, script_file_pattern));
            
            %run each analysis file
            for ii = 1:length(AnalysisFiles)
                FPath = fullfile(analysis_dir, AnalysisFiles(ii).name);
                fprintf('Starting set %d of %d \n', ii, length(AnalysisFiles));
                fprintf('%s\n', AnalysisFiles(ii).name);
                %remove old variables, which should be defined in thescript
                clear DateCell Folders NormFolders Times PeriodGuess;
                clear ExcludePics ExcludePicsNorm;
                clear OnesFolder ThreesFolder ExcludeOnes ExcludeThrees;
                clear TempFitInitP SinglesFolder ExcludeSingles TempFitInitP;
                clear UseQuantityToFit CutoffR DisplayAll ComputeCorr;
                clear BinWidth CutoffR save_path BinEndPts;
                clear FixedPeriod
                %try
                    run(FPath);
                    if ~exist('NormFolders','var')
                        NormFolders = [];
                    end
                    if ~exist('ExcludePicsNorm','var')
                        ExcludePicsNorm = {};
                    end
                    if ~exist('DisplayAll','var')
                        DisplayAll = [];
                    end
                    if ~exist('BinWidth','var')
                        BinWidth = [];
                    end
                    if ~exist('ComputerCorr','var')
                        ComputeCorr = [];
                    end
                    if ~exist('CutoffR','var')
                        CutoffR = [];
                    end
                    if ~exist('n_low_cutoff','var')
                        n_low_cutoff = 0.5;
                    end
                    
                    if ~exist('TempFitInitP','var')
                        TempFitInitP = [];
                    end
                    
                    if ~exist('UseQuantityToFit','var')
                        UseQuantityToFit = [];
                    end
                    
                    if ~exist('save_path','var')
                        [~, fname, ~] = fileparts(FPath);
                        save_path = fullfile(analysis_dir, fname);
                    end
                    
                    if ~exist(save_path,'dir')
                        mkdir(save_path);
                    end
                    
                    if ~exist('BinEndPts', 'var')
                        BinEndPts = [];
                    end
                    
                    if ~exist('FixedPeriod', 'var')
                        FixedPeriod = 0;
                    end
                    % do analysis for many k-vector sets
                    
                    dset = DiffuseSet();
                    if analyze_data_folders
                        dset.analyzeDataFolders(DateCell, Folders, NormFolders, Times,...
                            ExcludePics, ExcludePicsNorm, ...
                            Cx, Cy, OnesFolder, ThreesFolder, ExcludeOnes,...
                            ExcludeThrees, BinWidth, ComputeCorr,...
                            CutoffR, SinglesFolder, ExcludeSingles,...
                            save_path, [])
                    else
                        saved_file_pattern = fullfile(save_path, '*_struct.mat');
                        saved_files = dir(saved_file_pattern);
                        if length(saved_files) == 1
                            a = load(fullfile(saved_files.folder, saved_files.name));
                            dset.loadStruct(a.AsStruct);
                        else
                            fprintf('Did not find any fields matching pattern %s, but analyze_data_folders set to zero. Skipped script %s',...
                                saved_file_pattern, FPath);
                            if ii == length(AnalysisFiles)
                                return
                            end
                        end
                        
                    end
                    dset.analyzeSingleWavelength(PeriodGuess, FixedPeriod,...
                                   n_low_cutoff, UseQuantityToFit, TempFitInitP, save_path)
                    
                    close all;
            end
        end
        
        function dset_stack = analyzeSingleTemp(obj, analysis_dir, min_lambda, max_amp,...
                save_results,...
                save_dir, summary_file_dir, load_all_folders,...
                wavelength_dir_pattern, file_pattern)
            % analyzeSingleTemp(obj, analysis_dir, min_lambda, max_amp,...
            %    save_results,...
            %    save_dir, summary_file_dir, load_all_folders,...
            %    wavelength_dir_pattern, file_pattern)
            %
            % Load saved DiffuseSet data and process it.
            %
            % Extract temperature and density information from these sets.
            % Fit modulation decay to various models. Export this data.
            %
            % analysis_dir is the directory where the save DiffuseSet
            % objects are stored.If this is not supplied, the current
            % directory is assumed.
            %
            % min_lambda is the minimum modulation wavelength to be
            % included in the analysis. Shorter wavelengths will be ignored
            % when fitting the hydrodynamic models.
            %
            % save_results
            %
            % save_dir
            %
            % summary_file_dir
            %
            % load_all_folders is a boolean specifying whether or not to
            % load all folder data for the DiffuseSet objects. This is not
            % necessary for fitting the hydrodynamic models, but for other
            % analysis it might be required.
            %
            
            if ~exist('analysis_dir', 'var') || isempty(analysis_dir)
                analysis_dir = pwd;
            end
            
            if ~exist('min_lambda', 'var') || isempty(min_lambda)
                min_lambda = 0;
            end
            
            if ~exist('max_amp', 'var') || isempty(max_amp)
                max_amp = inf;
            end
            
            if ~exist('save_results', 'var') || isempty(save_results)
                save_results = 0;
            end
            
            if ~exist('save_dir', 'var') || isempty(save_dir)
                save_dir = pwd;
            end
            
            if ~exist('summary_file_dir', 'var') || isempty(summary_file_dir)
                summary_file_dir = save_dir;
            end
            
            if ~exist('load_all_folders', 'var') || isempty(load_all_folders)
                load_all_folders = 0;
            end
            
            if ~exist('file_pattern', 'var') || isempty(file_pattern)
                file_pattern = '*_struct.mat';
            end
            
            if ~exist('wavelength_dir_pattern', 'var') || isempty(wavelength_dir_pattern)
                wavelength_dir_pattern = 'set*';
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %load all saved analysis files
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            analysis_files = []; %dir(fullfile(analysis_dir, file_pattern));
            
            % load DiffuseSets from subdirectories of analysis folder with
            % the correct pattern.
            subdirs_to_use = dir( fullfile( analysis_dir, wavelength_dir_pattern) );
            subdirs_to_use = subdirs_to_use([subdirs_to_use.isdir]);
            for ii = 1:length(subdirs_to_use)
                subdir = subdirs_to_use(ii);
                files_temp = dir( fullfile(analysis_dir, subdir.name, file_pattern) );
                analysis_files = vertcat(analysis_files, files_temp);
            end
            
            dset_stack = [];
            for ii = 1:length(analysis_files)
                fname = fullfile(analysis_files(ii).folder, analysis_files(ii).name);
                diffset = DiffuseSet();
                diffset.loadStruct(fname, '', load_all_folders);
                dset_stack = vertcat(dset_stack, diffset);
            end
            
            % remove any undesirable files...
            % cutoff wavelengths
            dset_stack = dset_stack( horzcat(dset_stack.Period) > min_lambda);
            % TODO: implement amplitude cutoff...
%             dset_stack = dset_stack()
            
            %Sort by period
            [~, order] = sort( horzcat(dset_stack.Period) );
            dset_stack = dset_stack(order);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extract temperature and density info from sets
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            periods = horzcat(dset_stack.Period);
            AllTemps = horzcat(dset_stack.Temp);
            AllTempsUnc = horzcat(dset_stack.TempUnc);
            AllPeriods = horzcat(dset_stack.Period);
            AllPeriodsUnc = horzcat(dset_stack.PeriodUnc);
            AllDens = horzcat(dset_stack.DensAvg);
            AllDensUnc = horzcat(dset_stack.DensAvgUnc);
            AllDensGrads = horzcat(dset_stack.DensGrads);
            AllDensGradsUnc = horzcat(dset_stack.DensGradsUnc);
            AllSineFitParams = horzcat(dset_stack.fp_modulation);
            AllSineFitErrs = horzcat(dset_stack.stderr_modulation);
            % model fits to single wavevectors
            all_fit_params_no_bg = vertcat(dset_stack.Fpdo);
            all_std_errs_no_bg = vertcat(dset_stack.SEdo);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute averages related to this set (Density, Temp, etc)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Temp = mean(AllTemps(AllTemps~=0));
            TempSD = std(AllTemps(AllTemps~=0));
            TempSDM = std(AllTemps(AllTemps~=0)) / ...
                sqrt(length(AllTemps(AllTemps~=0)));
            
            Dens = mean(AllDens(AllDens~=0));
            DensSD = std(AllDens(AllDens~=0));
            DensSDM = std(AllDens(AllDens~=0)) / ...
                sqrt(length(AllDens(AllDens~=0)));
            
            DescriptStr = sprintf('DiffuseSummary_T=%0.2f_n=%0.2f',...
                Temp, Dens);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Combined fits to different fit models (i.e. fitting all
            % k-vectors at fixed temperature)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            num_sets = length(dset_stack);
            times = {dset_stack.Times};
            sinusoid_fp = {dset_stack.fp_modulation};
            sinusoid_stderr = {dset_stack.stderr_modulation};
            slicing_fn = @(c) c(2, :);
            amplitudes = cellfun(slicing_fn, sinusoid_fp, 'UniformOutput', false);
            amplitude_uncs = cellfun(slicing_fn, sinusoid_stderr, 'UniformOutput', false);
            
            
            % Fit model one:
            % dt^2n + (Gamma + Gamma_2*k^2)*dtn + (Gamma*D + Omega_2^2*k^2)*k^2*n = 0
            % P = [Gamma, D, Gamma_2, Omega_2]
            % M = [Amplitude, Offset, Modulation Period]
            % dampedOsc1D([InitPos, InitVel, Omega0^2, DecayLen, Cx, Bg]
            model_fn = @(P, M, T) dampedOsc1D([M(1), 0,...
                sqrt( P(2) * P(1) + P(4)^2 * (2 * pi / M(3))^2 ) * (2 * pi / M(3)),...
                2 / ( P(1) + P(3) * (2 * pi / M(3))^2 ), 0, M(2)], T);
            num_model_params = 4;
            init_model_params = [7e-3, 20e-3, 0, 0];
            fixed_model_params = [0, 0, 1, 1];
            model_lbs = [0, 0, 0, 0];
            model_ubs = [100e-3, 100e-3, inf, inf];
            
            num_curve_params = 3;
            init_curve_params = [0.05 * ones(1, num_sets); 0 * ones(1, num_sets); reshape(periods, [1, num_sets])];
            init_curve_params = init_curve_params(:);
            fixed_curve_params = [zeros(1, num_sets); ones(1, num_sets); ones(1, num_sets)];
            fixed_curve_params = fixed_curve_params(:);
            curve_lbs = [zeros(1, num_sets); -1 * ones(1, num_sets); zeros(1, num_sets)];
            curve_lbs = curve_lbs(:);
            curve_ubs = [ones(1, num_sets); ones(1, num_sets); inf * ones(1, num_sets)];
            curve_ubs = curve_ubs(:);
            
            [model_fit_params, model_std_errs, curve_fit_params, curve_fit_std_errs, chi_sqr, exitflag] ...
                = fit_model_simultaneous(times, amplitudes, amplitude_uncs,...
                model_fn,...
                init_model_params, fixed_model_params, model_lbs, model_ubs,...
                init_curve_params, fixed_curve_params, curve_lbs, curve_ubs, num_curve_params);     
            
            fighandle = obj.plotDiffuseModel(times,...
                amplitudes, amplitude_uncs, periods,...
                model_fn, model_fit_params, model_std_errs,...
                curve_fit_params, curve_fit_std_errs, chi_sqr);
            
             % save results, model 1
             if save_results
                 analysis_fname_stem = sprintf('minlambda=%0.1f_maxamp=%0.2f',...
                     min_lambda, max_amp);
                 % export results to text file
                 text_fname = sprintf('expt_hydro_fit_summary.txt');
                 export_dat = [Temp, TempSD, TempSDM, Dens, DensSD, DensSDM,...
                     obj.interleave_vects(model_fit_params, model_std_errs)];
                 titles_model = {'Gamma(Mhz)', 'GammaUnc(MHz)', 'D(latt^2*MHz)',...
                     'DUnc(latt^2*MHz)','Gamma_2(latt^2*MHz)','Gamma_2Unc(latt^2*MHz)',...
                     'omega_2(latt^2*Mhz)', 'omega_2Unc(latt^2*MHz)'};
                 export_titles = cat(2,{'Temp','TempSD','TempSDM','Dens','DensSD','DensSDM'}, ...
                     titles_model);           
                 append_dat = exist(fullfile(summary_file_dir, text_fname), 'file');
                 if append_dat
                     export_titles = cell(1, length(export_titles));
                 end
                 save_data_file(export_dat, export_titles, '\t', '', summary_file_dir, text_fname, append_dat);

                 % save figure
                 fig_fname = sprintf('%s_hydro_fit_%s.fig', DescriptStr, analysis_fname_stem);
                 fpath = fullfile(save_dir, fig_fname);
                 savefig(fighandle, fpath);
                 
                 % save png file
                 % png
                 [f, n, ~] = fileparts(fpath);
                 png_fpath = fullfile(f, [n, '.png']);
                 print(png_fpath, '-dpng');
             end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % bootstrap model 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            do_bootstrap = 1;
            ntrials = 1000;
            if do_bootstrap
                % two possible approaches to this error analysis: take
                % dependent variables with errorbars, and select normally
                % distributed points around the actual datapoint with
                % standard deviation given by the errorbar. Fit many of
                % these and look at the distribution of results.
                % Other approach: first fit the points, then use the fitted
                % values to generate the normally distributed points.
                
                [model_fit_params_bs, model_std_errs_bs, model_fit_params_distribution,...
                curve_fit_params_bs, curve_std_errs_bs, curve_fit_params_distribution,...
                chi_sqr_bs, num_failed_fits] ...
                = bootstrap_model_simultaneous(times, amplitudes, amplitude_uncs,...
                model_fn,...
                init_model_params, fixed_model_params, model_lbs, model_ubs,...
                init_curve_params, fixed_curve_params, curve_lbs, curve_ubs, num_curve_params,...
                ntrials, 'fit');
            
               fighandle = obj.plotDiffuseModel(times,...
                amplitudes, amplitude_uncs, periods,...
                model_fn, model_fit_params_bs, model_std_errs_bs,...
                curve_fit_params_bs, curve_std_errs_bs, chi_sqr_bs);
            
                % look at distribution of model parameters
                fighandle_model_dist = figure;
                ncols = ceil(sqrt(num_model_params));
                nrows = ceil(num_model_params / ncols);
                for kk = 1 : num_model_params
                    ax = subplot(nrows, ncols, kk);
                    h = histogram(model_fit_params_distribution(:, kk));
                    title(sprintf('Model Param %d', kk));
                end
                suptitle(sprintf('%s\n model parameters', DescriptStr));
                
                % llok at distribution of curve parameters
                fighandle_curve_dist = figure;
                ncols = ceil(sqrt(num_curve_params * num_sets));
                nrows = ceil(num_curve_params * num_sets / ncols);
                for kk = 1 : num_curve_params * num_sets
                    ax = subplot(nrows, ncols, kk);
                    h = histogram(curve_fit_params_distribution(:, kk));
                    title(sprintf('Curve Param %d', kk));
                end
                suptitle(sprintf('%s\n curve parameters', DescriptStr));
                
                if obj.Saving
                    analysis_fname_stem = sprintf('minlambda=%0.1f_maxamp=%0.2f_bootstrap',...
                        min_lambda, max_amp);
                    
                    text_fname = sprintf('expt_hydro_fit_bootstrap_summary.txt');
                    export_dat = [Temp, TempSD, TempSDM, Dens, DensSD, DensSDM,...
                        obj.interleave_vects(model_fit_params_bs, model_std_errs_bs)];
                    titles_model = {'Gamma(Mhz)', 'GammaUnc(MHz)', 'D(latt^2*MHz)',...
                        'DUnc(latt^2*MHz)','Gamma_2(latt^2*MHz)','Gamma_2Unc(latt^2*MHz)',...
                        'omega_2(latt^2*Mhz)', 'omega_2Unc(latt^2*MHz)'};
                    export_titles = cat(2,{'Temp','TempSD','TempSDM','Dens','DensSD','DensSDM'}, ...
                        titles_model);

                    append_dat = exist(fullfile(summary_file_dir, text_fname), 'file');
                    if append_dat
                        export_titles = cell(1, length(export_titles));
                    end
                    save_data_file(export_dat, export_titles, '\t', '', summary_file_dir, text_fname, append_dat);

                    % save figure
                    fig_fname = sprintf('%s_hydro_fit_%s.fig', DescriptStr, analysis_fname_stem);
                    fpath = fullfile(save_dir, fig_fname);
                    savefig(fighandle, fpath);
                    % save png file
                    [f, n, ~] = fileparts(fpath);
                    png_fpath = fullfile(f, [n, '.png']);
                    print(png_fpath, '-dpng');
                    
                    % save histogram, model fit params
                    fig_fname = sprintf('%s_%s_model_param_histogram.fig',...
                                        DescriptStr, analysis_fname_stem);
                    fpath = fullfile(save_dir, fig_fname);
                    savefig(fighandle_model_dist, fpath);
                    % png
                    [f, n, ~] = fileparts(fpath);
                    png_fpath = fullfile(f, [n, '.png']);
                    print(png_fpath, '-dpng');
                    
                    % save histogram, curve fit params
                    fig_fname = sprintf('%s_%s_curve_param_histogram.fig',...
                                        DescriptStr, analysis_fname_stem);
                    fpath = fullfile(save_dir, fig_fname);
                    savefig(fighandle_curve_dist, fpath);
                    % png
                    [f, n, ~] = fileparts(fpath);
                    png_fpath = fullfile(f, [n, '.png']);
                    print(png_fpath, '-dpng');  
                end
            
            end
            
            %compare sets
            if length(dset_stack) > 1 && load_all_folders
                ds = dset_stack(1);
                ds.compareDiffuseSets([dset_stack(2:end)], save_results, save_dir);
            end
        end
        
        function analyzeEverything(obj, main_dir, temp_dir_pattern, script_file_pattern, reanalyze_data_folders)
            %TODO: integrate this here!
            % Do full analysis of dataset using DiffuseSet and variable defining
            % script files. Assume all sub-directories of
            % this main directory house all the data for a given temperature.
            %
            % DiffuseSetSummary will enter each directory. It will then run DiffuseSet
            % for each k-vector. All the excluded pictures etc for a given k-vector are
            % defined in a matlab script file of the form set*.m. As the analysis is
            % done, most of the results will be saved in a subdirectory set*.
            %
            % Before performing a new analysis, run "clean_old_analysis" to move the
            % older analysis to the 'hide' subfolders.

            if ~exist('main_dir', 'var') || isempty(main_dir)
                main_dir = pwd;
            end
            
            if ~exist('temp_dir_pattern', 'var') || isempty(temp_dir_pattern)
                temp_dir_pattern = 'data_*';
            end
            
            if ~exist('wavelength_dir_pattern', 'var') || isempty(wavelength_dir_pattern)
                wavelength_dir_pattern = 'set*';
            end
            
            if ~exist('script_file_pattern', 'var') || isempty(script_file_pattern)
                script_file_pattern = 'set*.m';
            end
            
            if ~exist('reanalyze_data_folders', 'var') || isempty(reanalyze_data_folders)
                reanalyze_data_folders = 0;
            end
            
            %this option will redo analsis for every folder.
            %Don't use this if you only want to refit hydro params
            min_lambda = 0;
            max_amp = 1;
            save_results = 1;
            
            % get folders to use
            dir_form = fullfile(main_dir, temp_dir_pattern);
            all_files = dir(dir_form);
            % only want directories
            sub_dirs = all_files([all_files.isdir]);
            % remove any hidden directories
            ishidden = @(c) ~strcmp(c(1), '.');
            sub_dirs = sub_dirs(cellfun(ishidden, {sub_dirs.name}));
            
            summary_file_dir = main_dir;
            
            % loop through immediate subdirectories of the directory containing this
            % file. These directories should each contain all the data for a fixed
            % temperature. Run analysis.
            %
            
            for ii = 1:length(sub_dirs)
                fprintf('Directory %s\n', sub_dirs(ii).name);
                % check that our directory contains necessary components for analysis
                if ~isempty(dir(fullfile(sub_dirs(ii).folder, sub_dirs(ii).name, script_file_pattern)))
                    
                    sub_dir_path = fullfile(sub_dirs(ii).folder, sub_dirs(ii).name);
                    
                    % process all folders and save results
                    if reanalyze_data_folders
                        %TODO
%                         obj.analyzeMultiWavelengths();
                    end
                    
                    obj.analyzeSingleTemp(sub_dir_path, min_lambda, max_amp,...
                        save_results, sub_dir_path, summary_file_dir, 0,...
                        wavelength_dir_pattern, '*_struct.mat');
                        
                end
                close all;
            end
        end
        
        function displayFolder(obj, folder_path)
            %display all figures in a certain folder path
            %
            % folder_path
            %
            
            if ~exist('folder_path') || isempty(folder_path)
                folder_path = pwd;
            end
            
            figs = dir(fullfile(folder_path, '*.fig'));
            for ii = 1:length(figs)
                open(fullfile(figs(ii).folder, figs(ii).name));
            end
        end
        
        function clean_old_analysis(obj, root_dir_path, data_folder_prefix, hide_data_folders)
            % This function searches subdirectories for analysis files and
            % moves them into a dated folder hidden from future analysis.
            % This way future analysis will not destroy previous analysis
            % data.
            PatternsToHide = {'*.fig', '*_struct.mat', '*.txt', '*.dat', '*.png'};
            
            if ~exist('hide_data_folders', 'var')
                hide_data_folders = 0;
            end
            
            if hide_data_folders
                PatternsToHide{end + 1} = '*.mat';
            end
            
            hide_subdir_name = 'hide';
            
            if ~exist('dir_path', 'var') || isempty(root_dir_path)
                root_dir_path = pwd;
            end
            
            if ~exist('data_folder_prefix', 'var')
                data_folder_prefix = 'data_*';
            end
            
%             if ~exist('include_root_dir', 'var')
%                 include_root_dir = 1;
%             end
            
            dir_form = fullfile(root_dir_path, data_folder_prefix);
            
            % get all files
            all_files = dir(dir_form);
            % restrict to directories
            directories = all_files([all_files.isdir]);
            folders = {directories.folder};
            names = {directories.name};
            
            % first remove any hidden files
            fun = @(c) ~strcmp(c(1), '.');
            indices_hidden_folders = cellfun(fun, names);
            
            folders = folders(indices_hidden_folders);
            names = names(indices_hidden_folders);
            
            % now get full directory names
            dir_paths = cell(1, length(folders));
            for ii = 1:length(folders)
                dir_paths{ii} = fullfile(folders{ii}, names{ii});
            end
            
            include_root_dir = 1;
            if include_root_dir
                dir_paths{end + 1} = root_dir_path;
                folders{end + 1} = root_dir_path;
                names{end + 1} = '.';
            end
            
            for ii = 1:length(dir_paths)
                % avoid 'hidden' directories
                curr_dir = dir_paths{ii};
                
                fprintf('Directory %s\n', curr_dir);
                files = [];
                for jj = 1:length(PatternsToHide)
                    pattern = fullfile(curr_dir, PatternsToHide{jj});
                    files = vertcat(files, dir(pattern));
                end
                
                % create hide subdirectory if it does not exist
                hide_dir = fullfile(curr_dir, hide_subdir_name);
                if ~exist(hide_dir, 'dir')
                    mkdir(hide_dir);
                end
                
                % create subdirectory labelled by todays date and time
                today_dir = fullfile(hide_dir, datestr(now, 'yyyy-mm-dd_HH_MM'));
                if ~exist(today_dir, 'dir')
                    mkdir(today_dir);
                end
                
                % hide files
                for jj = 1:length(files)
                    loc = fullfile(files(jj).folder, files(jj).name);
                    
                    dest = fullfile(today_dir, files(jj).name);
                    movefile(loc, dest);
                    %delete(fullfile(files(jj).folder,files(jj).name));
                    fprintf('Moved %s\n to %s\n',loc, dest);
                end
            end
        end
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fitting and display functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function fighandle = showCuts(obj)
            % showCuts(obj, save_results, save_dir)
            % show "cuts", which are the bin averages along the direction
            % orthogonal to the density modulation.
            % 
            % save_results is a boolean specifying whether or to save the 
            % figure
            % 
            % save_dir is the directory these results will be saved in 
            %
            % fname is the name of exported data figure
                
            fighandle = figure('name', 'Show Cuts');
            errorbar(obj.AllPos, obj.AllCuts, obj.AllUnc);
            grid on;
            xlabel('Pos')
            ylabel('Avg')
            ylim([0,1]);
            title(obj.DescriptStr);
            Leg = {};
            for ii = 1:length(obj.Times)
                Leg = cat(1, Leg, sprintf('%dus', obj.Times(ii)));
            end
            legend(Leg);
            
        end
        
        function fighandle = showAzAvg(obj)
            % showAzAvg(obj, save_results, save_dir, fname)
            % 
            % Show the azimuthal average of all DataFolder sets defined in
            % the Folders, Ones, Threes, and Singles variables.
            % Normalization folders are not shown (TODO show these).
            % Also show useful sums and differences of these sets
            % e.g. density.
            %
            % save_results
            %
            % save_dir
            %
            % fname
            
            %TODO: add display of normalization states?
            fighandle = figure('name', 'AzAvg');
            nrows = 2;
            ncols = 2;
            
            % plot azimuthal avgs, modulation data
            subplot(nrows, ncols, 1);
            Leg = {};
            num_times = length( obj.dfs_index_folders);
            for ii = 1:num_times
                errorbar(obj.BinAvgDist{ii}, obj.occ_vstime{ii}, obj.occ_vstime_unc{ii});
                hold on;
                Leg{ii} = sprintf('%0.0f us', obj.Times(ii));
            end
            grid on;
            ylim([0,1]);
            legend(Leg);
            xlabel('Position (Latt Sites)');
            title('Az Avg, Modulation Times');
            
            
            % special states
            sp_sets = {obj.dfs(obj.dfs_index_normfolders(1)),...
                obj.Singles, obj.Ones, obj.Threes};
            s1 = sprintf('normalization = %03d', obj.dfs_index_normfolders(1));
            s2 = sprintf('<n_s> = %03d', obj.dfs_index_singles);
            s3 = sprintf('<n_1> = %03d', obj.dfs_index_ones);
            s4 = sprintf('<n_3> = %03d', obj.dfs_index_threes);
            leg = {s1, s2, s3, s4};
            %if any of these are empty, remove them
            fn = @(x) ~isempty(x);
            inds = cellfun(fn,sp_sets);
            sp_sets = sp_sets(inds);
            leg = leg(inds);
            
            % az avg, special folders
            subplot(nrows,ncols,2)
            for ii = 1:length(sp_sets)
                set = sp_sets{ii};
                errorbar(set.BinAvg,set.Occs_AzAvg,set.Occs_AzAvgUnc);
                hold on;
            end
            if ~isempty(obj.Ones)&&~isempty(obj.Threes)
                errorbar(obj.Threes.BinAvg,obj.Density,obj.DensityUnc);
                legend(cat(2,leg,'<n>'));
            else
                legend(leg);
            end
            
            grid on;
            ylim([0,1]);
            xlabel('Position (Latt Sites)');
            title('Az Avg, Static');
            
            % plot correlators, special folders
            subplot(nrows,ncols,3)
            if ~isempty(sp_sets)
                hold on;
                CorrIndices = [0,1];
                for ii = 1:length(sp_sets)
                    obj2 = sp_sets{ii};
                    NNCorr = squeeze( obj2.Density_Corr_AzAvg(...
                        obj2.CenterIndex_CorrMatrix + CorrIndices(1),...
                        obj2.CenterIndex_CorrMatrix + CorrIndices(2), :) );
                    NNCorrUnc = squeeze(obj2.Density_Corr_AzAvgUnc(...
                        obj2.CenterIndex_CorrMatrix + CorrIndices(1),...
                        obj2.CenterIndex_CorrMatrix + CorrIndices(2), :));
                    errorbar(obj2.RadialPos, NNCorr, NNCorrUnc);
                end
                grid on;
                legend(leg);
                xlabel('Lattice Sites')
                title('Correlators');
            end
            
            % plot atom numbers
            subplot(nrows,ncols,4)
            hold on;
            for ii = 1:length(obj.dfs)
                plot(obj.dfs(ii).PictureNumberInFolder,...
                    obj.dfs(ii).AtomNumbers, 'o');
            end
            grid on;
            ax = gca;
            ax.YLim(1) = 0;
            ax.XLim(1) = 1;
            xlabel('Shot');
            ylabel('Atom Num');
            
            fn = @(x) {sprintf('Folder %d',x)};
            legend(arrayfun(fn, obj.dfs_folder_numbers));
            title('Atom Num Vs. Shot');
            
            suptitle(sprintf('%s',obj.DescriptStr));
        end
               
        function [fighandle, Fp,SE] = fitSinglesVsDens(...
                obj, nLowCutoff, Method,...
                UseQuantityToFit, InitP, FixedP,...
                export_data, save_dir)
            % TODO implement this using fit_dqmc_ns_corrupup
            % [fighandle, Fp,SE] = fitSinglesVsDens(...
            %    obj, nLowCutoff, Method,...
            %    UseQuantityToFit, InitP, FixedP,...
            %    save_results, save_dir, export_data) 
            %
            % Trap free fit of singles density and ones correlator versus
            % total density
            %
            % nLowCutoff is the lowest density to be included in the fit.
            % The DQMC data set does not have low density data for all
            % temperatures
            %
            % Method may be "DQMC" or "HTSE". HTSE ignores the correlator
            % data. Not sure if HTSE mode is working currently.
            %
            % UseQuantityToFit = [UseDensity,UseCorrelator], booleans
            % specifying which quantities to use in the fitting
            %
            % InitP = [Temp, ImagingFidelity]
            %
            % FixedP = [FixTemp, FixImagingFidelity]. Typically this should
            % be set to [0, 1] with ImagingFidelity = 0.97.
            %
            % save_results
            %
            % save_dir
            %
            % export_data. Boolean specifying whether or not to export
            % experimental data and fit data as text files. Exports singles
            % density and ones correlator versus density.
            %
            % Fp = [Temp]
            % 
            % SE = [TempUnc]
            
            if ~exist('nLowCutoff','var')
                nLowCutoff = 0.66;
            end
            
            if ~exist('nHighCutoff','var')
                nHighCutoff = InitP(2);
            end
            
            if ~exist('Method','var')
                Method = 'DQMC';
            end
            
            if ~exist('UseQuantityToFit','var')
                UseQuantityToFit = [1,0];
            end
            
            if isequal(UseQuantityToFit,[0,0])
                UseQuantityToFit = [1,0];
                disp('UseQuantityToFit cannot be [0,0]. Set to [1,0]');
            end

            if ~exist('save_dir','var')
                save_dir = pwd;
            end
            
            %experimental data
            n = transpose(obj.Density);
            nerr = transpose(obj.DensityUnc);
            nerr = nerr(n>nLowCutoff & n<nHighCutoff);
            
            ns = transpose(obj.Singles.Occs_AzAvg);
            nserr = transpose(obj.Singles.Occs_AzAvgUnc);
            nserr = nserr(n>nLowCutoff & n<nHighCutoff);
            ns = ns(n>nLowCutoff & n<nHighCutoff);
            
            try
                CenInd = obj.Ones.CenterIndex_CorrMatrix;
                cup = transpose(squeeze(obj.Ones.Density_Corr_AzAvg(CenInd+1,CenInd,:))); %obj.Ones.
                cuperr = transpose(squeeze(obj.Ones.Density_Corr_AzAvgUnc(CenInd+1,CenInd,:)));
                cup = cup(n>nLowCutoff & n<nHighCutoff);
                cuperr = cuperr(n>nLowCutoff & n<nHighCutoff);
                
                cdn = transpose(squeeze(obj.Ones.Density_Corr_AzAvg(CenInd+1,CenInd,:)));
                cdnerr = transpose(squeeze(obj.Ones.Density_Corr_AzAvgUnc(CenInd+1,CenInd,:)));
                cdn = cdn(n>nLowCutoff & n<nHighCutoff);
                cdnerr = cdnerr(n>nLowCutoff & n<nHighCutoff);
                
                cavg = 0.5*(cup+cdn);
                cavgerr = 0.5*sqrt(cuperr.^2+cdnerr.^2);
            catch
                %if there is some problem here, do not use the correlators
                %for fitting.
                cavg = zeros(size(ns));
                cavgerr = ones(size(ns));
                UseQuantityToFit(2) = 0;
            end
            
            n = n(n>nLowCutoff & n<nHighCutoff);
            
            %Fitting...theoretical data...
            if strcmp(Method,'DQMC')
                %Must insure that the fitting functions produce shape 1 x length(n)
                DQMC = load(obj.DQMCPath);
                SinglesFn = @(P,n) transpose(squeeze(DQMC.DQMC_Grid.nsfn_n_T(n,P(1))));
                CupFn = @(P,n) transpose(squeeze(DQMC.DQMC_Grid.cupupNNfn_n_T(n,P(1))));
                CdnFn = @(P,n) transpose(squeeze(DQMC.DQMC_Grid.cupupNNfn_n_T(n,P(1))));
                try
                    DQMCSignFn = @(P,n) transpose(squeeze(DQMC.DQMC_Grid.Signfn_n_T(n,P(1))));
                catch
                    DQMCSignFn = @(P,n) zeros(size(n));
                end
                %                 DQMCSignFn = @(P,n) zeros(size(n));
            elseif strcmp(Method,'HTSE')
                hte = load('hte_interp.mat');
                SinglesFn = @(P,n) squeeze(hte.hte_interp.nsfn_n_T_U(n,P(1),8));
                CupFn = @(P,n) zeros(size(n));
                CdnFn = @(P,n) zeros(size(n));
                DQMCSignFn = @(P,n) zeros(size(n));
            else
                disp('Method is not supported. Fit functions set to zero')
                SinglesFn = @(P,n) zeros(size(n));
                CupFn = @(P,n) zeros(size(n));
                CdnFn = @(P,n) zeros(size(n));
                DQMCSignFn = @(P,n) zeros(size(n));
            end
            
            
            FgFn_ns = @(P,n) ((n/P(2)) - 0.5*(n/P(2)).^2)*P(2);
            %             FitFn_ns = @(P,n) SinglesFn(P,n/P(3))*P(2);
            FitFn_ns = @(P,n) SinglesFn(P,n/P(2))*P(2);
            DiffFn_ns = @(P) (FitFn_ns(P,n)-ns)./nserr;
            %             FitFn_cavg = @(P,n) CupFn(P,n/P(3))*P(3)^2;
            FitFn_cavg = @(P,n) CupFn(P,n/P(2))*P(2)^2;
            DiffFn_cavg = @(P) (FitFn_cavg(P,n)-cavg)./cavgerr;
            
            %Select what to fit. Either singles density, correlator, or
            %both
            if isequal(UseQuantityToFit,[1,1])
                FullDiffFn = @(P) [DiffFn_ns(P),DiffFn_cavg(P)];
            elseif isequal(UseQuantityToFit,[1,0])
                FullDiffFn = @(P) DiffFn_ns(P);
            elseif isequal(UseQuantityToFit,[0,1])
                FullDiffFn = @(P) DiffFn_cavg(P);
            else
                exception
            end
            
            %implement ability to fix parameters
            %Params = [Temp,Imaging Singles (i.e. both spin states) Efficiency, Imaging One Spin State Efficiency]
            if ~exist('InitP','var')
                InitP = [8,0.97];
            end
            if ~exist('FixedP','var')
                FixedP = [0,1];
            end
            FullDiffFn = @(P) FullDiffFn(P.*(1-FixedP) + InitP.*FixedP);
            
            LowBounds = [0,0.9];
            UpBounds = [15,1];
            
            opts1=  optimset('display','off');
            [Fp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                lsqnonlin(FullDiffFn,InitP,LowBounds,UpBounds,opts1);
            try
                CI = nlparci(Fp,residual,'jacobian',jacobian);
                SE = transpose(CI(:,2)-CI(:,1))/3.92;
            catch
                SE = nan(size(Fp));
            end
            
            %plotting
            fighandle = figure('name','DQMC Fit');
            nrows = 1;
            ncols = 3;
            subplot(nrows,ncols,1)
            ph1 = errorbar(n, ns, nserr, nserr, nerr, nerr);
            hold on;
            nInterp = linspace(0, 1, 100);
            ph2 = plot(nInterp, FitFn_ns(Fp, nInterp), 'r');
            plot(nInterp,FitFn_ns([Fp(1) + SE(1), Fp(2)], nInterp), 'r--');
            plot(nInterp,FitFn_ns([Fp(1) - SE(1),Fp(2)], nInterp), 'r--');
            ph3 = plot(nInterp, FgFn_ns(Fp,nInterp),'g');
            xlabel('<n>')
            ylabel('<n^s>');
            grid on;
            ylim([0,1]);
            legend([ph1, ph2, ph3], {'Expt', Method, 'NonIntFG'})
            title('n^s Vs. n');
            
            subplot(nrows, ncols, 2)
            ph1 = errorbar(n, cavg, cavgerr, cavgerr, nerr, nerr);
            hold on;
            ph2 = plot(nInterp, FitFn_cavg(Fp,nInterp), 'r');
            plot(nInterp,FitFn_cavg([Fp(1) + SE(1), Fp(2)], nInterp), 'r--');
            plot(nInterp,FitFn_cavg([Fp(1) - SE(1), Fp(2)], nInterp), 'r--');
            xlabel('<n>');
            ylabel('Corr(0,1), Avg');
            grid on;
            legend([ph1, ph2], {'Expt', Method});
            title('NN Corr Vs. n');
            
            try
                subplot(nrows, ncols, 3);
                plot(nInterp, DQMCSignFn(Fp,nInterp), 'r');
                xlabel('<n>');
                ylabel('Avg Sign');
                grid on;
                title('DQMC Sign Vs. n');
                ylim([0, 1]);
            catch err
                disp(err.message);
            end
            
            ttl = sprintf('%s \n T = %0.2f +/- %0.2f t \n Efficiency = %0.2f +/- %0.2f (applied to theory fit, not expt data) \n Used Density = %d, Used Correlator = %d',...
                obj.DescriptStr, Fp(1), SE(1), Fp(2), SE(2),...
                UseQuantityToFit(1), UseQuantityToFit(2));
            suptitle(ttl);
            
%             if save_results
%                 fname = sprintf('%s_fit_temp_DQMC.fig', obj.DescriptStr);
%                 fpath = fullfile(save_dir, fname);
%                 savefig(fighandle, fpath);
%             end
            
            if export_data
                fit_params_string = sprintf('#T=%f+/-%f_Fidelity=%f+/-%f',...
                    Fp(1), SE(1),Fp(2), SE(2));
                
                %raw data
                fname = sprintf('%s_temp_DQMC_exptdat.txt', obj.DescriptStr);
                datarray = [n', nerr', ns', nserr', cavg', cavgerr'];
                names_cell = {'<n>', '<n>unc',  '<n^s>', '<n^s>unc',...
                    '<n_sn_s>_c', '<n_sn_s>_c_unc'};
                save_data_file(datarray, names_cell, '\t',...
                    fit_params_string, save_dir, fname, 0);
                
                %fit results
                fname = sprintf('%s_temp_DQMC_fitdat.txt', obj.DescriptStr);
                datarray = [nInterp', FitFn_ns(Fp,nInterp)',...
                    FitFn_cavg(Fp,nInterp)'];
                names_cell2 = {'<n>', '<n^s>', '<n_sn_s>_c'};
                save_data_file(datarray, names_cell2, '\t',...
                    fit_params_string, save_dir, fname, 0);
            end
        end
        
        function [fighandle, fp_modulation, stderr_modulation,...
                AllCuts_NoBgPos, AllCuts_NoBg, AllCuts_NoBgUnc] = ...
                fitCuts(obj, PeriodGuess, PeriodFixed, Times, XPositions, ...
                FitCuts, FitCutsUnc, NormCuts, NormCutsUnc,...
                save_dir, export_results, export_fname_stem)
            % [fp_modulation, stderr_modulation, AllCuts_NoBgPos, AllCuts_NoBg, AllCuts_NoBgUnc] = ...
            %    fitCuts(obj, PeriodGuess, PeriodFixed, Times, XPositions, ...
            %    FitCuts, FitCutsUnc, NormCuts, NormCutsUnc,...
            %    display_results,...
            %    save_results, save_dir, export_results)
            %
            % Given a collection of density versus position data taken at a
            % number of times, fit the spatial sinusoid patterns. Determine
            % the amplitude of the modulation at each time. Fix the phase
            % to the phase of the initial modulation.
            %
            % PeriodGuess. Initial guess for period of modulation.
            %
            % PeriodFixed. Boolean, whether or not to force fit to use
            % PeriodGuess as the period.
            %
            % Times. Times associated with each 'cut'
            %
            % XPositions. Mean x-position of bins along direction of 
            % modulation.
            %
            % FitCuts. Mean density of bins along direction of modulation.
            %
            % FitCutsUnc. Uncertainty in mean density of bins.
            %
            % NormCuts. Mean density of bins along direction of modulation
            % used to normalize FitCuts.
            %
            % NormCutsUnc. Uncertainties associated with NormCuts.
            %
            % save_dir. Directory where figures and data are to be saved.
            % What is saved is determined by save_results and
            % export_results
            %
            % export_results. Boolean specifying whether or not to export
            % XPositions versus FitCuts to text files.
            %
            % TODO try simultaneous fit with same phase for all sets,
            % instead of fixing phase from zero time set.
            
            if ~exist('Times','var') || isempty(Times)
                Times = 1:size(FitCuts, 2);
            end
            
            if ~exist('XPositions','var') || isempty(XPositions)
                  XPositions = 1:size(FitCuts, 1);
            end
                
            if ~exist('NormCuts','var') || isempty(NormCuts)
                NormCuts = zeros(size(FitCuts));
                NormCutsUnc = zeros(size(FitCutsUnc));
            end    
            
            % TODO get rid of this setting...can implement it by either
            % supplying or not supplying NormCuts
            
            if ~exist('save_dir','var') || isempty(save_dir)
                save_dir = pwd;
            end
            
            if ~exist('export_results', 'var') || isempty(export_results)
                export_results = 0;
            end
            
            if ~exist('export_fname_stem', 'var') || isempty(export_fname_stem)
                export_fname_stem = sprintf('%s_1d_averages', obj.DescriptStr);
            end
            
            % begin analysis
            ntimes = size(FitCuts, 2);
            nparams = 4;
            
            fp_modulation = zeros(nparams, ntimes);
            stderr_modulation = zeros(nparams, ntimes);
            AllFitFnHandles = cell([1, ntimes]);
            
            AllXs = [];
            AllYs = [];
            AllYErrs = [];
            
            for ii = 1:ntimes
                %expt data for this time
                X = XPositions(:,ii);
                Y = FitCuts(:,ii) - NormCuts(:,ii);
                Err = sqrt(FitCutsUnc(:,ii).^2 + NormCutsUnc(:,ii).^2);
                
                if ii == 1
                    InitP = [1/PeriodGuess, 0.1, 0, mean( Y(~isnan(Y)) )];
                    FixedP = [PeriodFixed, 0, 0, 0];
                else
                    InitP = [fp_modulation(1,1), -0.2, fp_modulation(3, 1), mean( Y(~isnan(Y)) )];
                    FixedP = [1, 0, 1, 0];
                end
                
                lbs = [0, -1, -inf, -inf];
                ubs = [1, 1, inf, inf];
                
                %remove nans and fit
                X = X(~isnan(Y));
                Err = Err(~isnan(Y));
                Err(Err == 0) = 1;
                Y = Y(~isnan(Y));
                [Fp, ~, FFH, SE, red_chi_sqr] = fit1D(X, Y, 1./Err.^2,...
                    {'sinusoid1D'},...
                    InitP, FixedP, lbs, ubs, []);
                
                %store results
                fp_modulation(:, ii) = Fp;
                stderr_modulation(:, ii) = SE;
                AllXs = cat(2, AllXs, X);
                AllYs = cat(2, AllYs, Y);
                AllYErrs = cat(2, AllYErrs, Err);
                AllFitFnHandles{ii} = FFH;
            end
            
            %ensure first Fit param is positive. change phase if
            %necessary
            if fp_modulation(2, 1) < 0
                fp_modulation(2, :) = -fp_modulation(2, :);
                fp_modulation(3, :) = fp_modulation(3, :) + pi;
            end
            
            AllCuts_NoBgPos = AllXs(:, 1);
            AllCuts_NoBg = AllYs;
            AllCuts_NoBgUnc = AllYErrs;
            
            xinterp = linspace(min(AllXs(:)), max(AllXs(:)), 100);
            
            %Display cuts with fits. Each time is offset vertically from
            %the others for clarity.
            fighandle = figure('name','Fit Cuts');
                
            % plot all curves offset vertically
                subplot(1, 2, 1);
                Offsets = -0.5 * ((1:ntimes) - 1);
                for ii = 1:ntimes
                    Offset = -0.5 * (ii - 1);
                    errorbar(AllXs(:,ii), AllYs(:,ii) + Offset,...
                        AllYErrs(:, ii), 'ro');
                    hold on;
                    plot(xinterp, Offset + AllFitFnHandles{ii}(xinterp),'b');
                    % plot with maximum and minimum amplitudes from
                    % errorbars
                    plot(xinterp, Offset + ...
                        sinusoid1D([fp_modulation(1,ii),...
                        fp_modulation(2,ii) + stderr_modulation(2,ii),...
                        fp_modulation(3,ii), fp_modulation(4,ii)],...
                        xinterp), 'b--');
                    plot(xinterp, Offset + ...
                        sinusoid1D([fp_modulation(1, ii),...
                        fp_modulation(2,ii) - stderr_modulation(2,ii),...
                        fp_modulation(3,ii), fp_modulation(4,ii)],...
                        xinterp), 'b--');
                    grid on;
                    xlabel('Position (Latt Sites)');
                    ylabel('Filling');
                    ylim([Offset - 0.5, 0.5]);
                end
                ax = gca;
                [SortedOffsets, I] = sort(Offsets);
                ax.YTick = SortedOffsets;
                RedTimes = Times;
                ax.YTickLabels = RedTimes(I);
                
                subplot(1, 2, 2)
                Pos = AllXs(:, 1);
                ColNorm = NormCuts(:, 1);
                ColNormUnc = NormCutsUnc(:, 1);
                ColNormUnc = ColNormUnc( ~isnan(ColNorm) );
                ColNorm = ColNorm( ~isnan(ColNorm) );
                ColDens = AllYs(:, 1) + ColNorm;
                ColDensUnc = sqrt(AllYErrs(:,1).^2 - ColNormUnc.^2);
                errorbar(Pos, ColDens, ColDensUnc, 'r.-');
                hold on;
                errorbar(Pos, ColNorm, ColNormUnc, 'b.-');
                grid on;
                ylim([0,1]);
                xlabel('Position (Latt Sites)')
                ylabel('Filling')
                legend({'0us','long time norm'});
                title('Avg Column Density');
                
                Period = 1 / fp_modulation(1,1);
                PeriodUnc = stderr_modulation(1,1) / fp_modulation(1,1)^2;
                ttl = sprintf('%s\nTimes = %0.0f-%0.0fus \n Period = %0.1f +/- %0.1f, Initial Amp = %0.3f +/- %0.3f, Phase = %0.2f +/- %0.2f',...
                    obj.DescriptStr, Times(1), Times(end),...
                    Period, PeriodUnc,...
                    fp_modulation(2,1), stderr_modulation(2,1),...
                    fp_modulation(3,1), stderr_modulation(3,1));
                suptitle(ttl); 
            
            % export results to text files
            if export_results
                
                names_cell = {'xpos(a)','<n>','<n>_unc','<n>_nobg',...
                    '<n>_nobg_unc','<n>_no_long_time_norm',...
                    '<n>_no_long_time_norm_unc'};
                for ii = 1:ntimes
                    array = zeros(size(AllXs,1), 7);
                    array(:,1) = AllXs(:,ii);
                    array(:,2) = AllYs(:,ii);
                    array(:,3) = AllYErrs(:,ii);
                    array(:,4) = AllYs(:,ii) - fp_modulation(4,ii);
                    array(:,5) = AllYErrs(:,ii);
                    array(:,6) = FitCuts(:,ii);
                    array(:,7) = FitCutsUnc(:,ii);
                        
                    fname = sprintf('%s_time=%0.0fus_period=%0.1f.txt',...
                        export_fname_stem, Times(ii), Period);
                     save_data_file(array, names_cell, '\t', '', save_dir, fname);
                    
                    array = zeros(length(xinterp),3);
                    array(:,1) = xinterp;
                    array(:,2) = AllFitFnHandles{ii}(xinterp);
                    array(:,3) = AllFitFnHandles{ii}(xinterp) - fp_modulation(4,ii);
                    
                    names_cell_fit = {'xpos(a)','<n>','<n>_nobg'};
                    fname = sprintf('%s_time=%0.0fus_period=%0.1f.txt',...
                        export_fname_stem, Times(ii), Period);
                     save_data_file(array, names_cell_fit, '\t', '', save_dir, fname);
                end
            end
             
        end
            
        function fig_handle = plotDiffuseModel(obj, times,...
                amplitudes, amplitude_uncs, periods,...
                model_fn, model_fit_params, model_std_errs,...
                curve_fit_params, curve_fit_std_errs, chi_sqr)
            % plot results of fitDiffuseModel
            
             %%% argument checking
            if ~iscell(times)
                times = {times};
            end
            
            if ~iscell(amplitudes)
                amplitudes = {amplitudes};
            end
            
            if ~iscell(amplitude_uncs)
                amplitude_uncs = {amplitude_uncs};
            end
            
            %interpolated times
            min_time = min([0, cellfun(@min, times)]);
            max_time = max([500, cellfun(@max, times)]);
            interp_times = linspace(min_time, max_time, 200);
            
            num_sets = length(amplitudes);
            num_model_params = length(model_fit_params);
            num_extra_params = round( length(curve_fit_params) / num_sets );
            
            % plot diffuse times and fits
            fig_handle = figure;
            ax1 = subplot(1, 2, 1);
            ax2 = subplot(1, 2, 2);
            
            leg = cell(1, num_sets);
            plot_handles = [];
            for ii = 1 : num_sets
                % plot data on left plot, and data normalize to fit
                % amplitude on right plot
                ph = errorbar(ax1, times{ii}, amplitudes{ii}, amplitude_uncs{ii}, 'o');
                norm = curve_fit_params(num_extra_params * (ii-1) + 1);
                errorbar(ax2, times{ii}, amplitudes{ii} / norm, amplitude_uncs{ii} / norm, 'o', 'color', ph.Color);
                hold(ax1, 'on');
                hold(ax2, 'on');
                
                % plot fits
                temp_curve_params = curve_fit_params(num_extra_params * (ii-1) + 1 : num_extra_params * ii);
                plot(ax1, interp_times, model_fn(model_fit_params, temp_curve_params, interp_times), 'color', ph.Color);
                plot(ax2, interp_times, model_fn(model_fit_params, temp_curve_params, interp_times) / norm, 'color', ph.Color);
                
                leg{ii} = sprintf('period = %0.2f', periods(ii));
                plot_handles(ii) = ph;
            end
            legend(ax1, plot_handles, leg);
            
            % labels
            grid(ax1);
            grid(ax2)
            xlabel(ax1, 'Times (us)');
            xlabel(ax2, 'Times (us)');
            ylabel(ax1, 'Sinusoid Amplitudes');
            ylabel(ax2, 'Sinusoid Amplitudes');
            
            % create title
            ttl = '';
            for ii = 1:num_model_params
                ttl = sprintf('%s\n Param %d = %0.4f +/- %0.4f',...
                    ttl, ii, model_fit_params(ii), model_std_errs(ii));
            end
            ttl = sprintf('%s\n chi sqr = %0.2f', ttl, chi_sqr);
            suptitle(ttl);
        end 
        
        function interleaved_v = interleave_vects(obj, v1, v2)
            if numel(v1) ~= numel(v2)
                error('length of vectors to be interleaved must be the same');
            end
            
            v1 = reshape(v1, [1, length(v1)]);
            v2 = reshape(v2, [1, length(v2)]);
            interleaved_v = [v1; v2];
            interleaved_v = reshape(interleaved_v, [1, numel(interleaved_v)]);
        end
        
        function [offset, offset_unc] = estimate_offset(obj,...
                times, amplitudes, amplitudes_unc, periods, gammas, Ds)
            % Estimate offset for amplitude versus time curves
            % TODO: account for over- versus under-damped behavior. Right
            % now only meainginful for underdamped case
            if ~iscell(times)
                times = {times};
            end
            
            if ~iscell(amplitudes)
                amplitudes = {amplitudes};
            end
            
            if ~iscell(amplitudes_unc)
                amplitudes_unc = {amplitudes_unc};
            end
            
            num_sets = length(amplitudes);
            
            offset = zeros(1, num_sets);
            offset_unc = zeros(1, num_sets);
            
            for ii = 1:num_sets
                curr_times = times{ii};
                curr_amplitudes = amplitudes{ii};
                curr_amplitudes_unc = amplitudes_unc{ii};
                period = periods(ii);
                gamma = gammas(ii);
                D = Ds(ii);
                
                omega = sqrt(gamma * D) * (2*pi/period);
                % look at if we are over- or under-damped
                if omega^2 > gamma^2 / 4
                    osc_period = 2 * pi /omega;
                    cutoff_t = 1.5 * osc_period;
                else
                    tau = 1/(D * (2 * pi / period)^2);
                    cutoff_t = 2 * tau;
                end
                
                offset(ii) = mean(curr_amplitudes( curr_times > cutoff_t));
                npts = length(curr_amplitudes(curr_times > cutoff_t));
                offset_unc(ii) = sqrt( sum( curr_amplitudes_unc( curr_times > cutoff_t ).^2 ) ) / npts;
                
                if isnan(offset(ii)) || npts < 4
                    offset_unc(ii) = 0;
                    offset(ii) = 0;
                end
            end
        end
        
        function fighandle = compareDiffuseAnalysis(obj)
            % compare with and without subtracting background
            fighandle = figure;
            errorbar(obj.Times, obj.fp_modulation(2,:), obj.stderr_modulation(2,:), 'o');
            hold on;
            errorbar(obj.Times, obj.fp_modulation_nobgsub(2,:), obj.stderr_modulation_nobgsub(2,:), 'o');
            
            grid on;
            xlabel('Time (us)');
            ylabel('Amplitude');
            legend({'bg subtraction', 'no bg subtraction'});
            title(sprintf('Amplitude Vs. Diffusion Time\n with and without bg subtraction'));
            
        end
        
        function compareDiffuseTimes(obj, diffuse_sets, save_results, save_dir)
            %
            % compareDiffuseTimes(obj,diffuse_sets,save_results,save_dir)
            %
            % Compare DiffuseTimes plots for different sets.
            %
            % diffuse_sets is a stack of DiffuseSet instances to be
            % compared. It is not necessary to include the current object
            % in this stack. It is automatically added.
            
            if ~exist('save_results','var') || isempty(save_results)
                save_results = 0;
            end
            
            if ~exist('save_dir','var') || isempty(save_dir)
                save_dir = pwd;
            end
            
            OscFn = @(P,T) dampedOsc1D([P(3), 0,...
                sqrt( P(2) * P(1) ) * (2*pi/obj.Period),...
                2 / (P(1)), 0, 0], T);
            
            DatHandles = [];
            FitHandles = [];
            
            if size(diffuse_sets,2)>1
                diffuse_sets = cat(2, obj, diffuse_sets);
            else
                diffuse_sets = cat(1, obj, diffuse_sets);
            end
            
            %TODO this function probably doesn't work anymore because of
            %various changes to the class. So T and Y aren't defined
            %correctly.
            fighandle = figure('name', 'Diffuse Set Comparison');
            subplot(1, 2, 1);
            leg = cell([1,length(diffuse_sets)]);
            for ii = 1:length(diffuse_sets)
                ds = diffuse_sets(ii);
                
                T = ds.Times(1:end-1);
                Y = ds.fp_modulation(2, :);
                YErr = ds.stderr_modulation(2, :);
                DatH = errorbar(T, Y, YErr, 'o');
                hold on;
                
                InterpT = linspace(min(T), max(T), 100);
                FitH = plot(InterpT, OscFn(ds.Fpdo, InterpT), 'color', DatH.Color);
                
                DatHandles = cat(1, DatHandles, DatH);
                FitHandles = cat(1, FitHandles, FitH);
                leg{ii} = sprintf('%s_Period=%0.2f', ds.DescriptStr, ds.Period);
            end
            grid on;
            xlabel('Time (us)');
            ylabel('Amp');
            title('No Normalization');
            legend(DatHandles, leg, 'Interpreter', 'none');
            
            subplot(1, 2, 2);
            DatHandles2 = [];
            for ii = 1:length(diffuse_sets)
                ds = diffuse_sets(ii);
                
                T = ds.Times(1:end-1);
                Y = ds.fp_modulation(2, :) / ds.Fpdo(3);
                YErr = ds.stderr_modulation(2 ,:) / ds.Fpdo(3);
                DatH = DatHandles(ii);
                
                DatH2 = plot(T, Y, '.-', 'color', DatH.Color);
                hold on;
                InterpT = linspace(min(T), max(T), 100);
                DatHandles2 = cat(1, DatHandles2, DatH2);
            end
            grid on;
            xlabel('Time (us)');
            ylabel('Amp');
            title('Normalization');
            legend(DatHandles, leg, 'Interpreter', 'none');
            
            if save_results
                fname = sprintf('%s_compare-diffuse-times.fig', obj.DescriptStr);
                fpath = fullfile(save_dir, fname);
                savefig(fighandle, fpath);
            end
            
        end
        
        function compareDiffuseSets(obj, diffuse_sets, save_results, save_dir, fname)
            % 
            % compareDiffuseSets(obj,diffuse_sets,save_results,save_dir, fname)
            %
            % diffuse_sets is a stack of DiffuseSet objects
            %
            % save_results is bool
            %
            % save_dir is the directory to save figure in
            %
            % fname is the file name of figure to save.
            % 

            if size(diffuse_sets, 2) > 1
                diffuse_sets = cat(2, obj, diffuse_sets);
            else
                diffuse_sets = cat(1, obj, diffuse_sets);
            end
            
            if ~exist('save_results', 'var') || isempty(save_results)
                save_results = 0;
            end
            
            if ~exist('save_dir','var') || isempty(save_dir)
                save_dir = pwd;
            end
            
            if ~exist('fname', 'var') || isempty(fname)
                fname = sprintf('compare-diffuse-results_%s.fig', obj.DescriptStr);
            end
            
            %create comparison figure
            fighandle = figure('name','Compare Diffuse Sets');
            NRows = 2;
            NCols = 3;
            
            leg = cell([1,length(diffuse_sets)]);
            subplot(NRows,NCols,1)
            hold on;
            for ii = 1:length(diffuse_sets)
                ds = diffuse_sets(ii);
                if ~isempty(ds.Singles)
                    rpos = ds.Singles.BinAvg;
                    filling = ds.Singles.Occs_AzAvg;
                    fillingUnc = ds.Singles.Occs_AzAvgUnc;
                    errorbar(rpos,filling,fillingUnc);
                    leg{ii} = sprintf('%s_Period=%0.1f',ds.DescriptStr,ds.Period);
                end
            end
            grid on;
            ylim([0,1]);
            xlabel('Pos (Latt Sites)');
            title('Long Time Filling');
            try
                legend(leg);
            catch
            end
            
            subplot(NRows,NCols,2)
            hold on;
            for ii = 1:length(diffuse_sets)
                ds = diffuse_sets(ii);
                if ~isempty(ds.Density)
                    rpos = ds.dfs(1).BinAvg;
                    n = ds.Density;
                    nunc = ds.DensityUnc;
                    errorbar(rpos,n,nunc);
                end
            end
            grid on;
            ylim([0,1]);
            xlabel('Pos (Latt Sites)');
            title('Density');
            
            subplot(NRows,NCols,3)
            hold on;
            for ii =1:length(diffuse_sets)
                ds = diffuse_sets(ii);
                if ~isempty(ds.Singles) && ~isempty(ds.Density)
                    n = ds.Density;
                    nunc = ds.DensityUnc;
                    ns = ds.Singles.Occs_AzAvg;
                    nsunc = ds.Singles.Occs_AzAvgUnc;
                    
                    errorbar(n,ns,nsunc,nsunc,nunc,nunc);
                end
            end
            ylim([0,1]);
            xlim([0,1]);
            xlabel('Density <n>');
            ylabel('<n^s>');
            title('Density Vs. Singles Density')
            grid on;
            
            subplot(NRows,NCols,4)
            hold on;
            for ii = 1:length(diffuse_sets)
                ds = diffuse_sets(ii);
                if ~isempty(ds.Ones)
                    rpos = ds.dfs(1).BinAvg;
                    Ind = ds.Ones.CenterIndex_CorrMatrix;
                    corr1 = squeeze(ds.Ones.Density_Corr_AzAvg(Ind,Ind+1,:));
                    corr1Unc = squeeze(ds.Ones.Density_Corr_AzAvgUnc(Ind,Ind+1,:));
                    if ~isempty(corr1)
                        errorbar(rpos,corr1,corr1Unc);
                    end
                end
            end
            grid on;
            ylim([-0.06,0.02]);
            xlabel('Pos (Latt Sites');
            ylabel('<n_1n_1>_c');
            title('n1 Correlators')
            
            subplot(NRows,NCols,5)
            hold on;
            for ii = 1:length(diffuse_sets)
                ds = diffuse_sets(ii);
                if ~isempty(ds.Threes)
                    rpos = ds.dfs(end).BinAvg;
                    Ind = ds.Threes.CenterIndex_CorrMatrix;
                    corr3 = squeeze(ds.Threes.Density_Corr_AzAvg(Ind,Ind+1,:));
                    corr3Unc = squeeze(ds.Threes.Density_Corr_AzAvgUnc(Ind,Ind+1,:));
                    if ~isempty(corr3)
                        errorbar(rpos,corr3,corr3Unc);
                    end
                end
            end
            grid on;
            ylim([-0.06,0.02]);
            xlabel('Pos (Latt Sites');
            ylabel('<n_3n_3>_c');
            title('n3 Correlators')
            
            subplot(NRows,NCols,6)
            hold on;
            for ii =1:length(diffuse_sets)
                ds = diffuse_sets(ii);
                if ~isempty(ds.Ones) && ~isempty(ds.Density)
                    rpos = ds.dfs(end).BinAvg;
                    n = ds.Density;
                    nunc = ds.DensityUnc;
                    Ind = ds.Ones.CenterIndex_CorrMatrix;
                    corr1 = squeeze(ds.Ones.Density_Corr_AzAvg(Ind,Ind+1,:));
                    corr1Unc = squeeze(ds.Ones.Density_Corr_AzAvgUnc(Ind,Ind+1,:));
                    
                    errorbar(n,corr1,corr1Unc,corr1Unc,nunc,nunc);
                end
            end
            plot([0,1],[0,0],'k--');
            ylim([-0.06,0.02]);
            xlim([0,1]);
            xlabel('Density <n>');
            ylabel('<n_3n_3>_c');
            title('Correlators')
            grid on;
            
            if save_results
                fpath = fullfile(save_dir, fname);
                savefig(fighandle, fpath);
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %saving and loading routines
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
            if ~exist('save_dir','var') || isempty(save_dir)
                save_dir = pwd;
            end
            
            if ~exist('fname','var') || isempty(fname)
                fname = fullfile([sprintf('%s_struct',obj.DescriptStr),'.mat']);
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
                if ~isa(obj.(fields{ii}), 'DataFolder')
                    AsStruct.(fields{ii}) = obj.(fields{ii});
                else
                    if save_all_folders
                        AsStruct.(fields{ii}) = [];
                        for jj = 1:length(obj.(fields{ii}))
                            obj.(fields{ii})(jj).saveStruct(save_dir);
                        end
                    end
                end
            end
            
            fprintf('Saved to %s\n',fpath);
            save(fpath,'AsStruct');
        end
        
        function loadStruct(obj, StructToLoad, save_dir, load_all_folders)
            % loadStruct(obj, StructToLoad, save_dir, load_all_folders)
            %
            % load class from structure.
            % 
            % StructToLoad may either be the actual structure, or it may be
            % the file path to the structure
            %
            % save_dir is the directory the structure is saved in
            %
            % load_all_folders
            
            % if StructToLoad is a string
            if ischar(StructToLoad)
                if ~exist('save_dir', 'var') || isempty(save_dir)
                    [save_dir, ~, ~] = fileparts(StructToLoad);
                end
                a = load(StructToLoad);
                StructToLoad = a.AsStruct;
            end
            
            % now StructToLoad should be our desired structure
            if ~isstruct(StructToLoad)
                error('StructToLoad was not a structure');
            end
            
            if ~exist('save_dir', 'var') || isempty(save_dir)
                save_dir = pwd;
            end
            
            if ~exist('load_all_folders', 'var')
                load_all_folders = 0;
            end
            
            ObjFields = fieldnames(obj);
            for ii = 1:length(ObjFields)
                obj.(ObjFields{ii}) = [];
            end
            
            StructFields = fieldnames(StructToLoad);
            for ii = 1:length(StructFields)
                if ~any(strcmp({'dfs','Ones','Threes','Singles'},StructFields{ii}))
                    %if ~strcmp(StructFields{ii},'dfs') && ~strcmp(StructFields{ii},'Ones') && ~strcmp(StructFields{ii},'Threes')
                    try
                        obj.(StructFields{ii}) = StructToLoad.(StructFields{ii});
                    catch
                        fprintf('field %s was present in loaded struct, but it is not a field of DiffuseSet class. Skipped this field. \n',StructFields{ii});
                    end
                end
            end
            
            if load_all_folders
                DataFolderStack = [];
                for ii = 1:length(obj.dfs_folder_numbers)
                    dfTemp = DataFolder();
                    fname = sprintf('%04d-%02d-%02d-Folder=%03d.mat',...
                        obj.DateCell{1}, obj.DateCell{2}, obj.DateCell{3},...
                        obj.dfs_folder_numbers(ii));
                    fpath = fullfile(save_dir, fname);
                    struct = load(fpath);
                    dfTemp.loadStruct(struct.DataFolderAsStruct);
                    DataFolderStack = cat(1,DataFolderStack,dfTemp);
                end
                obj.dfs = DataFolderStack;
                %
                if ~isempty(obj.dfs_index_ones)
                    obj.Ones = obj.dfs(obj.dfs_index_ones);
                end
                if ~isempty(obj.dfs_index_threes)
                    obj.Threes = obj.dfs(obj.dfs_index_threes);
                end
                if ~isempty(obj.dfs_index_singles)
                    obj.Singles = obj.dfs(obj.dfs_index_singles);
                end
            end
        end
         
    end
end
