classdef DataFolder < handle
    %DataFolder...class object storing a single data set of Fluorescence
    % images.
    %
    % Some general points of philosophy about using this code:
    % 1. Each function should do one thing. So functions shouldn't both
    %    process and save data, or process and display data, etc.
    %
    %    Functions which are of general interest and which don't require
    %    any DataFolder fields should be stored in an external library. In
    %    some cases, the best approach is to have an external library
    %    functional like this, and then a wrapper function in this class to
    %    make it easier to call. E.g. for fitting data to DQMC, it may be
    %    easiest to write a DQMC fitting function which takes
    %    singles-density, density, and correlations as inputs and returns
    %    the temperature. In the class, there can be a wrapper function
    %    which automatically extracts these different densities from the
    %    fields of the class and does the fitting. The two advantages of
    %    this approach are that the library function can be easily reused,
    %    and the wrapper function saves the trouble of remembering where
    %    the different types of data are stored and cobbling it back
    %    together manually each and every time.
    %
    % 2. With the exception of the initialize function, most other
    %    functions shouldn't store data directly in the class. The idea
    %    here is that you create a class instance, and then you look at it.
    %    So all analysis tasks should be done in the initialize function
    %
    % 3. DataFolder classes each represent a folder of analyzed files. If
    %    you have different types of folders (e.g. ones, threes, singles),
    %    then it may be useful to create a different class which stores
    %    many DataFolder objects. It only makes sense to create a new type
    %    of container when the objects in the container are different. So
    %    in generally it is not necessary to createa a third level of
    %    objects to organize these other objects described here. 
    %       
    %    For example, with the diffusion/transport experiment I created a
    %    DiffuseSet class, which is a container to organize all of the
    %    different DataFolder sets needed to determine temperatures, time
    %    series sets, etc. But since all DiffuseSet objects can be analyzed
    %    in the same way, it should not be necessary to create another
    %    class to organize DiffuseSet objects. I initially tried this, but
    %    it seems there is not much advantage. All of the code can be
    %    written in DiffuseSet without much additional headache, and to the
    %    advantage of simplifying the code.
    %
    % 4. Any operations one identical kinds of DataFolders can be written
    %    in the DataFolder class. E.g. comparisons between DataFolders
    %    should be written in this class. These can then be called by a
    %    variety of organizing classes like 3.
    
    properties        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Analysis settings
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         lib_path = fullfile('..', 'lib')
        ImgCropSize = 100;
        CroppedPicStartCoords % may be used to force use of a certain window if CenterStyle is set to 'external'
        CenterStyle = 'All_CoM' % 'Individual_CoM','All_CoM','Fixed'. Determines
        %how individual pictures are cropped. Either the crop window is fixed, 
        %or it is selected for each image so that the center of mass is at the center of the cropped region.
        XStarts_ROI %Defines the coordinate of the upper left pixel in each cropped image, in terms of the coordinates of the uncropped image.
        YStarts_ROI
        
        % exclude pictures settings
        exclude_on_com = 0
        max_com_stdevs_x = 1
        max_com_stdevs_y = 1
        
        exclude_on_skew = 0
        max_skew_x = 1
        max_skew_y = 1
        
        exclude_on_atom_num = 0;
        max_anum_stdevs = 1;
        
        %AzAvg setings
        AzAvgCentering = '2Dfit' %2Dfit or find or external %TODO implement
        AzAvgDistanceStyle = 'mean' %mean, major, or minor
        AzAvgType = 'spatial' %spatial or density or external
%         UseEllipticalContours = 1 %if 0, forces contours to be circular
        Cx_AzAvg
        Cy_AzAvg
        
        %correlator settings
        NumNeighbors = 4; %number of neighbors to consider in correlator
        OnlyCorrelateSitesSameBin = 1 %TODO check implementation!
        UseSingleQuadCorrMat = 0; %TODO implement. Use only single quadrant of corr matrix w all unique values to save space.
        ErrorAnalysis = 'VarOfSampleCovariance' %'UncProp', 'VarOfSampleCovariance', 'Bootstrap'
        
        %file settings
        PathExp = '.*(?<year>\d{4})[\\\/](?<month>\d{2})[\\\/](?<day>\d{2})[\\\/](?<FolderIndex>\d{3})_.*' 
        %TODO is there a nicer way to deal with different path separators??? Maybe shouldn't be doing this with regexp...
        RecFileExp = '(?<year>\d{4})_(?<month>\d{2})_(?<day>\d{2})_(?<FolderIndex>\d{3})_(?<FileIndex>\d{3})_(?<PicNum>\d{3})_reconstr.mat' 
        FlFileExp =  '(?<hour>\d+)_(?<minute>\d{2})_(?<second>\d{2})_fl.fits'
        RecFileFormat = '*_reconstr.mat' %e.g.2017_08_15_059_001_002_reconstr.mat
        FlFileFormat = '*_fl.fits'            

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Dataset properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dataset %{Year,Month,Day,Dataset,FileNum,PictureNum}
        Date %{Year,Month,Day}
        ExcludedPictures = [] %Pictures to ignore in the analysis
        PictureNumberInFolder = []
        DatasetString = '' % want to replace this with identifier, but have some external dependencies...
        identifier = ''
        creation_time
        CenterIndex_CorrMatrix
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Constants
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kb = 1.3806e-23;
        hbar = 1.0546e-34;
        h = 6.6261e-34;
        mLi = 9.9883e-27;
        a = 1064e-9/sqrt(2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Derived Quantites
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Hubbard Params
        Uint
        tx
        ty
        tdiag
        UoverT
        
        %chemical potential
        ChemPot_0
        ChemPot
        ChemPotUnc
        %Trap Parameters
        Omegas
        OmegaMean
        OmegaMeanUnc
        %compressibilities
        Compressibility
        CompressibilityUnc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %general picture information
        nimgs
        ReconstrInfoStack %NImgs stack of structures
        ReconFPaths %Paths to the reconstruction files
        FlFPaths
        PicTimes
        RecErrRates %Error rates in identifying whether an atom is present or not from reconstruction
        %independant variables
        IndependentVariable %varying over images.
        IndependentVariableUnc
        %atom number data
        MeanAtomNum
        AtomNumSD
        AtomNumbers
        %cloud position data in lattice space
        X_CoM_fullpic %Center of mass of the cloud, in uncropped reconstruction picture coordinates
        Y_CoM_fullpic
        XFromCenter %real space x-distance from azavg center
        YFromCenter %real space y-distance from azavg center
        RFromCenter %real space distance from azavg center
        
        %occupation matrices
        Occs_Stack %Ny x Nx x NImgs
        Occs_ImgAvg %Ny x Nx
        Occs_ImgAvgUnc
        Occs_AzAvg %NBins x 1
        Occs_AzAvgUnc %NBins x 1
        Occs_AzAvgStack %NBins x NImgs
        Occs_AzAvgStackUnc %NBins x NImgs
        Occs_Expanded_AzAvg %Nbins x Ncorr x Ncorr
        Occs_Expanded_AllAvg %NBins x Ncorr x Ncorr...as azimuthal avg, no intermediate image averaging step...
        Occs_Expanded_AllAvgUnc
        %shifted occupations       
        Occs_Shifted_AzAvg %Nbins x Ncorr x Ncorr %only these ones are interesting...because only after doing azimuthal averaging does it become impossible to reconstruct one of these from a single one of the others.
        Occs_Shifted_AzAvgUnc
        Occs_Shifted_AllAvg %
        Occs_Shifted_AllAvgUnc
        
        %correlators = <n_i n_j>, i.e. without subtracting average
        Corr_ImgAvg %Ny x Nx x Ncorr x Ncorr
        Corr_AzAvg %NBins x Ncorr x Ncorr
        Corr_AzAvgUnc
        Corr_AllAvg %NBins x Ncorr x Ncorr
        Corr_AllAvgUnc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %final correlator...i.e. after subtracting average
        %there are three different ways of doing this (at least!) in terms
        %of error analysis...hence all of the different variables.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Method 1. in this case the AzAvg is not being calculated from the image
        %average...but from azimuthal average of individual quantities...
        %This is the main/primary version of the analysis
        Density_Corr_AzAvg_NotSymmetrized
        Density_Corr_AzAvgUnc_NotSymmetrized
        Density_Corr_AzAvg %Ncorr x Ncorr x Nbins
        Density_Corr_AzAvgUnc
        nnc2D
        nnc2Dunc
        %can also calculate g2 function
        g2_FullCorr_AzAvg
        g2_FullCorr_AzAvgUnc
        
        %Method 2. Average images to get <SzSz> and <Sz>. Then compute C_z
        %for each pixel in the image. I.e. get C_z as averaged over images.
        %Finally azimuthally average this C_z.
        FullCorr_ImgAvg
        FullCorr_ImgAvg_AzAvg
        FullCorr_ImgAvg_AzAvgUnc
        
        %Method 3. finally, can compute by treating images and azimuthal regions on
        %the same footing, and averaging all at once
        FullCorr_AllAvg
        FullCorr_AllAvgUnc
        
        %also interested in momentum space.
        Static_Structure_Factor
        Static_Structure_FactorUnc
        Static_Structure_Factor_BootstrapUnc
        Qxs
        Qys
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Error analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Higher order correlators, for error analysis
        Delta2 %E[(n_i-<n_i>)^2 * (n_j-<n_j>)^2]
        Sqrt_Var_SampleCovariance
        ErrorProp_Unc %naive error propogation
        Bootstrap_Density_Corrs
        Bootstrap_Unc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %azimuthal average binning parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        BinAvg %average 'distance' of all points in bin.
        BinUnc %SDM of points in bin
        BinEdges %'distance' endpoints of bin
        NBins %number of bins
        NPtsInBin %points in each bin
        NPtsInBin_Times_NPics
        %average quantities in bin...
        BinOcc %average occupation in a given bin.
        RadialPos %average distance from 2D gaussian fit center in real space.
        RadialPosUnc
        %distance grid/bin information
        DistGrid %TODO: rename AzAvgGrid
        MaskStack
        
        %2D fit results
        GaussFit_ImgAvg %TODO get rid of this too...
        GaussFitParams
        GaussFitUnc
        CloudAspectRatio %TODO get rid of this...
        
        %Miscellaneous
        ColorMap
        
        %fitting data
        NonIntFG = []
    end
    
    %methods intended to be used externally
    methods (Access = public)
        function obj = DataFolder(Dataset, pics_to_exclude, bin_edges, AzAvgType, DistGrid)
            %DataFolder Class for analysis of a single folder of
            %fluorescence images. Can handle a variety of different
            %situations. The two main ones are 1) a folder contains
            %identical shots to be averaged or 2) a folder contains a
            %series of shots where some independent variable is changed
            %versus shot.
            %
            %   obj = DataFolder(Dataset,PicturesToExclude,BinEdges,AzAvgType,DistGrid)
            %
            %%% Arguments:
            %%% Dataset: {Year, Month, Day, Dataset, FileNumber, PictureNumber}
            % or an n x 6 cell array where each row is a different datasets
            %
            %%% PicturesToExclude is a list of the file numbers within this
            %   dataset to exclude. So to exclude the 3rd and 5th files, e.g.
            %   PicturesToExclude = [3,5]
            %
            %%% BinEdges is a list of the bin edges to use in the azimuthal
            %   average
            %
            %%% AzAvgType = 'spatial' or 'density' or 'external'
            %
            %%% DistGrid: If AzAvgType is 'external', use this DistGrid for
            %   the azimuthal average instead of one determined by Gaussian fitting. 
            
            %if you supply arguments, run analysis. If not, instantiate an
            %empty member of the class.
            if ~exist('Dataset', 'var') || isempty(Dataset)
                return;
            end
            
%             if exist('Dataset', 'var') && ~isempty(Dataset)
                
            % number of datasets to process
            nsets = size(Dataset, 1);

            %argument checking.
            if ~exist('pics_to_exclude', 'var') || isempty(pics_to_exclude)
%                     pics_to_exclude = {[]};
                pics_to_exclude = cell(1, nsets);
            end

            if ~iscell(pics_to_exclude)
                pics_to_exclude = {pics_to_exclude};
            end

            if ~exist('bin_edges', 'var')
                bin_edges = [];
            end

            if ~exist('AzAvgType', 'var')
                AzAvgType = [];
            end

            if ~exist('DistGrid', 'var')
                DistGrid = [];
            end     

            %function that does the real work
            % if multiple datasets are given, loop over these but fix
            % region of interest and azimuthal average settings to
            % match the first datasets               

            obj = [];
            for ii = 1 : nsets
                df = DataFolder();

                if ii == 1
                    df.initialize(Dataset(ii, :), pics_to_exclude{ii}, bin_edges, AzAvgType, DistGrid);
                else   
                    % Fix ROI
                    df.CroppedPicStartCoords = obj(1).CroppedPicStartCoords;
                    df.CenterStyle = 'Fixed';
                    % Fix azimuthal average grid
                    DistGrid = obj(1).DistGrid;
                    df.initialize(Dataset(ii, :), pics_to_exclude{ii}, bin_edges, 'external', DistGrid);
                end

                obj = vertcat(obj, df); 
            end
                    
%             end
        end
        
        function initialize(obj, Dataset, PicturesToExclude, BinEdges, AzAvgType, DistGrid)
                %initialize(obj,Dataset,PicturesToExclude,BinEdges,AzAvgType,DistGrid)
                %This function actually imports data and initializes the
                %class. Separated this from the main function so I can
                %instantiate a class without doing this.
                
                %argument checking
                if length(Dataset) ~= 6 %&& length(Dataset) ~= 4
                    error('Dataset was not the correct length. Should be {yyyy, mm, dd, fff, III, ppp} or {yyyy, mm, dd, fff}');
                end
                if ~exist('PicturesToExclude','var')
                    PicturesToExclude = [];
                end
                if ~exist('AzAvgType','var') || isempty(AzAvgType)
                    AzAvgType = 'spatial';
                end
                
                if ~exist('BinEdges','var') || isempty(BinEdges)
                    if strcmp(AzAvgType, 'spatial')
                        BinEdges = sqrt(linspace(0, 45^2, 30));
                    elseif strcmp(AzAvgType, 'density')
                        BinEdges = 0 : 0.1 : 1;
                    else
                        error('Unsupported AzAvgType in DataFolder.m');
                    end
                end
                
                if ~exist('DistGrid', 'var')
                    DistGrid = [];
                end     
                
                %begin processing
                obj.creation_time = now;
                obj.DatasetString = sprintf('%04d-%02d-%02d-Folder=%03d',...
                    Dataset{1}, Dataset{2}, Dataset{3}, Dataset{4});
                obj.identifier = obj.DatasetString;
                fprintf('######################################################\n');
                fprintf('Initialized %s\n', obj.identifier);
                fprintf('######################################################\n');
                
                % if already exists, load from file ...
                saved_fname = sprintf('%s.mat', obj.identifier);
                if exist(saved_fname, 'file')
                    fprintf('Found saved folder data in %s\n', saved_fname);
                    % get last updated timestamp for folder
                    [~, folder_path] = findFileGeneral(cell2mat(Dataset(1:5)), '*.mat');
                    files = dir(folder_path);
                    % first file is directory itself
                    last_update_time = datenum(files(1).date);
                    
                    % load
                    obj_temp = DataFolder();
                    obj_temp.loadStruct(saved_fname);
                    if obj_temp.creation_time > last_update_time
                        obj.loadStruct(saved_fname);
                        fprintf('Loaded data\n');
                        return;
                    else
                        fprintf('Saved data last updated at , but folder was updated at so reanalysis required\n');
                    end
                    
                end
                
                %read in data
                obj.Dataset = Dataset;
                obj.Date = {Dataset{1:3}};
                obj.ExcludedPictures = PicturesToExclude;
                obj.CenterIndex_CorrMatrix = obj.NumNeighbors + 1;
                obj.AzAvgType = AzAvgType;
                if strcmp(AzAvgType, 'external')
                    %not used in this case, but should be set to external 
                    % so settings appear consistent
                    obj.AzAvgCentering = 'external';
                    obj.DistGrid = DistGrid;
                    if isempty(DistGrid)
                        error('AzAvgType set to "external", but DistGrid was empty');
                    end
                end
                
                %load pictures
                CropSize = [obj.ImgCropSize, obj.ImgCropSize];               
                [obj.Occs_Stack, obj.XStarts_ROI, obj.YStarts_ROI,...
                    obj.CroppedPicStartCoords,...
                    obj.X_CoM_fullpic, obj.Y_CoM_fullpic,...
                    obj.PictureNumberInFolder, obj.ReconstrInfoStack,...
                    obj.ReconFPaths, obj.FlFPaths] ...
                    = obj.loadPics(obj.Dataset, CropSize, obj.CenterStyle, obj.CroppedPicStartCoords);
               
                fprintf('ROI center style set to %s\n', obj.CenterStyle);
               
               if isnan(obj.XStarts_ROI)
                   fprintf('No analyzed pictures.\n');
                   return;
               end
               
                % remove excluded pictures, which are passed by argument to instance
                use_index_argument = ~ismember(obj.PictureNumberInFolder, obj.ExcludedPictures);
                if ~isempty(obj.ExcludedPictures)
                     fprintf('Excluding picture %d because instructed to by instance argument\n',...
                        obj.PictureNumberInFolder(~use_index_argument));

                    pictures_excluded_but_not_present = ...
                        obj.ExcludedPictures( ~ismember(obj.ExcludedPictures, obj.PictureNumberInFolder) );
                    if ~isempty(pictures_excluded_but_not_present)
                    fprintf('Picture %d should be excluded, but has not been analyzed or does not exist\n',...
                         pictures_excluded_but_not_present);
                    end
                end
                               
               % remove pictures with atom numbers that deviate too much
               if obj.exclude_on_atom_num
                   % only include pictures that were not already excluded
                   atom_nums = squeeze(sum(sum(obj.Occs_Stack, 1), 2));
                   atom_nums_not_excluded = squeeze(sum(sum(obj.Occs_Stack(:, :, use_index_argument), 1), 2));
                   anum_mean = mean(atom_nums_not_excluded);
                   anum_std = std(atom_nums_not_excluded);
                   
                   anum_max = anum_mean + obj.max_anum_stdevs * anum_std;
                   anum_min = anum_mean - obj.max_anum_stdevs * anum_std;
                   
                   use_index_anum = atom_nums < anum_max & atom_nums > anum_min;
                   fprintf('Excluding picture %d because its atom number deviated too much from the mean\n',...
                        obj.PictureNumberInFolder(~use_index_anum));
               end
                               
               if obj.exclude_on_com
                   % only include pictures we have not already excluded
                   % in means and std deviations
                    x_com_mean = mean(obj.X_CoM_fullpic(use_index_argument)); 
                    xcom_std = std(obj.X_CoM_fullpic(use_index_argument));
                    y_com_mean = mean(obj.Y_CoM_fullpic(use_index_argument));
                    ycom_std = std(obj.Y_CoM_fullpic(use_index_argument));

                    x_com_max = x_com_mean + xcom_std * obj.max_com_stdevs_x;
                    x_com_min = x_com_mean - xcom_std * obj.max_com_stdevs_x;
                    y_com_max = y_com_mean + ycom_std * obj.max_com_stdevs_y;
                    y_com_min = y_com_mean - ycom_std * obj.max_com_stdevs_y;

                    use_index_com = ...
                    obj.X_CoM_fullpic < x_com_max & obj.X_CoM_fullpic > x_com_min &...
                    obj.Y_CoM_fullpic < y_com_max & obj.Y_CoM_fullpic > y_com_min;

                    fprintf('Excluding picture %d because its center of mass deviated too much from the mean\n',...
                        obj.PictureNumberInFolder(~use_index_com));
               end

               if obj.exclude_on_skew
                   [mx_two, my_two, ~] = ...
                       get_moment(obj.Occs_Stack(:, :, use_index_argument), 2, [], []);
                   [mx_three, my_three, ~] = ...
                        get_moment(obj.Occs_Stack(:, :, use_index_argument), 3, [], []);
                   skewness_x = mx_three ./ sqrt(mx_two) .^3;
                   skewness_y = my_three ./ sqrt(my_two) .^3;

                   use_index_skew = ...
                   abs(skewness_x) < obj.max_skew_x & ...
                   abs(skewness_y) < obj.max_skew_y;

                   fprintf('Excluding picture %d because its skew was too large\n',...
                       obj.PictureNumberInFolder(~use_index_skew));
               end

               indices_to_use = use_index_argument; 
               if obj.exclude_on_atom_num
                   indices_to_use = indices_to_use .* use_index_anum;
               end
               if obj.exclude_on_com
                   indices_to_use = indices_to_use .* use_index_com;
               end
               if obj.exclude_on_skew
                   indices_to_use = indices_to_use .* use_index_skew;
               end
               indices_to_use = logical(indices_to_use);
               
                obj.Occs_Stack = obj.Occs_Stack(:, :, indices_to_use);
                obj.XStarts_ROI = obj.XStarts_ROI(indices_to_use);
                obj.YStarts_ROI = obj.YStarts_ROI(indices_to_use);
                obj.PictureNumberInFolder = obj.PictureNumberInFolder(indices_to_use);
                obj.X_CoM_fullpic = obj.X_CoM_fullpic(indices_to_use);
                obj.Y_CoM_fullpic = obj.Y_CoM_fullpic(indices_to_use);
                try
                    obj.ReconstrInfoStack = obj.ReconstrInfoStack(indices_to_use);
                catch
                end
                obj.ReconFPaths = obj.ReconFPaths(indices_to_use);
                obj.FlFPaths = obj.FlFPaths(indices_to_use);
                if numel(obj.IndependentVariable) > 0
                    obj.IndependentVariable = obj.IndependentVariable(indices_to_use);
                end
                    
                    
                %some info about reconstructions
                try
                    Errs1 = transpose([obj.ReconstrInfoStack.errorType1]);
                    Errs2 = transpose([obj.ReconstrInfoStack.errorType2]);
                    obj.RecErrRates = [Errs1,Errs2];
                catch Err
                    disp(Err);
                end
                
                % store number of images
                obj.nimgs = size(obj.Occs_Stack, 3);
                %keep track of some basic data about the atom number
                obj.AtomNumbers = squeeze(sum(sum(obj.Occs_Stack, 1), 2));
                obj.MeanAtomNum = mean(obj.AtomNumbers);
                obj.AtomNumSD = std(obj.AtomNumbers);
                %and binning
                obj.BinEdges = BinEdges;
                obj.NBins = length(obj.BinEdges) - 1;     
                
                %2D fit to average image
                [obj.GaussFitParams, obj.GaussFitUnc, obj.CloudAspectRatio, obj.GaussFit_ImgAvg] = ...
                    obj.get2DFit(mean(obj.Occs_Stack, 3));
                cx = obj.GaussFitParams(1);
                cy = obj.GaussFitParams(2);
                sx = obj.GaussFitParams(3); 
                sy = obj.GaussFitParams(4);
                Theta = obj.GaussFitParams(6);
                AspectRatio = sx / sy;
                obj.RFromCenter = ellipticalGrid(...
                    [obj.ImgCropSize, obj.ImgCropSize],...
                    [cx, cy, AspectRatio, Theta], 'mean');
                
                %azimuthal average grid
                if strcmp(obj.AzAvgCentering, '2Dfit')
                    obj.Cx_AzAvg = obj.GaussFitParams(1);
                    obj.Cy_AzAvg = obj.GaussFitParams(2);
                elseif strcmp(obj.AzAvgCentering,'find')
                    %TODO: implement
                    error('AzAvgCentering = "find" not implemented');
                elseif strcmp(obj.AzAvgCentering, 'external')
                    %in this case, obj.Cx_AzAvg should already have been
                    %set.
                    if strcmp(obj.AzAvgType, 'external')
                        %in this case, get an estimate of the center from
                        %the supplied DistGrid
                        [xx, yy] = meshgrid(1:obj.ImgCropSize, 1:obj.ImgCropSize);
                        [~, min_index] = min(obj.DistGrid(:));
                        obj.Cx_AzAvg = xx(min_index);
                        obj.Cy_AzAvg = yy(min_index);
                    end
                    if isempty(obj.Cx_AzAvg) || isempty(obj.Cy_AzAvg)
                        error('AzAvgCentering = "external", but Cx_AzAvg or Cy_AzAvg was empty');
                    end
                else
                    error('Error, AzAvgCentering not supported');
                end
                
                % TODO; don't like how these variables are named.
                % XFromCenter, YFromCenter, RFromCenter should be
                % calculated from 2D fit, instead of azAvg.
%                 [obj.DistGrid, obj.XFromCenter, obj.YFromCenter, obj.RFromCenter]...
                [obj.DistGrid, ~, ~, ~] = obj.getDistGrid(obj.AzAvgType,...
                    obj.Cx_AzAvg, obj.Cy_AzAvg,...
                    [size(obj.Occs_Stack,1), size(obj.Occs_Stack,2)] );
                fprintf('AzAvgType = %s\n', obj.AzAvgType);
                fprintf('AzAvgCentering = %s\n', obj.AzAvgCentering);
                
                %density and correlator analysis
                n_stack = obj.Occs_Stack;
                AzAvgGrid = obj.DistGrid;
                BinEndPts = obj.BinEdges;
            
                %compute everything we are going to need later.
                % density image and azimuthal averages
                [n, nsdm, r, rsdm, nptsbin, nI, nIsdm, nA, nAsdm, nIsdmA, nAsdmI] = ...
                    obj.getDensMoments(n_stack, AzAvgGrid, BinEndPts);
                
                % compute correlator azimuthal average components
                [nn, nnsdm, ni, nisdm, nj, njsdm, npts] = ...
                    obj.getCorrMoments(n_stack, AzAvgGrid, BinEndPts, obj.NumNeighbors, obj.OnlyCorrelateSitesSameBin);
                
                % compute correlator azimuthal average with error
                [nnc, nncunc] = obj.getCorrWithErr(nn, ni, nj, npts);
                
                % compute full g2 function with error
                [g2, g2unc] = obj.getG2WithErr(nn, ni, nj, npts);
                
                % compute 2D correlator components
                [nnI, nnIsdm, niI, niIsdm, njI, njIsdm, npts] = ...
                    obj.getCorrMomentsImgAvg(n_stack, AzAvgGrid, BinEndPts, obj.NumNeighbors, obj.OnlyCorrelateSitesSameBin);
                
                % compute full 2D correlator
                [obj.nnc2D, obj.nnc2Dunc] = obj.getCorrWithErr(nnI, niI, njI, npts);
                
                % ???
                [nnIcA, nnIcAsdm] = obj.getSupplCorr(nnI,njI,AzAvgGrid,BinEndPts,obj.OnlyCorrelateSitesSameBin);
                
                % compute static structure factor from correlator
                [obj.Qxs, obj.Qys, obj.Static_Structure_Factor, obj.Static_Structure_FactorUnc] = ...
                    obj.getStructureFactor( 0, permute(nnc, [2,3,1]), permute(nncunc, [2,3,1]) );
                
                % compute distance based on 2D fit
                [~, ~, obj.RadialPos, obj.RadialPosUnc, ~, ~, obj.MaskStack] = ...
                    azAvg_General(obj.RFromCenter, [], AzAvgGrid, BinEndPts);

                %assign values to class
                obj.Occs_Shifted_AzAvg = nj; %AzAvg(:,:,:,1);
                obj.Occs_Shifted_AzAvgUnc = njsdm; %AzAvgUnc(:,:,:,1);
                obj.Occs_Expanded_AzAvg = ni; %Occs_Expanded_AzAvg;
                obj.Corr_AzAvg = nn; %AzAvg(:,:,:,2);
                obj.Corr_AzAvgUnc = nnsdm; %AzAvgUnc(:,:,:,2);  
    %             obj.FullCorr_ImgAvg = FullCorr_ImgAvg;
                obj.FullCorr_ImgAvg_AzAvg = nnIcA;
                obj.FullCorr_ImgAvg_AzAvgUnc = nnIcAsdm;
                obj.Occs_ImgAvg = nI; %Occs_ImgAvg;
                obj.Occs_ImgAvgUnc = nIsdm;
                obj.Corr_ImgAvg = nnI; %Corr_ImgAvg;
                obj.Occs_AzAvgStack = nA; %Occs_AzAvgStack;
                obj.Occs_AzAvgStackUnc = nAsdm; %Occs_AzAvgStackUnc;
                obj.Occs_AzAvg = n; %Occs_AzAvg;
                obj.Occs_AzAvgUnc = nsdm; %Occs_AzAvgUnc;
                obj.BinAvg = r; %BinAvg;
                obj.BinUnc = rsdm; %BinUnc;
                obj.NPtsInBin = nptsbin;
                obj.ErrorProp_Unc = sqrt(nnsdm.^2 + (ni.*njsdm).^2 + (nj.*nisdm).^2);
                obj.Sqrt_Var_SampleCovariance = squeeze(permute(nnc, [2,3,1]));

                obj.Density_Corr_AzAvg_NotSymmetrized = squeeze(permute(nnc, [2,3,1]));
                obj.Density_Corr_AzAvgUnc_NotSymmetrized = squeeze(permute(nncunc, [2,3,1]));
                [obj.Density_Corr_AzAvg, obj.Density_Corr_AzAvgUnc] ...
                     = obj.getSymmetricCorrMat(obj.Density_Corr_AzAvg_NotSymmetrized, obj.Density_Corr_AzAvgUnc_NotSymmetrized);

                obj.g2_FullCorr_AzAvg = permute(squeeze(g2), [2, 3, 1]);
                obj.g2_FullCorr_AzAvgUnc = permute(squeeze(g2unc), [2, 3 ,1]);

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Basic analysis functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        function [OccStack,...
                roi_starts_x, roi_starts_y, roi_start_coords,...
                x_com,y_com, pic_num_infolder,...
                rec_info_stack, rec_fpaths, fl_fpaths] = ...
                loadPics(obj, data_set, crop_size, center_style, roi_start_coords)
            %Return a stack of binary
            %reconstruction images from a given folder.
            %  
            %   [OccStack,...
            %   roi_starts_x,roi_starts_y,roi_start_coords,...
            %   x_com,y_com, PicNumInFolder,...
            %   ReconstInfoStack,RecFilePaths,FlFilePaths] = ...
            %   loadPics(obj,Dataset,CropSize,CenterStyle,roi_start_coords)
            %   
            %   Returns OccStack, a 3D array where each slice along the
            %   third dimension is a binary reconstructed image. The
            %   
            %   dataset is specified by Dataset, which is a 1x6 cell array
            %   of numbers = {year month day data-set-number
            %   file-number-in-folder picture-number-in-fil}.
            %
            %   The size of
            %   each 2D image is specified by CropSize, which is a 1x2
            %   vector. 
            %
            %   CenterStyle is a string specifying how to crop the
            %   full reconstructed images to the final size. It can be
            %   either 'Fixed','Individual_CoM', or 'All_CoM'.
            %   'Individual_COM' centers each picture individually about its
            %   center of mass and different images in the stack will have
            %   different coordinates in absolute space. 'All_CoM' centers
            %   all the pictures about a common center of mass and
            %   maintains the absolute coordinates of the pictures (at
            %   least up to how well the reconstruction does this). 'Fixed'
            %   uses the region of interest with upper-most left coordinate
            %   specified by obj.CroppedPicStartCoords. 
            %
            %   The start point of
            %   the cropped region of interest in terms of the coordinates
            %   of the initial reconstruction image is given in XStarts_ROI
            %   and YStarts_ROI. 
            %
            %   Also return a stack of structures produced
            %   by the reconstruction code giving information about the
            %   reconstruction, ReconstInfoStack. 
            %
            %   RecFilePaths is a cell
            %   array giving the full path to each reconstruction file.
            %   
            %   
            %   FlFilePaths is a cell array giving the full path to the
            %   corresponding fluorescence image files which were used to
            %   generate the reconstructed file.
            %
            %   TODO: rewrite to limit use of class object. e.g. get rid of
            %   obj.CroppedPicStartCoords and make an argument.
            %   TODO: probably messiest function left. I am not very happy
            %   with this one. Cleanup!!!
           
            if ~strcmp(center_style, 'Individual_CoM') && ...
               ~strcmp(center_style, 'All_CoM') && ...
               ~strcmp(center_style, 'Fixed')
                error('Center style should to Individual_CoM, All_CoM, or Fixed');
            end
            
           LoadStartT = tic;
            
           %find the paths to all reconstructed fiels
            try
                Dsets = obj.get_FullDataset(data_set);
                PathList = sort(obj.getReconstructedFilenames(Dsets));
            catch Error
                disp(Error.message);
                PathList = '';
            end
            
            rec_info_stack = [];
            rec_fpaths = PathList;
            fl_fpaths = cell(length(PathList), 1);
                       
            NImgs = length(PathList);
            pic_num_infolder = ones([NImgs, 1]);
            OccStack = zeros([crop_size, NImgs]);
            roi_starts_x = ones([NImgs, 1]);
            roi_starts_y = ones([NImgs, 1]);
            x_com = zeros([NImgs, 1]);
            y_com = zeros([NImgs, 1]);
            
            w = warning ('off', 'all'); 
            
            %for each reconstruction file load information and crop picture
            for ii = 1:numel(PathList)
                [~, Fname, Ext] = fileparts(PathList{ii});
                [Tokens, ~] = regexp(fullfile('', [Fname, Ext]),...
                    obj.RecFileExp, 'tokens', 'match');
                pic_num_infolder(ii) = str2double(Tokens{1}{5});
                
                loadfile = load(PathList{ii}, 'occupationsRounded',...
                    'coordTrafo', 'displayROI', 'errorType1',...
                    'errorType2', 'runinfo');
                try
                       rec_info_stack = cat(1, rec_info_stack, rmfield(loadfile, 'occupationsRounded'));
                catch Error
                    disp(Error.message);
                    fprintf('Error concatenating structures at %s \n', PathList{ii});
                end
                
                try
                    fl_fpaths{ii} = loadfile.runinfo.inputFilename;
                catch
                    fl_fpaths{ii} = '';
                end
                
                % TODO: handle the case where there are no atoms and the
                % reconstruction fails. These files should be loaded to
                % prevent problems with numbering, but centering and etc
                % for these should be excluded.
                OccMatFull = double(loadfile.occupationsRounded > 0);
                if isempty(OccMatFull)
                    OccMatFull = zeros(crop_size);
                end
                
                %crop full reconstructed pictures to a region of interest.
                if ii == 1
                    if ~strcmp(center_style, 'Fixed')
%                         [OccMatFull, roi_start_x_centering, roi_start_y_centering] = obj.centerPics(OccMatFull);
                        [OccMatFull, roi_start_x_centering, roi_start_y_centering] = centerPics(OccMatFull);
%                         [OccMatFull, roi_start_x_resizing, roi_start_y_resizing] = obj.resizePics(OccMatFull, CropSize);
                        [OccMatFull, roi_start_x_resizing, roi_start_y_resizing] = resizePics(OccMatFull, crop_size);
                        roi_starts_x(ii) = roi_start_x_resizing + roi_start_x_centering - 1;
                        roi_starts_y(ii) = roi_start_y_resizing + roi_start_y_centering - 1;
                        roi_start_coords = [roi_start_x_resizing + roi_start_x_centering, roi_start_y_resizing + roi_start_y_centering] - 1;
                    
                    else
                        % if we are using Fixed centering style, fix
                        % centering from roi_start_coords
                        if isempty(roi_start_coords)
                           %TODO fix this centering... 
                            roi_start_coords = [floor((size(OccMatFull, 2) - crop_size(2)) / 2),...
                                floor((size(OccMatFull, 1) - crop_size(1)) / 2)];
                            warning('obj.CroppedPicStartcoords not set. Centered on reconstruction.');
                        end
                        
                        StartX = roi_start_coords(1);
                        StartY = roi_start_coords(2);
%                         OccMatFull = obj.getROI(OccMatFull, StartX, CropSize(2), StartY, CropSize(1));
                        OccMatFull = getROI(OccMatFull, StartX, crop_size(2), StartY, crop_size(1));  
                        roi_starts_x = roi_start_coords(1) * roi_starts_x; %XStarts_ROI*(1-obj.CroppedPicStartCoords(1));
                        roi_starts_y = roi_start_coords(2) * roi_starts_y; %YStarts_ROI*(1-obj.CroppedPicStartCoords(2));
                    end
                                        
                    OccStack(:,:,ii) = OccMatFull;
                else
                    try
                        if ~isempty(OccMatFull) && sum(sum(isnan(OccMatFull))) == 0
                            
                            if strcmp(center_style,'Individual_CoM')
                                %center each individual picture
                                [OccMatFull, roi_start_x_centering, roi_start_y_centering] = centerPics(OccMatFull);
                                [OccMatFull, roi_start_x_resizing, roi_start_y_resizing] = resizePics(OccMatFull, crop_size);
                                roi_starts_x(ii) = roi_start_x_resizing + roi_start_x_centering - 1;
                                roi_starts_y(ii) = roi_start_y_resizing + roi_start_y_centering - 1;
                            elseif strcmp(center_style,'All_CoM')
                                %use common centering
                                OccMatFull = getROI(OccMatFull, roi_starts_x(1), crop_size(2), roi_starts_y(1), crop_size(1));
                                roi_starts_x(ii) = roi_starts_x(1);
                                roi_starts_y(ii) = roi_starts_y(1);                               
                            elseif strcmp(center_style,'Fixed')
                                % fixed centering based on roi_starts
                                OccMatFull = getROI(OccMatFull, StartX, crop_size(2), StartY, crop_size(1));
                                roi_starts_x(ii) = StartX;
                                roi_starts_y(ii) = StartY;
                            end
                            
                            OccStack(:,:,ii) = OccMatFull;
                            
                        else
                            %deal with empty images or other problems by adding a slice of NaNs. We
                            %can remove these later...
                            OccStack(:, :, ii) = NaN(size(OccStack, 1), size(OccStack, 2));
                        end
                    catch Error
                        disp(Error.message);
                        OccStack(:, :, ii) = NaN(size(OccStack, 1), size(OccStack, 2));
                    end
                end
                [Xcom_roi, Ycom_roi, ~] = get_moment(OccMatFull, 1, [], []);
                x_com(ii) = Xcom_roi + roi_starts_x(ii) - 1;
                y_com(ii) = Ycom_roi + roi_starts_y(ii) - 1;
            end
            
            w = warning ('on','all');
            
            %find and remove NaN images.
            NanImages = squeeze(all(all(isnan(OccStack), 1), 2));
            ImgIndicesToRemove = find(NanImages);
            OccStack(:, :, ImgIndicesToRemove) = zeros(size(OccStack, 1), size(OccStack, 2), length(ImgIndicesToRemove));
            if ~isempty(ImgIndicesToRemove)
                fprintf('WARNING: Picture Number %d in %d/%d/%d Folder %03d replaced by zeros because had Nans. (loadPics.m) \n',...
                    ImgIndicesToRemove, data_set{1}, data_set{2}, data_set{3}, data_set{4});
            end
            
            LoadStopT = toc(LoadStartT);
            fprintf('Loading %d pictures took %0.2f s \n', numel(PathList), LoadStopT);
            
            if isempty(OccStack)
                OccStack = zeros(crop_size);
                roi_starts_x = NaN;
                roi_starts_y = NaN;
                roi_start_coords = [NaN, NaN];
                x_com = NaN;
                y_com = NaN;
                pic_num_infolder = NaN;
                rec_info_stack = {''};
                rec_fpaths = {''};
                fl_fpaths = {''};
            end
        end
         
        function getCorr_AllAvg(obj)
            %getCorr_AllAvg
            %
            %   this function assumes you've run getAzAvg() first...
            %   does a azimuthal averaging + image averaging at the same time,
            %   treating these on the same footing. Basically another way to
            %   do the error propogation to test that what we are doing makes
            %   sense...
            %
            %   TODO: deprecate this function

            StartT = tic;
            
            DGrid = repmat(obj.DistGrid,[1,1,size(obj.Occs_Stack,3)]);
                %since don't keep Corr around, need to recreate it...
             %[Corr,~,~] = obj.getFullCorrelatorData(obj.Occs_Stack,obj.NumNeighbors);
            
             if obj.OnlyCorrelateSitesSameBin
%                  TODO implement only in single bin functionality here!!!
                w = warning ('off','all');
                %implement averaging only points in same bin. Do this by
                %taking advantage of the fact azAvg_General will ignore
                %NaNs when it takes averages.
                MaskStack = obj.MaskStack;
                SameBinMat = obj.getCorrPtsInSameBin(MaskStack,obj.NumNeighbors);
                SameBinMat(SameBinMat==0)=nan;
                NImgs = size(obj.Occs_Stack,3);
                [~,~,Corr_AllAvg,Corr_AllAvgUnc,~,~,~]= azAvg_General(obj.getCorrelator(obj.Occs_Stack,obj.NumNeighbors).*permute(repmat(SameBinMat,[1,1,1,1,NImgs]),[1,2,5,3,4]),[],DGrid,obj.BinEdges);
             else
                [~,~,Corr_AllAvg,Corr_AllAvgUnc,~,~,~]= azAvg_General(obj.getCorrelator(obj.Occs_Stack,obj.NumNeighbors),[],DGrid,obj.BinEdges);
            end

            obj.Corr_AllAvg = Corr_AllAvg;
            obj.Corr_AllAvgUnc = Corr_AllAvgUnc;

            %no need to do expansion before averaging...
            [~,~,Occs_AllAvg,Occs_AllAvgUnc,~,obj.NPtsInBin_Times_NPics,~] = azAvg_General(obj.Occs_Stack,[],DGrid,obj.BinEdges);
            Occs_Expanded_AllAvg = repmat(Occs_AllAvg,[1,2*obj.NumNeighbors+1,2*obj.NumNeighbors+1]);
            Occs_Expanded_AllAvgUnc = repmat(Occs_AllAvgUnc,[1,2*obj.NumNeighbors+1,2*obj.NumNeighbors+1]);
            
            %shifted occupations
            Occs_Shifted = obj.getAllShiftedMats(obj.Occs_Stack,obj.NumNeighbors);
            [~,~,Occs_Shifted_AllAvg,Occs_Shifted_AllAvgUnc,~,~,~]= azAvg_General(Occs_Shifted,[],DGrid,obj.BinEdges);
            %probably worth storing these because can't get them any other
            %way...also they are very small after averaging.
            obj.Occs_Shifted_AllAvg = Occs_Shifted_AllAvg;
            obj.Occs_Shifted_AllAvgUnc = Occs_Shifted_AllAvgUnc;
            %no need to do shifting before averaging??? Actually, maybe
            %there is.
%             Occs_Shifted_AllAvg = squeeze(obj.getAllShiftedMats(Occs_AllAvg,obj.NumNeighbors));
%             Occs_Shifted_AllAvgUnc =squeeze(obj.getAllShiftedMats(Occs_AllAvgUnc,obj.NumNeighbors));
            
            %%%correlator
            FullCorr_AllAvg = Corr_AllAvg - Occs_Shifted_AllAvg.*Occs_Expanded_AllAvg;
            FullCorr_AllAvgUnc = sqrt((Corr_AllAvgUnc).^2+...
            (Occs_Expanded_AllAvg.*Occs_Shifted_AllAvgUnc).^2+...
            (Occs_Shifted_AllAvg.*Occs_Expanded_AllAvgUnc).^2);
            %%%permute and symmetrize
            FullCorr_AllAvg = permute(squeeze(FullCorr_AllAvg),[2,3,1]);
            FullCorr_AllAvgUnc = permute(squeeze(FullCorr_AllAvgUnc),[2,3,1]);
            [obj.FullCorr_AllAvg,obj.FullCorr_AllAvgUnc] = obj.getSymmetricCorrMat(FullCorr_AllAvg,FullCorr_AllAvgUnc);  
       
            EndT = toc(StartT);
            fprintf('Took %0.2f to run getCorr_AllAvg\n',EndT);
            w = warning ('on','all');
        end
        
        function [DistGrid, XFromCenter, YFromCenter, RFromCenter] = ...
                getDistGrid(obj, azavg_type, cx_azavg, cy_azavg, grid_size)
            %getDistGrid
            %
            %   DistGrid = getDistGrid(obj,AzAvgType,Size) Given Size, a
            %   1x2 array specifying the size of the output array and
            %   AzAvgType, a string which is either 'spatial' or
            %   'density',return a 2D array of distances to be used in 
            %   azimuthal averaging. If AzAvgType is 'spatial', DistGrid is
            %   a type of scaled distance. If AzAvgType is 'density',
            %   DistGrid is related to the density.
            %
            %   TODO: limit dependence on class functions???
            [xx,yy] = meshgrid(1:grid_size(1), 1:grid_size(2));
            XFromCenter = xx - cx_azavg;
            YFromCenter = yy - cy_azavg;
            RFromCenter = sqrt((XFromCenter).^2 + (YFromCenter).^2);
            
            if strcmp(azavg_type, 'spatial')
                Sx = obj.GaussFitParams(3); 
                Sy = obj.GaussFitParams(4);
                Theta = obj.GaussFitParams(6);
                AspectRatio = Sx / Sy;
%                 if ~obj.UseEllipticalContours
%                     AspectRatio = 1;
%                 end
                DistGrid = ellipticalGrid(grid_size,...
                    [cx_azavg, cy_azavg, AspectRatio, Theta],...
                    obj.AzAvgDistanceStyle);
            elseif strcmp(azavg_type, 'density')
                DistGrid = binImg(mean(obj.Occs_Stack,3), 4, 4, 'Same');
            elseif strcmp(azavg_type, 'external')
                %DistGrid was already set in this case.
                DistGrid = obj.DistGrid;
                if isempty(DistGrid)
                    error('DistGrid option was external, but no grid was supplied');
                end
            else
                error('Unsupported AzAvgType in DataFolder.m')
            end
		end
        
        function [FitParams, StdErrs, aspect_ratio, gaussfit_img] = get2DFit(obj, img)
            %get2DFit
            %
            %   [FitParams,StdErrs] = get2DFit(obj,Img) Given a 2D array
            %   Img, fit a 2D gaussian and return the fit parameters,
            %   FitParams and their standard errors, StdErrs. FitParams =
            %   [Center-X,
            %   Center-Y,Sigma-X,Sigma-Y,Amplitude,Theta,Background].
            
            StartT = tic;
            nx = size(img, 2);
            ny = size(img, 1);
            [xx, yy] = meshgrid(1:nx, 1:ny);
            %parameter guesses
            GaussInitParams = [nx/2 ,ny/2, 20, 20, 1, 0, 0];
            FixedP = zeros(1, 7);
            lbs = [0, 0, 0, 0, 0, -inf, -inf];
            ubs = [100, 100, inf, inf, inf, inf, inf];
            [FitParams, ~, GaussFnHandle, StdErrs] = fit2D(xx, yy, img, [],...
                'gaussian2D', GaussInitParams, FixedP,...
                lbs, ubs, [], 'fit');
            
            %extract fit parameters
            %ensure sigmas are positive numbers
            FitParams(3) = abs(FitParams(3)); 
            FitParams(4) = abs(FitParams(4));
            Theta = FitParams(6);
            %Ellipticity = 1;
            
            Sx = FitParams(3);
            Sy = FitParams(4);
            if Sx > Sy
                %make sure that the smaller direction is always Sx...
                %helpful for consistancy between data sets
                Sy = FitParams(3);
                SyErr = StdErrs(3);
                Sx = FitParams(4);
                SxErr = StdErrs(4);
                Theta = Theta + pi/2;
                
                FitParams(3) = Sx;
                FitParams(4) = Sy;
                FitParams(6) = Theta;
                StdErrs(3) = SxErr;
                StdErrs(4) = SyErr;
               
            end            
            
            %get image of fit.
            gaussfit_img = GaussFnHandle(xx, yy);
            aspect_ratio = Sy/Sx;
            EndT = toc(StartT);
            fprintf('get2DFit took %0.2fs \n', EndT);
        end
        
        function findAzAvgCenter(obj, Img, Size, GaussFitParams)
            %findAzAvgCenter. Test importance of azimuthal average
            %centering by performing azimuthal average for small shifts of
            %the center and plotting the results
            %
            %findAzAvgCenter(obj)
            CxStart = GaussFitParams(1);
            CyStart = GaussFitParams(2);
            Sx = GaussFitParams(3); 
            Sy = GaussFitParams(4);
            Theta = GaussFitParams(6);
            AspectRatio = Sy/Sx;
%             if ~obj.UseEllipticalContours
%                 AspectRatio = 1;
%             end
%             Size = [obj.ImgCropSize,obj.ImgCropSize];
            
            StepSize = 0.5;
            MaxSteps = 8;
            [XSteps,YSteps] = meshgrid(-MaxSteps/2:1:MaxSteps/2);
            XSteps = StepSize*XSteps;
            YSteps = StepSize*YSteps;
            
            BinEdges = [0,2:40];
            
            AzAvgs = zeros([length(BinEdges) - 1, numel(XSteps)]);
            Uncs = zeros([length(BinEdges) - 1, numel(XSteps)]);

            for ii = 1:numel(XSteps)
                DistGrid = ellipticalGrid(Size,...
                    [CxStart + XSteps(ii), CyStart + YSteps(ii), AspectRatio, Theta]);
                [~, ~, AzAvg, Unc, ~, ~, ~] = azAvg_General(Img, [], DistGrid, BinEdges);
                AzAvgs(:, ii) = AzAvg;
                Uncs(:, ii) = Unc;
            end
            
            [~, I] = min(AzAvgs(1, :));
            
            plot(AzAvgs);
            grid on;
            title(sprintf('Min Shift at dx = %0.1f, dy = %0.1f \n Cx = %0.1f, Cy = %0.1f',...
                XSteps(I), YSteps(I), CxStart + XSteps(I), CyStart + YSteps(I))); 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Functions for computing densities and density correlators
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        function [n, nsdm, r, rsdm ,npts, nI, nIsdm, nA, nAsdm, nIsdmA, nAsdmI] =...
                getDensMoments(obj, n_stack, AzAvgGrid, BinEndPts)
            % GETDENSMOMENTS Compute various average from a stack of images.
            %[n,nsdm,r,rsdm,npts,nI,nIsdm,nA,nAsdm,nIsdmA,nAsdmI] = getDensMoments(obj,n_stack,AzAvgGrid,BinEndPts)
            %
            %Three types of averages, azimuthal <>_A, image <>_I, and simultaneously <>_IA. 
            %n = <<n>_I>_A = <<n>_A>_I = <n>_IA, so only one type of first moment.
            %nsdm = [sqrt(<n^2>-<n>^2)*sqrt(n/n-1)]/sqrt(n)...the first
            %term in brackets is the standard deviation
            %r is the average distance of each bin, i.e. the average value
            %of AzAvgGrid on each bin
            %rsdm is the sdm
            %npts is the number of points in t
            %We also return the intermediate average
            %nI = nI(x,y) := <n>_I
            %nA = nA(b,m) := <n>_A
            %Second moments are all different...but none are very useful.
            %n2 := <n^2>_IA = <n>_IA. This is the most important second moment. Relevant
            %if we really consider all points in a single bin as equivalent.
            %nA2I := <<n>_A^2>_I this moment may be useful if the quantity we are
            %interested in is the average over a bin.
            %nI2A := <<n>_I^2>_A. This moment tells us about how much points vary between shots.
            %TODO:
            %nIsdmA
            %nAsdmI
            
            StartT = tic;
            %average over images
            num_imgs = size(n_stack, 3);
            nI = squeeze(mean(n_stack, 3)); 
            
            % use the fact images can only take the values zero and one,
            % which implies that we can get the mean square value from the
            % mean
            %1/N_m * \sum_m n(i, j, m)^2 = 1/N_m * \sum_m n(i, j, m)
            nIsdm = sqrt( num_imgs / (num_imgs - 1) ) * sqrt(nI - nI .^ 2) / sqrt(num_imgs);
            % equivalent expression for distribution that can take values
            % besides zero and one
            % nIsdm = std( n_stack, [], 3) / sqrt( num_imgs );
            
            %azimuthal average of image average
            [r, rsdm, n, nIsdmA, ~, npts, ~] = azAvg_General(nI, [], AzAvgGrid, BinEndPts);
            
            %typically define sdm without (n-1) factor. But have to correct for lack of (n-1) factor in SD first.  
            n2 = n;
            nsdm = sqrt( (num_imgs * npts) ./ (num_imgs * npts - 1)) .* sqrt(n2 - n.^2) ./ sqrt(num_imgs * npts);
%             nI2A = (nIsdmA.*sqrt(npts)).^2 + n.^2; %2nd moment
            
            %azimuthal average of each image
            [~, ~, nA, nAsdm, ~, ~, ~] = azAvg_General(n_stack, [], AzAvgGrid, BinEndPts);
            nAsdmI = std(nA, [], ndims(nA)) / sqrt(num_imgs);
%             nA2I = mean(nA.^2,ndims(nA));
            EndT = toc(StartT);
            fprintf('getDensMoments took %0.2fs \n', EndT);
        end
        
        function [nn, nnsdm, ni, nisdm, nj, njsdm, npts] =...
                getCorrMoments(obj, n_stack, AzAvgGrid, BinEndPts, NumNeighbors, RestrictSameBin)
            %GETCORRMOMENTS Compute correlation matrices from a stack of images.
            %
            %   [nn,nnsdm,ni,nisdm,nj,njsdm,npts] = 
            %   getCorrMoments(obj, n_stack, AzAvgGrid, BinEndPts, NumNeighbors, RestrictSameBin)
            %
            %   $$nn = nn(b,i_y,i_x) = <n(y,x,m)n(y+i_y,x+i_x,m)>_{AI} = 
            %   \sum_{y,x,m}
            %   n(y,x,m)*n(y+iy,x+ix,m)*Char_b(y,x)*Char_b(y+i_y,x+i_x)$$
            %   , where b runs over bins, m runs over images, and
            %   Char_b is the characteristic function of the bin b. i.e.
            %   Char_b(x,y) = 1 if (x,y) is in bin b, and 0 otherwise.
            %   nnsdm
            %   $$nj = nj(b,iy,ix) = <n(y+iy,x+ix,m)>_IA = \sum_{y,x,m} 
            %   n(y+iy,x+ix,m) * Char_b(y+iy, x+ix)$$
            %   njsdm
            %   $$ni = ni(b,i_y,i_x) = <n(y,x,m)>_IA = \sum_{y,x,m} 
            %   n(y,x,m)*Char_b(y,x)$$ if we are not restricting to the same
            %   bin.
            %   $$ni = ni(b,i_y,i_x) = <n(y-i_y,x-i_x,m)>_{IA} = \sum_{y,x,m}
            %   n(y-i_y,x-i_x,m)*Char_b(y-i_y,x-i_x)$$ if we are 
            %   restricting to the same bin.
            %   nisdm
            %   $$npts = npts(b,i_y,i_x) = \sum_{y,x,m} 
            %   Char_b(y,x)*Char_b(y+i_y,x+i_x)$$
            
            %nn2 = nn2(bin,i,j) := <nn^2>_AI
            %nnI = nnI(x,y,i,j) := <nn>_I
            %nnI2A = nnI2A(bin,i,j) := <<nn>_I^2>_A
            %nnA = nnA(i,j,m) := <nn>_A
            %nnA2I = nnA2I(bin,i,j) := <<nn>_A^2>_I
            
            StartT = tic;
            numimgs = size(n_stack,3);
            
            if RestrictSameBin
                w = warning ('off','all');
                %implement averaging only points in same bin. Do this by
                %taking advantage of the fact azAvg_General will ignore NaNs
                zos = zeros([size(n_stack, 1), size(n_stack, 2)]);
                [~, ~, ~, ~, ~, ~, MaskStack] = azAvg_General(zos, [], AzAvgGrid, BinEndPts);
                %$$SameBinMat(y, x, i_y, i_x) = 
                %\sum_b Char_b(y,x) * Char_b(y-(i_y-NumNeighbors-1),x-(i_x-NumNeighbors-1))$$
                SameBinMat = obj.getCorrPtsInSameBin(MaskStack, NumNeighbors);
                SameBinMat(SameBinMat == 0) = nan;

                %nj
                %$$nj(b,i_y,i_x) = (1/N_b/NImgs)* \sum_{x,y,m}
                % n(y-(i_y - NumNeighbors - 1), x - (i_x - NumNeighbors - 1),m)
                % * Char_b(y,x) * Char_b(y - (i_y - NumNeighbors - 1), x - (i_x - NumNeighbors - 1))$$
                [~, ~, nj, njIsdmA, ~, ~, ~] = ...
                    azAvg_General(obj.getAllShiftedMats(mean(n_stack, 3), NumNeighbors) .* SameBinMat,...
                    [], AzAvgGrid, BinEndPts);
                
                %ni
                %$$ni(b,i_y,i_x) = (1/N_b/NImgs)* \sum_{x,y,m} 
                % n(y-(-i_y - NumNeighbors - 1), x - (-i_x - NumNeighbors - 1),m)
                % * Char_b(y,x) * Char_b(y - (-i_y - NumNeighbors - 1), x - (i_x - NumNeighbors - 1))$$
%                 ni = flip(flip(nj,ndims(nj)-1),ndims(nj));
                [~, ~, ni, niIsdmA, ~, ~, ~] = ...
                    azAvg_General(repmat(mean(n_stack, 3),...
                    [ones(1, 2), 2 * NumNeighbors + 1, 2 * NumNeighbors + 1])...
                    .* SameBinMat, [], AzAvgGrid, BinEndPts);
                
                %nn
                %$$nn(b,i_y,i_x) = (1/N_b/NImgs)*\sum_{x,y,m} n(y,x,m)*n(y-(i_y -
                %NumNeighbors - 1),x-(i_x - NumNeighbors -
                %1),m)*Char_b(y,x)*Char_b(y-(i_y-NumNeighbors-1),
                %x-(i_x-NumNeighbors-1)),$$ where $N_b = \sum_{x,y} 
                %Char_b(y,x)*Char_b(y-(i_y-NumNeighbors-1),
                [~, ~, nn, nnIsdmA, ~, NPtsInBin, ~] = ...
                    azAvg_General(...
                    squeeze(mean(obj.getCorrelator(n_stack, NumNeighbors), 3)) .* SameBinMat,...
                    [], AzAvgGrid, BinEndPts);
                
                %sdm
                npts = NPtsInBin * numimgs;
                nptsMinusOne = npts;
                nptsMinusOne(nptsMinusOne == 0) = 1;
                njsdm =  sqrt(npts ./ nptsMinusOne) .* sqrt(nj - nj.^2) ./ sqrt(npts);
                nisdm = njsdm;
                nnsdm = sqrt(npts ./ nptsMinusOne) .* sqrt(nn - nn.^2) ./ sqrt(npts);
                
                w = warning('on','all');
            else
                %nj
                [~, ~, nj, njIsdmA, ~, NPtsInBin, ~] = azAvg_General(...
                    obj.getAllShiftedMats(mean(n_stack, 3), NumNeighbors), [], AzAvgGrid, BinEndPts);
                [~, ~, n, ~, ~, npts, ~] = azAvg_General(mean(n_stack, 3), [], AzAvgGrid, BinEndPts);
                %ni
                nsdm = sqrt(npts / (npts - 1)) * sqrt(n - n.^2) / sqrt(npts);
                ni = repmat(n, [ones(1, ndims(n)) ,2 * NumNeighbors + 1, 2 * NumNeighbors + 1]);
                %nn
                [~, ~, nn, nnIsdmA, ~, NPtsInBin, ~] = azAvg_General(...
                    squeeze(mean(obj.getCorrelator(n_stack,NumNeighbors), 3)),...
                    [], AzAvgGrid, BinEndPts);
                
                %sdms
                npts = NPtsInBin * numimgs;
                nptsMinusOne = npts;
                nptsMinusOne(nptsMinusOne == 0) = 1;
                njsdm =  sqrt(npts ./ nptsMinusOne) .* sqrt(nj - nj.^2) ./ sqrt(npts);
                nisdm = repmat(nsdm, [ones(1, ndims(n)) ,2 * NumNeighbors + 1, 2 * NumNeighbors + 1]);
                nnsdm = sqrt(npts ./ nptsMinusOne) .* sqrt(nn - nn.^2) ./ sqrt(npts);
            end
 
            EndT = toc(StartT);
            fprintf('getCorrMoments took %0.2fs \n',EndT);
        end
       
        function [nnI, nnIsdm, niI, niIsdm, njI, njIsdm, npts] = ...
                getCorrMomentsImgAvg(obj, n_stack, AzAvgGrid, BinEndPts, NumNeighbors, RestrictSameBin)
            %getCorrMomentsImgAvg Computes the average 2D correlation matrix over a stack of images.
            %   [nnI,nnIsdm,njI,njIsdm,npts] = getCorrMomentsImgAvg(obj,n_stack,AzAvgGrid,BinEndPts,NumNeighbors,RestrictSameBin)
            %
            %   nnI = nnI(y,x,iy,ix) = <n(x,y)n(x+ix,y+iy)>_I
            %   nnIsdm
            %   niI = niI(y,x,iy,ix) = <n(x,y)>_I
            %   niIsdm
            %   njI = njI(b,iy,ix) = <n(x+ix,y+iy)> Note that there is no point in returning an niI, as...
            %   njIsdm
            %   npts is an array the same size as nnI giving the number of
            %       points averaged to produce each point in nnI.
            
            StartT = tic;
            numimgs = size(n_stack, 3);
            numimgsMinusOne = max([numimgs - 1, 1]);
            nnI = squeeze(mean(obj.getCorrelator(n_stack, NumNeighbors), 3));
            nnIsdm = sqrt(numimgs / numimgsMinusOne) * sqrt(nnI - nnI .^2) / sqrt(numimgs);
             %niI = njI, and there is no point making a distinction until we
                %start doing azimuthal averages.
                
            niI = repmat(mean(n_stack, 3), [1, 1, 2 * NumNeighbors + 1, 2 * NumNeighbors + 1]);
            niIsdm = repmat(std(n_stack, [], 3) / sqrt(size(n_stack, 3)),...
                            [1, 1, 2 * NumNeighbors + 1, 2 * NumNeighbors + 1]);
            njI = squeeze(obj.getAllShiftedMats(mean(n_stack, 3), NumNeighbors));
            njIsdm = sqrt(numimgs / numimgsMinusOne) * sqrt(njI - njI .^2) / sqrt(numimgs);
            
            npts = numimgs * ones(size(nnI));
           
            EndT = toc(StartT);
            fprintf('getCorrMomentsImgAvg took %0.2fs \n', EndT);
        end
        
        function [nnIcA, nnIcAsdm] = getSupplCorr(obj, nnI, njI, AzAvgGrid, BinEndPts, RestrictSameBin)
            %getSupplCorr
            %[nnIcA,nnIcAsdm] = getSupplCorr(obj,nnI,njI,AzAvgGrid,BinEndPts,RestrictSameBin)
            %
            %TODO: finish!
            %nnc := <nn(a,b)>IA-<n(a)>_IA*<n(b)>_IA
            %nnAcI := <<nn(a,b)>_A-<n(a)>_A*<n(b)>_A>_I %this one takes
            %much longer to compute!
            %nnIcA := <<nn(a,b)>_I-<n(a)>_I*<n(b)>_I>_A
            NumNeighbors = (size(nnI,ndims(nnI))-1)/2;
            
            if RestrictSameBin
                w = warning ('off','all');
                %implement averaging only points in same bin. Do this by
                %taking advantage of the fact azAvg_General will ignore NaNs
                [~,~,~,~,~,~,MaskStack] = azAvg_General(zeros([size(nnI,1),size(nnI,2)]),[],AzAvgGrid,BinEndPts);
                SameBinMat = obj.getCorrPtsInSameBin(MaskStack,NumNeighbors);
                SameBinMat(SameBinMat==0)=nan;
                
                niI = flip(flip(njI,ndims(njI)-1),ndims(njI));
                [~,~,nnIcA,nnIcAsdm,~,~,~] = azAvg_General((nnI-niI.*njI).*SameBinMat,[],AzAvgGrid,BinEndPts);
            else
                niI = repmat(njI(:,:,NumNeighbors+1,NumNeighbors+1),[1,1,2*NumNeighbors+1,2*NumNeighbors+1]);
                [~,~,nnIcA,nnIcAsdm,~,~,~] = azAvg_General((nnI-niI.*njI),[],AzAvgGrid,BinEndPts);
            end
            
        end
        
        function [nnc2D, nncunc2D] = getCorrWithErr2D(obj, nnI, njI, npts)
            %Run getCorrMomentsImgAvg to get nnI, njI, and npts from the
            %image stack.
            NumNeighbors = (size(nnI, ndims(nnI)) - 1) / 2;
            niI = repmat(njI(:, :, NumNeighbors + 1, NumNeighbors + 1),...
                [1, 1, size(njI, ndims(njI) - 1), size(njI, ndims(njI))]);
            [nnc2D, nncunc2D] = obj.getCorrWithErr(nnI, niI, njI, npts);
        end
        
        function [nnc, nncunc] = getCorrWithErr(obj, nn, ni, nj, npts)
            %getCorrWithErr Compute density correlator with average value subtracted off and its uncertainty.
            %[nnc,nncunc] = getCorrWithErr(obj,nn,ni,nj,npts)
            %
            %nn = <n(x,y)n(x+i,y+j)>
            %ni = <n(x,y)>
            %nj = <n(x+i,y+j).
            %npts is the same size as nn, ni, and nj, giving the number of
            %points averaged to produce that particular value.
            %nn may be any size. So this function can produce a final g2 if
            %you supply  a fully averaged nn of size NBins x NCorrelators x
            %NCorrelators, or if you supply nn averaged over images of size
            %Nx x Ny x NCorrelators x NCorrelators.

            nnc = (nn - ni .* nj);
            Delta2 = nn - 2*nn.*(ni + nj) + 4*nn.*ni.*nj - 3*ni.^2.*nj.^2 + ni.^2.*nj + ni.*nj.^2;
            nncunc = sqrt((Delta2-nnc.^2)./npts);
            %We are measuring the density correlator Cov(n_i,n_j) = <n_i n_j> - <n_i><n_j>. This
            %correlator is exactly the covariance of the variables n_i. In
            %the language of statistics, we are calculating the Sample
            %Covariance (actually we are only using an estimator of the sample covariance since we
            % are using the population mean instead of the true mean).
            %The uncertainty in this quantity should then be
            %the Variance of the Sample Covariance.
            %One can show Var[Cov(n_i,n_j)] = 1/n*(Delta_2 - Cov(n_i,n_j)^2)
            %where, Delta_2(n_i,n_j) = <(n_i - <n_i>)^2*(n_j - <n_j>^2)> 
            % = <(n_i - 2*n_i<n_i> +<n_i>^2)*(n_j - 2*n_j<n_j> + <n_j>^2)>
            %Now, using n_i^2 = n_i
            % = <n_i n_j> - 2*<n_i n_j>(<n_i>+<n_j>) 
            %+ 4*<n_i n_j><n_i><n_j> -3 <n_i>^2<n_j>^2 
            %+ <n_i>^2<n_j> + <n_i><n_j>^2
            %If we want to correct this expression for the fact we are
            %using the population mean instead of the true mean, we find
            %Var[...] = 1/n*(Delta_2 + 1/(n-1)*Var[n_i]*Var[n_j] - (n-2)/(n-1)*Cov(n_i,n_j))
            %Recall Var[n_i] = <n_i^2> - <n_i>^2.
            %But for large n, this amounts to the same thing as before, so
            %it is not worth pursuing.
            %(see e.g. http://www.randomservices.org/random/sample/Covariance.html)   
        end
        
        function [g2, g2unc] = getG2WithErr(obj, nn, ni, nj, npts)
            %getG2WithErr Compute g2 correlation function with appropriate uncertainty.
            %[g2,g2unc] = getG2WithErr(obj,nn,ni,nj,npts)
            %
            %nn = <n(x,y)n(x+i,y+j)>
            %ni = <n(x,y)>
            %nj = <n(x+i,y+j).
            %npts is the same size as nn, ni, and nj, giving the number of
            %points averaged to produce that particular value.
            %nn may be any size. So this function can produce a final g2 if
            %you supply  a fully averaged nn of size NBins x NCorrelators x
            %NCorrelators, or if you supply nn averaged over images of size
            %Nx x Ny x NCorrelators x NCorrelators.
            
            g2 = nn./(ni.*nj);
            nnsdm = sqrt(npts./(npts-1)).*sqrt(nn-nn.^2)./sqrt(npts);
            nisdm = sqrt(npts./(npts-1)).*sqrt(ni-ni.^2)./sqrt(npts);
            njsdm = sqrt(npts./(npts-1)).*sqrt(nj-nj.^2)./sqrt(npts);
            g2unc = g2.*sqrt((nnsdm./nn).^2 + (nisdm./ni).^2 + (njsdm./nj).^2); 
            %this uncertainty is pretty hard to get. Can't find an expression for it so far. Using naive error propogation instead.
        end
        
        function [qxs, qys, sfact, sfactunc] = getStructureFactor(obj, UseFFT, nnc, nncunc)
            %getStructureFactor Compute static structure factor, which is fourier transform of 
            %correlation matrix.
            %[qxs,qys,sfact,sfactunc] = getStructureFactor(obj,UseFFT,nnc,nncunc)
            %
            %TODO: does fourier transform method give same result as other
            %method?

            if ~exist('UseFFT','var')
                UseFFT = 0;
            end

            %TODO remove these default values from this function
            %TODO get error in case where use FFT
            if ~exist('nnc','var')
                nnc = obj.Density_Corr_AzAvg;
                nncunc = obj.Density_Corr_AzAvgUnc;
            end
            
            if ~exist('nncunc','var')
                nncunc = zeros(size(nnc));
            end
            
            [qxs, qys, sfact, sfactunc] = get_structure_fact(nnc, nncunc);
            
%             NNeighbors = (size(nnc,1)-1)/2;

            %to match with FFT points...use this..
%             if UseFFT
%                 NumQs = size(nnc,1);
%             else
%                 NumQs = size(nnc,1);
%             end
%             [qxs,qys] = meshgrid(0:NumQs-1,0:NumQs-1);
%             qxs = 2*pi/NumQs*qxs;
%             qys = 2*pi/NumQs*qys;
%             qxs = fftshift(qxs);
%             qys = fftshift(qys);
%             %take care of potential numerical error problems
%             qys(abs(qys-pi)<1e-15) = pi; 
%             qxs(abs(qxs-pi)<1e-15) = pi;
%             %finish shifting
%             qys(qys>=pi) = qys(qys>=pi)-2*pi;
%             qxs(qxs>=pi) = qxs(qxs>=pi)-2*pi;
% 
%             sfact = zeros(NumQs,NumQs,obj.NBins);
%             sfactunc = zeros(NumQs,NumQs,obj.NBins);
%             
%             %Can write S(q) = sum_{ix,iy>=0) 2*cos(qx*ix +
%             %qy*iy)<n_i*n_j>_c
%             %=> S(q)_Unc = sqrt(sum_{ix,iy>=0} 2*cos(qx*ix + qy*iy)
%             %sigma_ij^2)
%             if UseFFT
%                 for ii = 1:size(nnc,3)
%                     sfact(:,:,ii) = fftshift(fft2(nnc(:,:,ii)));
%                     sfactunc(:,:,ii) = zeros(size(nnc(:,:,ii)));
%                 end
%             else
%                 Indices = (1:(2*NNeighbors+1))-obj.CenterIndex_CorrMatrix;
%                 [IndicesX,IndicesY] = meshgrid(Indices,Indices);
%                 [IndicesX_Exp,IndicesY_Exp,~] = meshgrid(Indices,Indices,1:obj.NBins);
% 
%                 for ii = 1:NumQs
%                     for jj = 1:NumQs
%                         ExpSingle = exp(1i*qxs(ii,jj)*IndicesX+1i*qys(ii,jj)*IndicesY);
%                         Exp = repmat(ExpSingle,[1,1,obj.NBins]);
%                         sfact(ii,jj,:) = squeeze(sum(sum(nnc.*Exp,2),1));
%     %                     StructFactUncTemp(ii,jj,:) = sqrt(squeeze(sum(sum(obj.Density_Corr_AzAvgUnc.^2,2),1)));
%                         sfactunc(ii,jj,:) = sqrt(squeeze(sum(sum(4*cos(qxs(ii,jj)*IndicesX + qys(ii,jj)*IndicesY).^2.*nncunc.^2.*(IndicesY_Exp>=0),2),1)));
%                     end
%                 end
%             end
            
%             obj.Static_Structure_Factor = sfact; %squeeze(sum(sum(obj.Density_Corr_AzAvg,2),1));
%             obj.Static_Structure_FactorUnc = sfactunc; %sqrt(squeeze(sum(sum(obj.Density_Corr_AzAvgUnc.^2,2),1)));
%             obj.Qxs = qxs;
%             obj.Qys = qys;
        end
        
        %need to move some of these functions to a different section
        function doBootstrapErrorAnalysis(obj, NSets, NTrials)
            %Divide images into NSets groups. From this set, choose at random
            %NTrials groups (allowing us to choose the same group multiple times).
            %Then from these, compute correlator. Repeat this many times
            %and look at the distribution of results.
            StartT = tic;
            
            if ~exist('NSets','var')
                NSets = 30;
            end
            
            if ~exist('NTrials','var')
                NTrials = 1000;
            end
            
            %set up groups of images.
            NImgs = size(obj.Occs_Stack,3);
            EdgeImages = zeros(NSets+1,1);
            EdgeImages(end) = NImgs;
            for ii = 2:(NSets)
                if floor(NImgs/NSets)==0
                    error('In DataFolder doBootstrapErrorAnalysis, chose too may sets for bootstrap');
                end
                EdgeImages(ii) = floor(NImgs/NSets)*(ii-1);
            end
            
            NImgsInBin = EdgeImages(2:end)-EdgeImages(1:end-1);
            RelativeWeights = NImgsInBin/(sum(NImgsInBin));
            %RelativeWeights = ones(length(NImgsInBin),1)/NImgs;
            
            BootStrap_Occs_Expanded_AzAvg = []; %NBins x NCorr x NCorr x NBootstrap
            BootStrap_Occs_Shifted_AzAvg = [];
            BootStrap_Occs_Corr_AzAvg = [];
                       
            %generate average values for each set. Then we can choose any
            %combinations of sets at will and compute the average of that
            %combination...
            for ii = 1:NSets
                fprintf('Set %d/%d \n',ii,NSets);
                ICurrentGroup = EdgeImages(ii)+1:EdgeImages(ii+1);
                Occs_Stack = obj.Occs_Stack(:,:,ICurrentGroup);
                [Corr,Occs_Expanded,Occs_Shifted] = obj.getFullCorrelatorData(Occs_Stack,obj.NumNeighbors);
                
                Occs_ImgAvg = squeeze(mean(Occs_Stack,3));
                Occs_Expanded_ImgAvg = repmat(Occs_ImgAvg,[1,1,2*obj.NumNeighbors+1,2*obj.NumNeighbors+1]);
                Occs_Shifted_ImgAvg = squeeze(mean(Occs_Shifted,3));
                Corr_ImgAvg = squeeze(mean(Corr,3));
                Corr = [];
                
                %Azimuthal averaging
                [~,~,Occs_Expanded_AzAvg,~,~,~,~]= azAvg_General(Occs_Expanded_ImgAvg,[],obj.DistGrid,obj.BinEdges);
                [~,~,Occs_Shifted_AzAvg,~,~,~,~]= azAvg_General(Occs_Shifted_ImgAvg,[],obj.DistGrid,obj.BinEdges);
                [~,~,Corr_AzAvg,~,~,~,~]= azAvg_General(Corr_ImgAvg,[],obj.DistGrid,obj.BinEdges);
                          
                BootStrap_Occs_Expanded_AzAvg = cat(4,BootStrap_Occs_Expanded_AzAvg,Occs_Expanded_AzAvg);
                BootStrap_Occs_Shifted_AzAvg = cat(4,BootStrap_Occs_Shifted_AzAvg,Occs_Shifted_AzAvg);
                BootStrap_Occs_Corr_AzAvg = cat(4,BootStrap_Occs_Corr_AzAvg,Corr_AzAvg);
            end
            
            obj.Bootstrap_Density_Corrs = [];
            BootStrap_StructFact = []; 
            for ii = 1:NTrials
                fprintf('Trial %d/%d \n',ii,NTrials);
                Blocks = ceil(rand(NSets,1)*NSets);
                
                %                 Trial_Occs_Expanded_AzAvg = mean(BootStrap_Occs_Expanded_AzAvg(:,:,:,Blocks),4);
                %                 Trial_Occs_Shifted_AzAvg = mean(BootStrap_Occs_Shifted_AzAvg(:,:,:,Blocks),4);
                %                 Trial_Corr_AzAvg = mean(BootStrap_Occs_Corr_AzAvg(:,:,:,Blocks),4);
                CurrentWeights = zeros(1,1,1,length(RelativeWeights));
                CurrentWeights(1,1,1,:) = RelativeWeights(Blocks)/sum(RelativeWeights(Blocks));
                Repeats = size(BootStrap_Occs_Expanded_AzAvg);
                Repeats = Repeats(1:end-1);
                CurrentWeights = repmat(CurrentWeights,[Repeats,1]);
                Trial_Occs_Expanded_AzAvg = sum(BootStrap_Occs_Expanded_AzAvg(:,:,:,Blocks).*CurrentWeights,4);
                Trial_Occs_Shifted_AzAvg = sum(BootStrap_Occs_Shifted_AzAvg(:,:,:,Blocks).*CurrentWeights,4);
                Trial_Corr_AzAvg = sum(BootStrap_Occs_Corr_AzAvg(:,:,:,Blocks).*CurrentWeights,4);
                Trial_Density_Corr_AzAvg = Trial_Corr_AzAvg - Trial_Occs_Expanded_AzAvg.*Trial_Occs_Shifted_AzAvg;
                
                
                obj.Bootstrap_Density_Corrs = cat(4,obj.Bootstrap_Density_Corrs,Trial_Density_Corr_AzAvg);
                
                [~,~,StructFactTemp,~] = obj.getStructureFactor(permute(Trial_Density_Corr_AzAvg,[2,3,1]));
                BootStrap_StructFact = cat(4,BootStrap_StructFact,StructFactTemp);
            end
            obj.Bootstrap_Unc = std(obj.Bootstrap_Density_Corrs,0,4);
            obj.Static_Structure_Factor_BootstrapUnc = std(BootStrap_StructFact,0,4);
            
            EndT = toc(StartT);
            fprintf('took %0.2f s to run bootstrap with %d sets and %d trials \n',EndT,NSets,NTrials);
            
        end
        
        function doBootstrapErrorAnalysis2(obj, NSets, NTrials)
            %TODO replace doBootstrapErrorAnalysis with this function
            %(seems to be working...)
            %TODO get errorbars for other things...like densities and etc.
            %from this method also...
            %Divide images into NSets groups. From this set, choose at random
            %NTrials groups (allowing us to choose the same group multiple times).
            %Then from these, compute correlator. Repeat this many times
            %and look at the distribution of results.
            StartT = tic;
            
            %any class objects should be assigned to local variables here.
            n_stack = obj.Occs_Stack;
            NumNeighbors = obj.NumNeighbors;
            AzAvgGrid = obj.DistGrid;
            BinEndPts = obj.BinEdges;
            RestrictSameBin = obj.OnlyCorrelateSitesSameBin;
            
            if ~exist('NSets','var')
                NSets = 30;
            end
            
            if ~exist('NTrials','var')
                NTrials = 1000;
            end
            
            %set up groups of images.
            NumImgs = size(n_stack,ndims(n_stack));
            EdgeImages = zeros(NSets+1,1);
            EdgeImages(end) = NumImgs;
            for ii = 2:(NSets)
                if floor(NumImgs/NSets)==0
                    error('In DataFolder doBootstrapErrorAnalysis, chose too may sets for bootstrap');
                end
                EdgeImages(ii) = floor(NumImgs/NSets)*(ii-1);
            end
            
            %generate average values for each set. Then we can choose any
            %combinations of sets at will and compute the average of that
            %combination...
            ni_AllSets = []; %NBins x NCorr x NCorr x NBootstrap
            nj_AllSets = [];
            nn_AllSets = [];
            for ii = 1:NSets
                fprintf('Set %d/%d \n',ii,NSets);
                ICurrentGroup = EdgeImages(ii)+1:EdgeImages(ii+1);
                n_stack_set = n_stack(:,:,ICurrentGroup);
                [nn_set,~,ni_set,~,nj_set,~,~] = obj.getCorrMoments(n_stack_set,AzAvgGrid,BinEndPts,NumNeighbors,RestrictSameBin);
          
                ni_AllSets = cat(4,ni_AllSets,ni_set);
                nj_AllSets = cat(4,nj_AllSets,nj_set);
                nn_AllSets = cat(4,nn_AllSets,nn_set);
            end
            
            %setup trials
            NImgsInBin = EdgeImages(2:end)-EdgeImages(1:end-1);
            RelativeWeights = NImgsInBin/(sum(NImgsInBin));
            nnc_bootstrap = [];
            sfact_bootstrap = []; 
            for ii = 1:NTrials
                fprintf('Trial %d/%d \n',ii,NTrials);
                Blocks = ceil(rand(NSets,1)*NSets);
                CurrentWeights = zeros(1,1,1,length(RelativeWeights));
                CurrentWeights(1,1,1,:) = RelativeWeights(Blocks)/sum(RelativeWeights(Blocks));
                Repeats = size(ni_AllSets);
                Repeats = Repeats(1:end-1);
                CurrentWeights = repmat(CurrentWeights,[Repeats,1]);
                
                ni_trial = sum(ni_AllSets(:,:,:,Blocks).*CurrentWeights,4);
                nj_trial = sum(nj_AllSets(:,:,:,Blocks).*CurrentWeights,4);
                nn_trial = sum(nn_AllSets(:,:,:,Blocks).*CurrentWeights,4);
                [nnc_trial,~] = obj.getCorrWithErr(nn_trial,ni_trial,nj_trial,ones(size(ni_trial)));
                
                nnc_bootstrap = cat(4,nnc_bootstrap,nnc_trial);
                
                [~,~,sfact_temp,~] = obj.getStructureFactor(permute(nnc_trial,[2,3,1]));
                sfact_bootstrap = cat(4,sfact_bootstrap,sfact_temp);
            end
            
            %everything passed to class members should be done after this
            %point.
            obj.Bootstrap_Density_Corrs = nnc_bootstrap;
            obj.Bootstrap_Unc = std(nnc_bootstrap,0,4);
            obj.Static_Structure_Factor_BootstrapUnc = std(sfact_bootstrap,0,4);
            
            EndT = toc(StartT);
            fprintf('took %0.2f s to run bootstrap with %d sets and %d trials \n',EndT,NSets,NTrials);     
        end
        
        function [NN_Errs] = testBootstrap(obj, MaxSets)
            %Test function for running bootstrap
            %MaxSets = 10;
            NN_Errs = zeros(MaxSets-1,1);
            for ii = 2:MaxSets
                obj.doBootstrapErrorAnalysis(ii,1000);
                NN_Errs(ii-1) = obj.Bootstrap_Unc(1,5,6);
            end
            plot([2:1:MaxSets],NN_Errs,'bo');
            xlabel('Number of sets divided images into')
            ylabel('Unc')
        end
        
        function ShiftMatStack = getQuadShiftedMats(obj, MatStack, NumNeighbors)
            %TODO replace ??? with what?
            [XShift,YShift] = meshgrid(0:1:NumNeighbors);
            ShiftMatStack = zeros([size(MatStack),size(XShift)]);
            if ndims(MatStack)==3
                for ii = 1:numel(XShift)
%                     [ShiftMatStack(:,:,:,ii),~,~] = obj.getShiftedMat(MatStack,XShift(ii),YShift(ii),0);
                    [ShiftMatStack(:, :, :, ii), ~, ~] = getShiftedMat(MatStack, XShift(ii), YShift(ii), 0);
                end
            elseif ndims(MatStack)==2
                for ii = 1:numel(XShift)
%                 	ShiftMatStack(:,:,ii) = obj.getShiftedMat(MatStack,XShift(ii),YShift(ii),0);
                    ShiftMatStack(:,:,ii) = getShiftedMat(MatStack, XShift(ii), YShift(ii), 0);
                end
            else
                disp('Congratulations, the n-d array size you used is not supported by getAllShiftedMats. Maybe write this function in a more clever way so it works with all sizes. Thank you!');
            end
        end
        
        function ShiftMatStack = getAllShiftedMats(obj, MatStack, NumNeighbors)
            %getAllShiftedMats Returns shifted versions of a matrix along extra dimensions.
            %
            %   ShiftMatStack = getAllShiftedMats(obj,MatStack,NumNeighbors)
            %   If MatStack is an nd array with dimensions n1 x ... x nm and
            %   NumNeighbors is an integer, then ShiftMatStack is the n1 x ...
            %   x nm x 2*NumNeighbors + 1 x 2*NumNeighbors + 1 array where
            %   each extra dimensions is the original nd array shifted in its
            %   first two coordinates by some amount. So,
            %   ShiftMatStack(x1,x2,x3,...xm,i1,i2) = 
            %   MatStack(x1 - (1i - NumNeighbors -1) , x2 - (i2 - NumNeighbors - 1),x3,...,xm).
            %   
            %Currently this function is only implemented for m = 2 and 3.
            [XShift,YShift] = meshgrid(-NumNeighbors:1:NumNeighbors);
            ShiftMatStack = zeros([size(MatStack),size(XShift)]);
            if ndims(MatStack)==3
                for ii = 1:numel(XShift)
%                     ShiftMatStack(:,:,:,ii) = obj.getShiftedMat(MatStack,XShift(ii),YShift(ii),0);
                    ShiftMatStack(:, :, :, ii) = getShiftedMat(MatStack, XShift(ii), YShift(ii), 0);
                end
            elseif ndims(MatStack)==2
                for ii = 1:numel(XShift)
%                 	ShiftMatStack(:,:,ii) = obj.getShiftedMat(MatStack,XShift(ii),YShift(ii),0);
                    ShiftMatStack(:, :, ii) = getShiftedMat(MatStack, XShift(ii), YShift(ii), 0);
                end
            else
                disp('Congratulations, the n-d array size you used is not supported by getAllShiftedMats. Maybe write this function in a more clever way so it works with all sizes. Thank you!');
            end
        end
        
        function [CorrStack] = getCorrelator(obj, ImageStack, NumNeighbors)
            %getCorrelator Compute 2D correlation matrices for every image in a stack of images.
            %
            %   [CorrStack] = getCorrelator(obj,ImageStack,NumNeighbors) If
            %   ImageStack is an Ny x Nx x NImgs array, then CorrStack is
            %   an Ny x Nx x NImgs x 2*NumNeighbors+1 x 2*NumNeighbors+1
            %   array. CorrStack(y,x,m,iy,ix) =  
            %   ImageStack(y,x,m)*ImageStack(y-(iy-NumNeighbors-1),x-(ix-NumNeighbors-1),m)
            %
            %   TODO: use this and previous functions in place of
            %   getFullCorrelatorData
           if ~isnumeric(ImageStack)
                disp('First argument not numeric')
           end
           
           [XShift,YShift] = meshgrid(-NumNeighbors:1:NumNeighbors);
           CorrStack = zeros([size(ImageStack),size(XShift)]);
           if ndims(ImageStack)==3
               for ii = 1:numel(XShift)
%                     CorrStack(:,:,:,ii) = ImageStack.*obj.getShiftedMat(ImageStack,XShift(ii),YShift(ii),0);
                    CorrStack(:, :, :, ii) = ImageStack .* getShiftedMat(ImageStack, XShift(ii), YShift(ii), 0);
               end
           elseif ndims(ImageStack)==2
               for ii = 1:numel(XShift)
%                     CorrStack(:,:,ii) = ImageStack.*obj.getShiftedMat(ImageStack,XShift(ii),YShift(ii),0);
                    CorrStack(:, :, ii) = ImageStack .* getShiftedMat(ImageStack, XShift(ii), YShift(ii), 0);
               end
           else
               error('unsupported ndims for ImageStack');
           end

        end
        
        function [CorrStack, OccStack, ShiftedOccStack] =...
                getFullCorrelatorData(obj, ImageStack, NumNeighbors)
            %%%TODO fully deprecate this function.
            %[SzSz,Occs,ShiftedOccs] = getFullCorrelatorData(ImageStack,NumNeighbors)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %ImageStack has size Ny x Nx x NImgs ... a stack of occupation matrices
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %NumNeighbors is the number of neighbors to consider for the correlation
            %matrix. The size of the correlation matrix will be 2*NumNeighbors+1 x
            %2*NumNeighbors +1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %OccStack is size Ny x Nx x NImgs x 2*NumNeighbors + 1 x 2*NumNeighbors +1
            %...this is simply ImageStack repeated along extra dimensions...
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %SzSzStack and ShiftedOccStack have size Ny x Nx x NImgs x 2*NumNeighbors+1 x 2*NumNeighbors+1
            %Think of them as stacks of correlation matrices...and images...
            %NOTE: this shifts matrices using circshift...possible introduces problems
            %at edges...
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %error checking on inputs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isnumeric(ImageStack)
                disp('First argument not numeric')
            end

            % if length(RadiiBinEdges)~=2
            %     disp('RadList must be a list of length 2')
            % end

            if ~isnumeric(NumNeighbors) || isnan(NumNeighbors) || isinf(NumNeighbors)
                disp('NumNeighbors was not a number')
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %real function begins.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            FullSize = 2*NumNeighbors + 1;

            Nx = size(ImageStack,2);
            Ny = size(ImageStack,1);
            NImgs = size(ImageStack,3);

            %describes shifts
            [ShiftIndexMatrix_X,ShiftIndexMatrix_Y] = meshgrid(-NumNeighbors:1:NumNeighbors);
            
            %CorrMatrixPoints = zeros(NSitesToAvg*NImgs,FullSize,FullSize); %Nx7x7 matrix...
            ShiftedOccStack = zeros(Ny,Nx,NImgs,FullSize,FullSize);
            OccStack = repmat(ImageStack,[1,1,1,FullSize,FullSize]);
            %Occs = repmat(ImageStack,[1,1,1,FullSize,FullSize]); %this is useful later...it is the occupation matrix, but same size as correlation matrix. defined next.
            CorrStack = zeros(Ny,Nx,NImgs,FullSize,FullSize); %thought this is an easier object to deal with...

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %compute correlators
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for ll = 1:NImgs
                OccMat = ImageStack(:,:,ll);
                for jj = 1:FullSize
                    for kk = 1:FullSize
                        %TODO replace circ shifts with
                        %obj.getShiftedMat...check they give the same
                        %results...
                        %think can do circ shift w/o loop over stack...???
                        ShiftedMat = circshift(OccMat, [ShiftIndexMatrix_Y(jj, kk), ShiftIndexMatrix_X(jj, kk)]);
                        CorrStack(:, :, ll ,jj, kk) =  OccMat .* ShiftedMat;
                        %CorrelatorAvg_SzSz(:,:,j,k) = CorrelatorAvg_SzSz(:,:,j,k) + OccMat.*ShiftedMat/StackSize;
                        ShiftedOccStack(:, :, ll, jj, kk) = ShiftedMat;
%                         SzSzStack(:,:,:,jj,kk) = obj.getCorrelator(ImageStack,ShiftIndexMatrix_X(jj,kk),ShiftIndexMatrix_Y(jj,kk));
                    end
                end
            end
        end     
             
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Less important functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [times] = getPicTimes(obj, fl_file_paths)
            %getPicTimes Return times all pictures used in this class
            %instance were taken at.
            %
            %   [Times] = getPicTimes(obj,FlFilePaths) Given FlFilePaths, a
            %   cell array of file paths, return a corresponding array of
            %   the times these pictures were taken.
            %
            %see also datenum
            
            if ischar(fl_file_paths)
                fl_file_paths = {fl_file_paths};
            end
            
            times = zeros(length(fl_file_paths), 1);
            
            for ii = 1:length(fl_file_paths)
                [dir_path, file_name, ext] = fileparts(fl_file_paths{ii}); 
                file_name = fullfile('', [file_name, ext]);
                [tokens, matches] = regexp(dir_path, obj.PathExp, 'tokens', 'match');
                year = str2double(tokens{1}{1});
                month = str2double(tokens{1}{2});
                day = str2double(tokens{1}{3});
                
                [tokens, matches] = regexp(file_name, obj.FlFileExp, 'tokens', 'match');
                hour = str2double(tokens{1}{1});
                minute = str2double(tokens{1}{2});
                second = str2double(tokens{1}{3});
                times(ii) = datenum(year, month, day, hour, minute, second);
            end
        end
        
        function [times, phi1s, phi2s, reordering_index] = getRecPhases(obj, rec_info_stack)
                % Determine lattice phase drift over time
            
                if ~exist('RecInfoStack', 'var') || isempty(rec_info_stack)
                    rec_info_stack = obj.ReconstrInfoStack;
                end
                
                phis = [];
                times = [];

                phi_temp = [rec_info_stack.coordTrafo];
                fl_fnames = {};
                for jj = 1:length(phi_temp)
                    fl_fnames{jj} = rec_info_stack(jj).runinfo.inputFilename;
                    
                    temp_phis = [phi_temp(jj).gridParameters.phi1,...
                                 phi_temp(jj).gridParameters.phi2];
                    
                    phis = cat(1, phis, temp_phis);
                end
                
%                 try
                times = cat(1, times, obj.getPicTimes(fl_fnames));
                times = (times - times(1)) * 24 * 60 * 60; %number is # of days
%                 catch e
%                     disp(e);
%                     times = 1:length(phi_temp);
%                 end
                
                [times, reordering_index] = sort(times);
                %unwrap phases
                phi1s = unwrap(2 * pi * phis(reordering_index, 1)) / (2*pi);
                phi2s = unwrap(2 * pi * phis(reordering_index, 2)) / (2*pi);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Test consistancy of different methods.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function FigHandle = compareAveraging(obj)
            %compareAveraging Compare different types of averaging correlators. i.e. average
            %over pictures then azimuthally vs. average over pictures and
            %azimuthally simultaneously.
            %
            %   FigHandle = compareAveraging(obj) Returns a handle to the
            %   figure which displays a comparison of the different
            %   averaging methods for creating the full correlation matrix.
            
            FigHandle = figure('name','Compare Averaging');
            %%%AzAvg
            FullCorr = permute(obj.Density_Corr_AzAvg(obj.CenterIndex_CorrMatrix,...
                obj.CenterIndex_CorrMatrix + 1, :), [3, 1, 2]);
            FullCorrErr = permute(obj.Density_Corr_AzAvgUnc(obj.CenterIndex_CorrMatrix,...
                obj.CenterIndex_CorrMatrix + 1, :), [3, 1, 2]);
            errorbar(obj.BinAvg, FullCorr, FullCorrErr, 'r-o');
            hold on;
            %%%ImgAvg_AzAvg
            FullCorr = permute(obj.FullCorr_ImgAvg_AzAvg(obj.CenterIndex_CorrMatrix,...
                obj.CenterIndex_CorrMatrix + 1, :), [3, 1, 2]);
            FullCorrErr = permute(obj.FullCorr_ImgAvg_AzAvgUnc(obj.CenterIndex_CorrMatrix,...
                obj.CenterIndex_CorrMatrix + 1, :), [3, 1, 2]);
            errorbar(obj.BinAvg, FullCorr, FullCorrErr, 'b-o');
            %%%AllAvg
            if ~isempty(obj.FullCorr_AllAvg)
                FullCorr = permute(obj.FullCorr_AllAvg(obj.CenterIndex_CorrMatrix,...
                    obj.CenterIndex_CorrMatrix + 1, :), [3 ,1, 2]);
                FullCorrErr = permute(obj.FullCorr_AllAvgUnc(obj.CenterIndex_CorrMatrix,...
                    obj.CenterIndex_CorrMatrix + 1, :), [3, 1, 2]);
                errorbar(obj.BinAvg, FullCorr, FullCorrErr, 'g-o');
            end
            hold off;
            grid on;
            legend({'Azimuthal Average quantities, then compute correlator',...
                'Compute correlator, then azimuthal average',...
                'Average on same footing, then compute correlator'})
            suptitle(obj.identifier);
        end
        
        function FigHandle = compareUncertainty(obj)
            %compareUncertainty Produces a figure comparing different
            %methods for estimating the uncertainty of the density
            %correlator.
            %
            %   FigHandle = compareUncertainty(obj) Compares uncertainties 
            %   for different error propogation methods and returns a
            %   handle to the comparison figure, FigHandle. The three types
            %   of error propogation considered are a bootstrap,naive error
            %   propogation, and proper error propogation. WARNING: 
            %   bootstrap takes quite a while to run!
            
            if isempty(obj.Bootstrap_Unc)
                obj.doBootstrapErrorAnalysis(30,1000);
            end
            FigHandle = figure('name','Compare Uncertainty');
            CenterIndex = obj.CenterIndex_CorrMatrix;
            errorbar(obj.BinAvg,squeeze(obj.Density_Corr_AzAvg(CenterIndex, CenterIndex + 1 ,:)),...
                obj.Sqrt_Var_SampleCovariance(:, CenterIndex, CenterIndex + 1), 'b');
            hold on;
            errorbar(obj.BinAvg,squeeze(obj.Density_Corr_AzAvg(CenterIndex, CenterIndex + 1, :)),...
                squeeze(obj.Bootstrap_Unc(:, CenterIndex, CenterIndex + 1)), 'go');
            errorbar(obj.BinAvg,squeeze(obj.Density_Corr_AzAvg(CenterIndex, CenterIndex + 1, :)),...
                squeeze(obj.ErrorProp_Unc(:, CenterIndex, CenterIndex + 1)), 'r');
            %errorbar(obj.BinAvg,squeeze(obj.Density_Corr_AzAvg(5,6,:)),squeeze(obj.Density_Corr_AzAvgUnc(5,6,:)),'r');
            grid on;
            legend({'Variance of Sample Covariance', 'Bootstrap', 'Error Propogation'});
            suptitle(obj.identifier);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Fit to various profiles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TODO cleanup these functions. Add errobars. Possibly introduce a
        %more streamlined way of doing things. e.g. have a single main
        %function that does most of the fitting stuff, and the individual
        %functions are mostly for interpretation?

        function fitAtomicLimit(obj, UKHz, mode)
            %fitAtomicLimit
            %
            %   [Fp,FFH,SE] = fitAtomicLimit(obj,U,mode) Given U in KHz and
            %   a handle to axes, Axes, returns fit parameters Fp, FFH an
            %   anonymous function handle to the fit function, and SE, the
            %   standard error of the fit. Fp =
            %   [Beta,Omega,AvgMu,DeltaMu,DetectionEfficiency,Sign of U]
            %
            %   mode = 'singles', 'doubles', 'ups', or 'dns'   
            %
            %   Note that this depends on functions from the RealTimeAnalysis
            %   code base for fitting.
            %
            %   TODO: Reduce dependency on class by either making densities
            %   arguments or defining them at the top of the function in
            %   terms of class members (so it is clear how to remove all
            %   the class dependency from the function).
            %   TODO: everything in units of U.
            %   TODO: cleanup fitting code, reduce dpeendencies, etc.
            %   TODO: Merge this function with fitSinglesAtomicLimit. Add a
            %   text-based argument to select singles vs. doubles.

            % TODO: write this function elsewhere ... could make it similar
            % to fit_fg.m or fit_dqmc_attractive.m
            if ~exist('mode','var')
                mode = 'singles';
            end
            
            if strcmp(mode,'singles')
                fn = 'atomicLimitRadial_SinglesDensity1D';
            elseif strcmp(mode,'doubles')
                fn = 'atomicLimitRadial_DoublesDensity1D';
            elseif strcmp(mode,'ups')||strcmp(mode,'dns')
                error('modes "ups" and "dns" not yet implemented.');
            else
                error('unsupported value for mode variables');
            end
            
            %guess fit parameters
            U = obj.h*UKHz*1e3;
            T = 0.1; %units of U
            Beta = 1/T;
            Omega = 2*pi*400*sqrt(obj.mLi/abs(U))*(750e-9/obj.a); %units of sqrt(m/U)
            AvgMu = 1;
            DeltaMu = 0;
            
            %fitting
            InitP = [Beta,Omega,AvgMu,DeltaMu,0.96,sign(UKHz)];
            FixedP = [0,0,0,1,1,1];
            [Fp,~,FFH,SE] = fit1D(obj.BinAvg*obj.a,obj.Occs_AzAvg,1./obj.Occs_AzAvgUnc.^2,fn,InitP,FixedP);
            Xinterp = linspace(min(obj.BinAvg),max(obj.BinAvg),100);
            %parse results
            TU = 1/Fp(1);
            TnK = abs(U)/Fp(1)/obj.kb/1e-9;
            TnK_Unc = (1/Fp(1))*SE(1)/Fp(1)*abs(U)/obj.kb/1e-9;
            OmegaMean = Fp(2)/sqrt(obj.mLi/abs(U))*obj.a/750e-9;
            OmegaMeanUnc = SE(2)/sqrt(obj.mLi/abs(U))*obj.a/750e-9;
            AspectRatio = max([obj.GaussFitParams(3)/obj.GaussFitParams(4),obj.GaussFitParams(4)/obj.GaussFitParams(3)]);
            OmegaStrong = Fp(2)/sqrt(obj.mLi/abs(U))*obj.a/750e-9*sqrt(AspectRatio); 
            OmegaWeak = Fp(2)/sqrt(obj.mLi/abs(U))*obj.a/750e-9/sqrt(AspectRatio);
            AvgMu = Fp(3);
            DeltaMu = Fp(4);
            Ps = tanh(Fp(1)*Fp(4));
                      
            %store info
            obj.OmegaMean = OmegaMean;
            obj.ChemPot_0 = AvgMu*abs(U);
            obj.ChemPot = AvgMu*abs(U) - 0.5 * obj.mLi * OmegaMean^2 * (obj.BinAvg * obj.a).^2;
        
            %display fit results
            fprintf('U = %0.1f KHz \n', U/obj.h/1e3)
            fprintf('T = %0.3f U = %0.1f +/- %0.1f nK \n', TU, TnK, TnK_Unc);
            fprintf('OmegaMean = (2pi) %0.1f +/- %0.1f Hz \n', OmegaMean/(2*pi), OmegaMeanUnc/(2*pi));
            fprintf('OmegaWeak = (2pi) %0.1f Hz \n', OmegaWeak/(2*pi));
            fprintf('OmegaStrong = (2pi) %0.1f Hz \n', OmegaStrong/(2*pi));
            fprintf('AvgMu/U = %0.3f +/- %0.3f \n', AvgMu,SE(3));
            fprintf('DeltaMu/U = %0.3f +/- %0.3f \n', DeltaMu,SE(4));
            fprintf('P^s = %0.3f \n', Ps);
            fprintf('Density Scale Factor = %0.3f +/- %0.3f \n\n', Fp(5), SE(5));

            %Create figure
            fh = figure('name','Fit to atomic limit');
            subplot(2,2,1)
                errorbar(obj.BinAvg, obj.Occs_AzAvg, obj.Occs_AzAvgUnc, 'ro');
                hold on;
                %plot(Xinterp,atomicLimitRadial_SinglesDensity1D([Beta,Omega,AvgMu,DeltaMu,0.96],Xinterp),'b');
                plot(Xinterp, FFH(Xinterp*obj.a),'b');
                hold off;
                grid on;
                xlabel('Lattice Sites')
                ylabel('Filling')
            subplot(2,2,2)
                errorbar(obj.ChemPot/abs(U),obj.Occs_AzAvg,obj.Occs_AzAvgUnc,'ro');
                hold on;
                plot(AvgMu - 0.5 * obj.mLi * OmegaMean^2 * (Xinterp * obj.a) .^2 / abs(U), FFH(Xinterp * obj.a), 'b');
                grid on;
                xlabel('ChemPot/|U|')
                ylabel('Filling')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Specialized analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Compressibility, CompressibilityUnc,...
                nInBetween, nInBetweenUnc, rInBetween, rInBetweenUnc] ...
                = getCompressibility(obj, display_results,...
                save_results, save_dir, n, nunc, r, runc)
            %getCompressibility Compute density compressibility assuming a
            %harmonic trapping potential.
            %
            %   [Compressibility,CompressibilityUnc] = getCompressibility(obj)
            %
            %   TODO: Make density and radius and etc. arguments?
            
            %read in values
            if ~exist('display_results', 'var') || isempty(display_results)
                display_results = 1;
            end
            
            if ~exist('save_results','var') || isempty(save_results)
                save_results = 0;
            end
            
            if ~exist('save_dir','var') || isempty(save_dir)
                save_dir = pwd;
            end
            
            if ~exist('n','var')
                n = obj.Occs_AzAvg;
            end
            
            if ~exist('nunc','var')
                nunc = obj.Occs_AzAvgUnc;
            end
            
            if ~exist('r','var')                  
                r = obj.BinAvg;
            end
            
            if ~exist('runc','var')
                runc = obj.BinUnc;
            end
            %interpolate points where we will be calculating the
            %compressiblity
            nInBetween = 0.5*(n(2 : end) + n(1 : end - 1));
            nInBetweenUnc = 0.5 * sqrt(nunc(2 : end).^2 + nunc(1 : end - 1).^2);
            ndiff = n(2:end) - n(1:end - 1);
            ndiffunc = 2 * nInBetweenUnc;

            rInBetween = 0.5 * (r(2:end) + r(1:end-1));
            rInBetweenUnc = 0.5 * (runc(2:end).^2 + runc(1:end-1).^2);
            rdiff = (r(2:end) - r(1:end-1));
            rdiffunc = 2 * rInBetweenUnc;
            %numerical derivative
            dndr = ndiff ./ rdiff;
            dndrunc = dndr .* sqrt( (ndiffunc./ndiff).^2 + (rdiffunc./rdiff).^2 );

            %dn/du = - dn/dr * 1/r * 1/(m*omega^2)
            Compressibility =-dndr .* 1 ./ (rInBetween);
            CompressibilityUnc = abs(Compressibility) .* ...
                sqrt((dndrunc./dndr).^2 + (rInBetweenUnc./rInBetween).^2);
            
            if display_results
                fh = figure('name','Compressibility');
                subplot(1,2,1)
                errorbar(rInBetween, Compressibility, CompressibilityUnc,...
                    CompressibilityUnc, rInBetweenUnc, rInBetweenUnc,'ro');
                grid on;
                xlabel('Position (Lattice Sites)');
                ylabel('\kappa =dn/d\mu = -dn/dr 1/r (1/m\omega^2 a^2)');
                
                subplot(1,2,2)
                errorbar(nInBetween, Compressibility, CompressibilityUnc,...
                    CompressibilityUnc, rInBetweenUnc, rInBetweenUnc,'ro');
                grid on;
                xlabel('<n>');
                ylabel('\kappa =dn/d\mu = -dn/dr 1/r (1/m\omega^2 a^2)');
                ttl = sprintf('%s\n Compressibility, assuming harmonic trap', obj.identifier);
                suptitle(ttl);
                if save_results
                    fpath = fullfile(save_dir,sprintf('%s_compressibility.fig',obj.identifier));
                    savefig(fh,fpath);
                    fprintf('Saved compressibility data to %s\n',fpath);
                end
            end
        end
        
        function [DomainNumber, DomainParity] = identifyDomains(obj, Img)
            %identifyDomains
            %
            %   [DomainNumber,DomainParity] = identifyDomains(obj,Img)
            %   Identify domains for a single picture. DomainNum gives a
            %   unique number for all points in a given domain. It is an array
            %   the same size as Img. DomainParity gives the parity for each
            %   domain. +1 corresponds to domains that would have a 1 in the
            %   position (y,x) = (1,1), -1 domains would have a zero there.
            %   DomainParity is the same size as Img.
            %
            %   TODO: Work in progress. Not sure what state I left this in.
            
            [X,Y] = meshgrid(1:size(obj.Occs_Stack, 2), 1:size(obj.Occs_Stack, 1));
            SiteParity = mod(X+Y, 2);
%             Img = obj.Occs_Stack(:,:,1);
            DomainNumber = zeros(size(Img)); %give each domain a unique number
             
            %loop over all points in the picture, beginning with the upper
            %left corner.
            for ii = 1:size(Img, 2)
                for jj = 1:size(Img, 1)
                    if ii == 1 && jj == 1
                        %first pixel.
                        DomainNumber(jj,ii) = 1;
                    else
                        %check if the next pixel is in the same domain as
                        %the previous y-pixel. If so, give it the same
                        %domain number.
                        if jj~=1 && (Img(jj,ii) ~= Img(jj-1, ii))
                            DomainNumber(jj, ii) = DomainNumber(jj-1, ii);
                            %now need to check the previous x-pixel. It may
                            %also be in the same domain, and now that we
                            %know that we must make sure those domain
                            %numbers are the same.
                            if ii~=1 &&(Img(jj, ii) ~= Img(jj, ii-1))
                                %if this is the case, reassign the domain
                                %number.
                                DomainNumber(DomainNumber==DomainNumber(jj, ii)) = ...
                                    DomainNumber(jj, ii-1);
                            end
                        %if it was not the same domain as the previous
                        %y-pixel, check the previous x-pixel.
                        elseif ii~=1 && (Img(jj, ii) ~= Img(jj, ii-1))
                            DomainNumber(jj, ii) = DomainNumber(jj, ii-1);
                        %if pixel was not in the same domain as any
                        %previously processed pixels, give it a new number.
                        else
                            DomainNumber(jj, ii) = max(DomainNumber(:)) + 1;
                        end
                    end
                end
            end
            DomainParity = 2*xor(Img,SiteParity) - 1; %+1 if corresponds to having a 1 in the ULC at position (y,x) = (1,1), -1 otherwise
            
            %remove domains with only a single element %no real need to do
            %this...will be taken care of by azimuthal average.
%             [NRepeats,Centers] = hist(DomainNumber(:),unique(DomainNumber));
%             UniqueMultiplePtDomains = unique(DomainNumber);
%             UniqueMultiplePtDomains = UniqueMultiplePtDomains(NRepeats>1);
%             
%             DomainNumberMod = DomainNumber;
%             for ii = 1:numel(DomainNumberMod)
%                 if sum(DomainNumberMod(ii)==UniqueMultiplePtDomains)==0
%                     DomainNumberMod(ii) = 0;
%                 end
%             end
%             
%             DomainNum = DomainNumberMod + min(DomainNumberMod(DomainNumberMod~=0));
        end
        
        function BinContainer = domainAnalysis(obj,DistGrid,BinEdges)
            %domainAnalysis
            %
            %   BinContainer = domainAnalysis(obj,DistGrid,BinEdges)
            %   Domain analysis for all pictures. BinContainer is a cell
            %   array, with Ncells = NBins. Each cell contains an Nx4 vector
            %   containing information about the domains in a given bin for
            %   all pictures in the dataset. Each line of this vector has the
            %   form [X_CenterOfMass,YCom,Domain_Size,NDomainPts,Domain_Parity]. See
            %   identifyDomains() for a description of domain parity.
            %   
            %   TODO: Work in progress. Not sure what state I left this in.
            
            tic;
            BinContainer = {};
            for ii = 1:(length(BinEdges)-1)
                BinContainer{ii} = [];
            end
            
            
            [X,Y] = meshgrid(1:size(obj.Occs_Stack,2), 1:size(obj.Occs_Stack,1)); %except don't want real X,Y position, but scaled appropriately?
            for ii = 1:size(obj.Occs_Stack,3)
                Img = obj.Occs_Stack(:,:,ii);
                [DomainNum,DomainParity] = obj.identifyDomains(Img);
               %get center of mass of each domain...can do everything in a
               %single call to azAvg_General by stacking things. This may
               %be much harder to read. Essentially, think of as a separate
               %call for each layer of the stack, i.e. X,Y and
               %DomainParity. But faster because we only have to do the
               %overhead for azimuthal averaging once.
               [~,~,Pos,SDM,~,NPtsDomain,~] = azAvg_General(cat(3,X,Y,DomainParity),[],DomainNum,0.5:1:(max(DomainNum(:))+0.5));
               Xcom = Pos(:,1);
               Xsize = SDM(:,1).*sqrt(NPtsDomain);
               Xsize = Xsize(~isnan(Xcom));
               Ycom = Pos(:,2);
               Ycom = Ycom(~isnan(Xcom));
               Ysize = SDM(:,2).*sqrt(NPtsDomain);
               Ysize = Ysize(~isnan(Xcom));
               DomParity = Pos(:,3);
               DomParity = DomParity(~isnan(Xcom));
               NPtsDomain = NPtsDomain(~isnan(Xcom));
               Xcom = Xcom(~isnan(Xcom));
               
               %x-center
%                [D,~,Xcom,XcomSDM,~,PtsInBin,~] = azAvg_General(X,[],DomainNum,0.5:1:(max(DomainNum(:))+0.5));
%                NPtsDomain = PtsInBin(~isnan(Xcom));
%                Xsize = XcomSDM(~isnan(Xcom)).*sqrt(NPtsDomain);
%                Xcom = Xcom(~isnan(Xcom));
%                %y-center
%                [~,~,Ycom,YcomSDM,~,~,~] = azAvg_General(Y,[],DomainNum,0.5:1:(max(DomainNum(:))+0.5));
%                Ysize = YcomSDM(~isnan(Ycom)).*sqrt(NPtsDomain);
%                Ycom = Ycom(~isnan(Ycom));
%                %parity
%                [~,~,DomParity,~,~,~,~] = azAvg_General(DomainParity,[],DomainNum,0.5:1:(max(DomainNum(:))+0.5));
%                DomParity = DomParity(~isnan(DomParity));
               Rsize = sqrt(Xsize.^2+Ysize.^2);
            
               %store info
               for jj = 2:length(BinEdges)
                   LinearIndices = sub2ind(size(DistGrid),round(Ycom),round(Xcom));
                   ChosenPts = DistGrid(LinearIndices)<BinEdges(jj) & DistGrid(LinearIndices)>=BinEdges(jj-1);
                   NewInfo = [Xcom(ChosenPts),Ycom(ChosenPts),Rsize(ChosenPts),NPtsDomain(ChosenPts),DomParity(ChosenPts)];
                   BinContainer{jj-1} = cat(1,BinContainer{jj-1},NewInfo);
               end
               
               if 0
                   %display this info for debugging.
                    figure;
                    imagesc(DomainNum)
                    axis equal; axis image; hold on;
                    scatter(Xcom,Ycom,'rx')
                    
                    figure;
                    histogram(Rsize,linspace(0,3,15));
                end
            
            end
            EndT = toc;
            fprintf('Domain analysis for one folder to %0.2f s \n',EndT);
            
        end
        
        function [DomainMeanSize, DomainStdSize, MeanNumberInDomain, NumberDomains]...
                = showDomainsVsBin(obj, BinContainer)
            %showDomainsVsBin
            %
            %   [DomainMeanSize,DomainStdSize,MeanNumberInDomain,NumberDomains] = showDomainsVsBin(obj,BinContainer)
            %   Takes the same type of object output by domainAnalysis()
            %   function. should maybe clean this fn up. domainAnalysis() does
            %   not assume any specific type of DistGrid/BinEdges. This one
            %   assumes whatever ones we set when instantiating object.
            %
            %   TODO: Work in progress. Not sure what state I left this in.
            
            DomainMeanSize = zeros(obj.NBins,1);
            DomainStdSize = zeros(obj.NBins,1);
            MeanNumberInDomain = zeros(obj.NBins,1);
            NumberDomains = zeros(obj.NBins,1);
            for ii = 1:obj.NBins
                DomainMeanSize(ii) = mean(BinContainer{ii}(:,3));
                DomainStdSize(ii) = std(BinContainer{ii}(:,3));
                MeanNumberInDomain(ii) = mean(BinContainer{ii}(:,4));
                NumberDomains(ii) = size(BinContainer{ii},1);
            end
            figure;
            errorbar(obj.RadialPos,DomainMeanSize,DomainStdSize,'bo')
            grid on;
            xlabel('Bin Position (sites)')
            ylabel('Domain Size (sites)')
        end
        
        function DpParity_Corr_c = getDomainCorrelators(obj, DomainParityStack)
            %getDomainCorrelators Get domain correlator for each image.
            %
            %   DpParity_Corr_c = getDomainCorrelators(obj,DomainParityStack)
            %
            %   TODO: Work in progress. Not sure what state I left this in.
            
            [DomainParityCorr, DomainParityExpanded, DomainParityShifted] = ...
                obj.getFullCorrelatorData(DomainParityStack, 4);
            
            %azimuthally average each image. %TODO replace the redundant
            %averaging as I've done elsewhere.
            [~,~,DpCorr_AzAvg,DpCorr_AzAvgUnc,~,~,~] = ...
                azAvg_General(mean(DomainParityCorr,3),[],obj.DistGrid,obj.BinEdges);
            [~,~,DpExpanded_AzAvg,~,~,~,~] = ...
                azAvg_General(mean(DomainParityExpanded,3),[],obj.DistGrid,obj.BinEdges);
            [~,~,DpShifted_AzAvg,~,~,~,~] = ...
                azAvg_General(mean(DomainParityShifted,3),[],obj.DistGrid,obj.BinEdges);
%             Dp_AzAvg = azAvg_General(mean(DomainParityStack,3),[],obj.DistGrid,obj.BinEdges);
%             
%             DpShifted_AzAvg = repmat(Dp_AzAvg,[]);
%             DpExpande_AzAvg = 
            %get final correlator
%             obj.Density_Corr_AzAvg = obj.Corr_AzAvg - obj.Occs_Expanded_AzAvg.*obj.Occs_Shifted_AzAvg;
            DpParity_Corr_c = squeeze(DpCorr_AzAvg - DpExpanded_AzAvg.*DpShifted_AzAvg);
            DpParity_Corr_c = permute(DpParity_Corr_c,[2,3,1]);
            [DpParity_Corr_c,~] = df.getSymmetricCorrMat(DpParity_Corr_c,[]);
            
        end
        
        function [FigHandle, XPos, DensX, DensXUnc, ...
                NNCorrHorz_X, NNCorrHorzUnc_X,...
                LineFitParamsX, LineFitStdErrX,...
                YPos, DensY, DensYUnc,...
                NNCorrVert_Y, NNCorrVertUnc_Y,...
                LineFitParamsY, LineFitStdErrY,...
                Dens, DensSDM, DensSD] ...
                = sumDirection(obj, Angle, Cx, Cy, CropR, BinWidth, ComputeCorrelator)
            %sumDirection
            %
            %   [XPos,DensX,DensXUnc,CorrA_X,CorrAUnc_X,LineFitParamsX,LineFitStdErrX,...
            %   YPos,DensY,DensYUnc,CorrB_Y,CorrBUnc_Y,LineFitParamsY,LineFitStdErrY,...
            %   Dens,DensSDM,DensSD] = sumDirection(obj,Angle,Cx,Cy,CropR,BinWidth,ComputeCorrelator,Display)
            %   Rotate the average image and sum along the X direction. Can
            %   also compute correlator, and etc.
            %   Angle, is the angle in radians to rotate the iamge
            %   Cx and Cy are the coordinates of the center about which we
            %   rotate. These are `relative coordinates' i.e. in the space of
            %   our cropped images. To make them absolute between different
            %   data sets, you must also fix CroppedPicStartCoords.
            %   CropR is the maximum radius to be used in sums. Anything
            %   outside of this radius is ignored
            %   ComputeCorrelator is a boolean specifying whether or not to
            %   compute the correlators in this region. If it is zero, the
            %   correlator values that are returned will be zeros
            %   Display is a boolean specifying whether or not to plot the
            %   results
            %
            %   TODO: Cleanup and simplify
            %   TODO: There is some error (?) which causes apparent
            %   imbalance between nearest neighbor correlators along
            %   different directions in either binning scheme. I.e. C(0,1)
            %   and C(1,0) look very imbalanced. But they cannot be very
            %   imbalanced, because the 2D correlator images are rather
            %   balanced.

            if ~exist('Cx','var')
                Cx = round(obj.Cx_AzAvg);
            end
            
            if ~exist('Cy','var')
                Cy = round(obj.Cy_AzAvg);
            end
            
            if isempty(Cx)
                 Cx = round(obj.Cx_AzAvg);
            end
            
            if isempty(Cy)
                Cy = round(obj.Cy_AzAvg);
            end
            
            if ~exist('CropR','var')
                CropR = 14;
            end
            
            if ~exist('BinWidth','var')
                BinWidth = 1;
            end
            
            if ~exist('ComputeCorrelator','var')
                ComputeCorrelator = 0;
            end
            
            if ~exist('Display','var')
                Display = 1;
            end
            
            %construct X and Y coordinates
            [Xs,Ys] = meshgrid(1:size(obj.Occs_ImgAvg,2),1:size(obj.Occs_ImgAvg,1));
            Rs = sqrt((Xs-Cx).^2+(Ys-Cy).^2);
            
            %crop off everything outside of a certain radius
            CropMask = double(obj.DistGrid < CropR);
            NanMask = CropMask;
            NanMask(NanMask==0)=nan;
            ImgCrop = obj.Occs_ImgAvg .* NanMask;

            %get Density average and SDM over the cropped region
            [~, ~, Dens, DensSDM, ~, NPts, ~] = azAvg_General(ImgCrop,[],zeros(size(ImgCrop)),[-1,1]);
            DensSD = DensSDM.*sqrt(NPts);
            
            %rotate the coordinates
            RotXs = Xs*cos(Angle)-Ys*sin(Angle);
            RotYs = Ys*cos(Angle) + Xs*sin(Angle);
            
            RotXsCrop = RotXs;
            RotYsCrop = RotYs;
            
            %use rotated coordinates to construct bins
%             BinWidth = 1;
            BinEndPtsX = [floor(min(RotXs(:))):BinWidth:ceil(max(RotXs(:)))];
            BinCentersX = 0.5*(BinEndPtsX(2:end) + BinEndPtsX(1:end-1));
            BinEndPtsY = [floor(min(RotYs(:))):BinWidth:ceil(max(RotYs(:)))];
            
            w = warning ('off','all');
            [XPos, ~, DensX, DensXUnc, ~ ,NPtsInBin_X, MaskStackX] = ...
                azAvg_General(ImgCrop, [], RotXsCrop, BinEndPtsX);
            %Remove Nans
            NPtsInBin_X = NPtsInBin_X(~isnan(DensX));
            XPos = XPos(~isnan(DensX));
            DensXUnc = DensXUnc(~isnan(DensX));
            DensXUnc(DensXUnc == 0) = 1;
            DensX = DensX(~isnan(DensX));
            
            %fit to get average slope           
            linefn = @(P,X) P(1) + P(2)*X;
            difffnX = @(P) (linefn(P,XPos)-DensX)./DensXUnc;
            InitPX = [mean(DensX),0];
            LBsX = [0,-inf];
            UBsX = [2,inf];
            
            try
                opts1=  optimset('display','off');
                [LineFitParamsX,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                 lsqnonlin(difffnX,InitPX,LBsX,UBsX,opts1);
            catch
                LineFitParamsX = [0,0];
            end
            try
                CI = nlparci(LineFitParamsX,residual,'jacobian',jacobian); %95%Confidence intervals
                LineFitStdErrX = (CI(:,2)-CI(:,1))/3.92;
            catch
                LineFitStdErrX = zeros(size(LineFitParamsX));
            end
            
            if ComputeCorrelator
                %Compute X correlator                
                n_stack_crop = obj.Occs_Stack.*repmat(NanMask,[1,1,size(obj.Occs_Stack,3)]);
                AzAvgGrid = RotXsCrop;
                BinEndPts = BinEndPtsX;
                NumNeighbors = obj.NumNeighbors;
                RestrictSameBin = obj.OnlyCorrelateSitesSameBin;
                
                [n,nsdm,r,rsdm,nptsbin,nI,nIsdm,nA,nAsdm,nIsdmA,nAsdmI] = obj.getDensMoments(n_stack_crop,AzAvgGrid,BinEndPts);
                [nn,nnsdm,ni,nisdm,nj,njsdm,npts] = obj.getCorrMoments(n_stack_crop,AzAvgGrid,BinEndPts,NumNeighbors,RestrictSameBin);
                [nnc,nncunc] = obj.getCorrWithErr(nn,ni,nj,npts);
                
                XPos = r;
                DensX = n;
                DensXUnc = nsdm;
                Corr_AzAvg_X = permute(nnc,[2,3,1]);
                Sqrt_Var_SampleCovariance_X = permute(nncunc,[2,3,1]);
                
                [~,~,Corr_AzAvg_X_Test,Corr_AzAvg_X_TestUnc,~,~,~] = azAvg_General(obj.nnc2D,[],RotXsCrop,BinEndPtsX);
            else
                Corr_AzAvg_X = zeros([2*obj.NumNeighbors+1,2*obj.NumNeighbors+1,length(DensX)]);
                Sqrt_Var_SampleCovariance_X = Corr_AzAvg_X;
            end
            
            %Y (orthogonal direction) azimuthal average
            [YPos, ~, DensY, DensYUnc, ~, NPtsInBin_Y, MaskStackY] = ...
                            azAvg_General(ImgCrop,[],RotYsCrop,BinEndPtsY);
            %Remove Nans
            NPtsInBin_Y = NPtsInBin_Y(~isnan(DensY));
            DensYUnc = DensYUnc(~isnan(DensY));
            DensYUnc(DensYUnc == 0) = 1;
            YPos = YPos(~isnan(DensY));
            DensY = DensY(~isnan(DensY));
            
            %Average slope along Y           
            difffnY = @(P) (linefn(P,YPos)-DensY)./DensYUnc;
            InitPY = [mean(DensY),0];
            LBsY = [0,-inf];
            UBsY = [2,inf];
            
            try
                opts1=  optimset('display','off');
                [LineFitParamsY,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                    lsqnonlin(difffnY,InitPY,LBsY,UBsY,opts1);
            catch
                LineFitParamsY = [0,0];
            end
            
            try
                CI = nlparci(LineFitParamsX,residual,'jacobian',jacobian);
                LineFitStdErrY = (CI(:,2)-CI(:,1))/3.92;
            catch
                LineFitStdErrY = zeros(size(LineFitParamsX));
            end
            
            if ComputeCorrelator
                %Y correlators
                AzAvgGrid = RotYsCrop;
                BinEndPts = BinEndPtsY;
                NumNeighbors = obj.NumNeighbors;
                RestrictSameBin = obj.OnlyCorrelateSitesSameBin;
            
                w = warning ('on','all');
                [n, nsdm, r, rsdm, nptsbin, nI, nIsdm, nA, nAsdm, nIsdmA, nAsdmI] = ...
                    obj.getDensMoments(n_stack_crop, AzAvgGrid,BinEndPts);
                [nn, nnsdm, ni, nisdm, nj, njsdm, npts] = ...
                    obj.getCorrMoments(n_stack_crop, AzAvgGrid, BinEndPts, NumNeighbors, RestrictSameBin);
                [nnc, nncunc] = obj.getCorrWithErr(nn,ni,nj,npts);
                [g2, g2unc] = obj.getG2WithErr(nn,ni,nj,npts);
                [nnI, nnIsdm, niI, niIsdm, njI, njIsdm, ~] = ...
                    obj.getCorrMomentsImgAvg(n_stack_crop, AzAvgGrid, BinEndPts, NumNeighbors, RestrictSameBin);
                [nnIcA, nnIcAsdm] = obj.getSupplCorr(nnI, njI, AzAvgGrid, BinEndPts, RestrictSameBin);
                              
                YPos = r;
                DensY = n;
                DensYUnc = nsdm;
                Corr_AzAvg_Y = permute(nnc, [2,3,1]);
                Sqrt_Var_SampleCovariance_Y = permute(nncunc, [2,3,1]);  
                
                %another test correlator
                [~, ~, Corr_AzAvg_Y_Test, Corr_AzAvg_Y_TestUnc, ~, ~, ~] =...
                    azAvg_General(obj.nnc2D, [], RotYsCrop, BinEndPtsY);
            else
                Corr_AzAvg_Y = zeros([2*obj.NumNeighbors+1, 2*obj.NumNeighbors+1, length(DensY)]);
                Sqrt_Var_SampleCovariance_Y = Corr_AzAvg_Y;
            end
            
            %assign correlators for plotting
            NNCorrHorz_X = squeeze(Corr_AzAvg_X(obj.CenterIndex_CorrMatrix,obj.CenterIndex_CorrMatrix+1,:));
            NNCorrHorzUnc_X = squeeze(Sqrt_Var_SampleCovariance_X(obj.CenterIndex_CorrMatrix,obj.CenterIndex_CorrMatrix+1,:));
            NNCorrVert_X = squeeze(Corr_AzAvg_X(obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix,:));
            NNCorrVertUnc_X = squeeze(Sqrt_Var_SampleCovariance_X(obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix,:));
            NNCorrMean_X = 0.5*(NNCorrHorz_X + NNCorrVert_X);
            NNCorrMeanUnc_X = 0.5*sqrt(NNCorrHorzUnc_X.^2 + NNCorrVertUnc_X.^2);
            
            NNNCorr1_X = squeeze(Corr_AzAvg_X(obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix+1,:));
            NNNCorr1Unc_X = squeeze(Sqrt_Var_SampleCovariance_X(obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix+1,:));
            NNNCorr2_X = squeeze(Corr_AzAvg_X(obj.CenterIndex_CorrMatrix-1,obj.CenterIndex_CorrMatrix+1,:));
            NNNCorr2Unc_X = squeeze(Sqrt_Var_SampleCovariance_X(obj.CenterIndex_CorrMatrix-1,obj.CenterIndex_CorrMatrix+1,:));
            NNNCorrMean_X = 0.5*(NNNCorr1_X + NNNCorr2_X);
            NNNCorrMeanUnc_X = 0.5*sqrt(NNNCorr1Unc_X.^2 + NNNCorr2Unc_X.^2);
            
            NNCorrHorz_Y = squeeze(Corr_AzAvg_Y(obj.CenterIndex_CorrMatrix,obj.CenterIndex_CorrMatrix+1,:));
            NNCorrHorzUnc_Y = squeeze(Sqrt_Var_SampleCovariance_Y(obj.CenterIndex_CorrMatrix,obj.CenterIndex_CorrMatrix+1,:));
            NNCorrVert_Y = squeeze(Corr_AzAvg_Y(obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix,:));
            NNCorrVertUnc_Y = squeeze(Sqrt_Var_SampleCovariance_Y(obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix,:));
            NNCorrMean_Y = 0.5*(NNCorrHorz_Y + NNCorrVert_Y);
            NNCorrMeanUnc_Y = 0.5*sqrt(NNCorrHorzUnc_Y.^2 + NNCorrVertUnc_Y.^2);
            
            NNNCorr1_Y = squeeze(Corr_AzAvg_Y(obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix+1,:));
            NNNCorr1Unc_Y = squeeze(Sqrt_Var_SampleCovariance_Y(obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix+1,:));
            NNNCorr2_Y = squeeze(Corr_AzAvg_Y(obj.CenterIndex_CorrMatrix-1,obj.CenterIndex_CorrMatrix+1,:));
            NNNCorr2Unc_Y = squeeze(Sqrt_Var_SampleCovariance_Y(obj.CenterIndex_CorrMatrix-1,obj.CenterIndex_CorrMatrix+1,:));
            NNNCorrMean_Y = 0.5*(NNNCorr1_Y + NNNCorr2_Y);
            NNNCorrMeanUnc_Y = 0.5*sqrt(NNNCorr1Unc_Y.^2 + NNNCorr2Unc_Y.^2);
            
            FigHandle = figure('name','sumDirection');
            NRows = 4;
            NCols = 6;
            
            subplot(NRows, NCols, 1);
            imagesc(CropMask .* obj.Occs_ImgAvg, [-1,1]);
            axis equal;
            axis image;
            ttl = sprintf('Density = %0.02f +/- %0.02f\n SD = %0.02f',...
                Dens, DensSDM, DensSD);
            title(ttl);
            colorbar;
            
            subplot(NRows, NCols, NCols + 1 );
            imagesc(squeeze(obj.nnc2D(:,:,obj.CenterIndex_CorrMatrix,obj.CenterIndex_CorrMatrix+1)) .*...
                CropMask, [-0.07,0.07]);
            axis equal;
            axis image;
            if isempty(obj.ColorMap)
                obj.getColorMap();
            end
            colormap(obj.ColorMap);
            title('C(0,1) = H');
            colorbar;
            
            % display bins used for averaging
            subplot(NRows, NCols, 2)
            imagesc(obj.showMasks(MaskStackX .* NanMask));
            axis equal;
            axis image;
            title('X avg bins');
            
            subplot(NRows, NCols, NCols + 2)
            plot(XPos, NPtsInBin_X);
            title('Npts in X bin');
            
            subplot(NRows, NCols, 2*NCols + 2)
            imagesc(obj.showMasks(MaskStackY .* NanMask));
            axis equal;
            axis image;
            title('Y avg bins');
            
            subplot(NRows, NCols, 3*NCols + 2)
            plot(YPos, NPtsInBin_Y);
            title('Npts in Y bin');
            
            % display correlators
            subplot(NRows, NCols, 2*NCols+1);
            imagesc(squeeze(obj.nnc2D(:,:,obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix)).*CropMask,...
                                        [-0.07,0.07]);
            axis equal;
            axis image;
            colormap(obj.ColorMap);
            title('C(1,0) = V');
            colorbar;
            
            subplot(NRows, NCols, 3*NCols+1);
            CorrDiff = squeeze(obj.nnc2D(:,:,obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix)) ...
                - squeeze(obj.nnc2D(:,:,obj.CenterIndex_CorrMatrix,obj.CenterIndex_CorrMatrix+1));
            CorrDiffNans = CorrDiff.*NanMask;
            MeanDiff = mean(CorrDiff(~isnan(CorrDiffNans)));
            imagesc(CorrDiff.*CropMask,[-0.02,0.02]);
            axis equal;
            axis image;
            colormap(obj.ColorMap);
            colorbar;
            ttl = sprintf('V - H Mean=%0.4f',MeanDiff);
            title(ttl);
            
            subplot(NRows,NCols, [5,6]);
            errorbar(XPos,DensX,DensXUnc);
            hold on;
            XInterp = linspace(min(XPos),max(XPos),100);
            plot(XInterp,linefn(LineFitParamsX,XInterp));
            ttl = sprintf('X Density, Cut A\n Slope = %0.2f +/- %0.2f /20sites',LineFitParamsX(2)*20,LineFitStdErrX(2)*20);
            title(ttl)
            ylim([0,1]);
            grid on;
            
            subplot(NRows,NCols, NCols + [5,6]);
            errorbar(XPos, NNCorrHorz_X, NNCorrHorzUnc_X, 'r.');
            hold on;
            errorbar(XPos,NNCorrVert_X, NNCorrVertUnc_X, 'b.');
            errorbar(XPos,NNCorrMean_X, NNCorrMeanUnc_X, 'k.');
            %                 errorbar(XPos,squeeze(Corr_AzAvg_X_Test(:,obj.CenterIndex_CorrMatrix,obj.CenterIndex_CorrMatrix+1)),...
            %                     squeeze(Corr_AzAvg_X_Test(:,obj.CenterIndex_CorrMatrix,obj.CenterIndex_CorrMatrix+1)),'.m');
            %                  errorbar(XPos,squeeze(Corr_AzAvg_X_Test(:,obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix)),...
            %                     squeeze(Corr_AzAvg_X_Test(:,obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix)),'.');
            plot(XPos,zeros(size(XPos)),'k--');
            grid on;
            legend({'H = C(0,1)','V = C(1,0)','Mean'});
            ax = gca;
            title('NN Corr');
            ylim([-0.08,0.08])
            title('X');
            
            subplot(NRows,NCols,2*NCols+[5,6]);
            errorbar(XPos,NNNCorr1_X,NNNCorr1Unc_X,'r.');
            hold on;
            errorbar(XPos,NNNCorr2_X,NNNCorr2Unc_X,'b.');
            errorbar(XPos,NNNCorrMean_X,NNNCorrMeanUnc_X,'k.');
            plot(XPos,zeros(size(XPos)),'k--');
            grid on;
            legend({'H = C(0,1)','V = C(1,0)'});
            ax = gca;
            title('NNN Corr');
            ylim([-0.08,0.08])
            title('X');
            
            subplot(NRows,NCols,3*NCols+[5,6]);
            %                 errorbar(DensX,NNCorrHorz_X,NNCorrHorzUnc_X,NNCorrHorzUnc_X,DensXUnc,DensXUnc,'.');
            hold on;
            %                 errorbar(DensX,NNCorrVert_X,NNCorrVertUnc_X,NNCorrVertUnc_X,DensXUnc,DensXUnc,'.');
            errorbar(DensX,NNCorrMean_X,NNCorrMeanUnc_X,NNCorrMeanUnc_X,DensXUnc,DensXUnc,'.');
            errorbar(DensX,NNNCorrMean_X,NNNCorrMeanUnc_X,NNNCorrMeanUnc_X,DensXUnc,DensXUnc,'.');
            grid on;
            legend({'V'});
            title('Corr Vs. <n>');
            xlabel('<n>');
            ylabel('Corr');
            legend({'NN','NNN'});
            ylim([-0.08,0.08]);
            title('X');
            
            subplot(NRows,NCols, [3,4])
            errorbar(YPos,DensY,DensYUnc);
            hold on;
            YInterp = linspace(min(YPos),max(YPos),100);
            plot(YInterp,linefn(LineFitParamsY,YInterp));
            
            
            ylim([0,1]);
            grid on;
            ttl = sprintf('Y Density, Cut A\n Slope = %0.2f +/- %0.2f /20sites',LineFitParamsY(2)*20,LineFitStdErrY(2)*20);
            title(ttl)
            
            subplot(NRows,NCols,NCols+[3,4]);
            errorbar(YPos,NNCorrHorz_Y,NNCorrHorzUnc_Y,'r.');
            hold on;
            errorbar(YPos,NNCorrVert_Y,NNCorrVertUnc_Y,'b.');
            errorbar(YPos,NNCorrMean_Y,NNCorrMeanUnc_Y,'k.');
            plot(YPos,zeros(size(YPos)),'k--');
            %                 errorbar(YPos,squeeze(Corr_AzAvg_Y_Test(:,obj.CenterIndex_CorrMatrix,obj.CenterIndex_CorrMatrix+1)),...
            %                     squeeze(Corr_AzAvg_Y_Test(:,obj.CenterIndex_CorrMatrix,obj.CenterIndex_CorrMatrix+1)),'.m');
            %                  errorbar(YPos,squeeze(Corr_AzAvg_Y_Test(:,obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix)),...
            %                     squeeze(Corr_AzAvg_Y_Test(:,obj.CenterIndex_CorrMatrix+1,obj.CenterIndex_CorrMatrix)),'.m');
            plot(YPos,zeros(size(YPos)),'k--');
            grid on;
            legend({'H = C(0,1)','V = C(1,0)','Mean'});
            ax = gca;
            ylim([-0.08,0.08])
            title('Y NN Corr');
            
            subplot(NRows,NCols,2*NCols+[3,4]);
            errorbar(YPos,NNNCorr1_Y,NNNCorr1Unc_Y,'r.');
            hold on;
            errorbar(YPos,NNNCorr2_Y,NNNCorr2Unc_Y,'b.');
            errorbar(YPos,NNNCorrMean_Y,NNNCorrMeanUnc_Y,'k.');
            plot(YPos,zeros(size(YPos)),'k--');
            grid on;
            ax = gca;
            ylim([-0.01,0.05]);
            title('Y NNN Corr');
            
            subplot(NRows,NCols,3*NCols+[3,4]);
            %                 errorbar(DensY,NNCorrHorz_Y,NNCorrHorzUnc_Y,NNCorrHorzUnc_Y,DensYUnc,DensYUnc,'.');
            hold on;
            %                 errorbar(DensY,NNCorrVert_Y,NNCorrVertUnc_Y,NNCorrVertUnc_Y,DensYUnc,DensYUnc,'.');
            errorbar(DensY,NNCorrMean_Y,NNCorrMeanUnc_Y,NNCorrMeanUnc_Y,DensYUnc,DensYUnc,'.');
            errorbar(DensY,NNNCorrMean_Y,NNNCorrMeanUnc_Y,NNNCorrMeanUnc_Y,DensYUnc,DensYUnc,'.');
            grid on;
            legend({'NN','NNN'});
            title('Y Corr Vs. <n>');
            xlabel('<n>');
            ylabel('Corr');
            ylim([-0.08,0.08])
            
            suptitle(obj.identifier);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Display functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fh = showAzAvg(obj, save_results, save_dir, file_name)
            %showAzAvg Display results of azimuthal averaging and image
            %averaging for a stack of images.
            %  
            %   fh = showAzAvg(obj, save_results, save_dir, file_name)
            %
            %   Displays averages over a stack
            %   of images stored in obj.Occs_Stack. Display both 2D image
            %   average and combined image and azimuthal average. Also 
            %   displays the 2D grid and bins used in the azimuthal
            %   average. Also displays a 2D gaussian fit to the average
            %   image. Shows atom number and center of mass versus image.
            %   Also shows reconstruction error vs. image.
            %
            % save_results
            %
            % save_dir
            %
            % file_name
            
            if ~exist('save_results','var')
                save_results = 0;
            end
            if ~exist('save_dir','var') || isempty(save_dir)
                save_dir = pwd;
            end
            if ~exist('file_name','var') || isempty(file_name)
                file_name = sprintf('%s_azavg_summary.fig',obj.identifier);
            end
            
            nrows = 3;
            ncols = 6;
            
            fh = figure('name', 'AzAvg');
            
            % Average image
            h = subplot(nrows, ncols, 1);
                        imagesc(obj.Occs_ImgAvg, [0,1]); 
            axis equal; 
            axis image;

            hold on;
            % plot azimuthal average center
            scatter(obj.Cx_AzAvg, obj.Cy_AzAvg, 'bx');
            % plot gaussian fit center
            scatter(obj.GaussFitParams(1), obj.GaussFitParams(2), 'rx');
            
            % show lattice axes
            axes_len = 30;
            % y-axis
            quiver(obj.Cx_AzAvg, obj.Cy_AzAvg, 0, axes_len, 'g', 'LineWidth', 2);
            % x-like
            quiver(obj.Cx_AzAvg, obj.Cy_AzAvg, axes_len, 0, 'r', 'LineWidth', 2);

            title(sprintf('Average occupation'));
            
            try
                lambda1 = obj.ReconstrInfoStack(1).coordTrafo.gridParameters.lambda(1);
                lambda2 = obj.ReconstrInfoStack(1).coordTrafo.gridParameters.lambda(2);
                theta1 = obj.ReconstrInfoStack(1).coordTrafo.gridParameters.theta1;
                theta2 = obj.ReconstrInfoStack(1).coordTrafo.gridParameters.theta2;
                offset1 = obj.ReconstrInfoStack(1).coordTrafo.indexOffset(1);
                offset2 = obj.ReconstrInfoStack(1).coordTrafo.indexOffset(2);
                beam_angle = 2 * atan(lambda2 / lambda1);
                % image transformed back to fluorescence space
                xform_params = [theta1, theta2, 0, 0, lambda1, lambda2, offset1, offset2, 0, 0];
                FlPic = readimg(char(obj.FlFPaths{1}));
                % approximate fluorescence center
                [x_fl_center, y_fl_center, ~] = latt2img_coord(xform_params, obj.Cx_AzAvg + obj.XStarts_ROI(1) - 1,...
                                                                             obj.Cy_AzAvg + obj.YStarts_ROI(1) - 1);
                % transform reconstruction axes into fluorescence space
                xr_axis = 0:axes_len;
                yr_axis = 0:axes_len;
                [xfl_xr_axis, yfl_xr_axis] = latt2img_coord(xform_params,...
                                                  xr_axis + obj.Cx_AzAvg + obj.XStarts_ROI(1) - 1,...
                                     ones(size(xr_axis)) * (obj.Cy_AzAvg + obj.YStarts_ROI(1) - 1));
                [xfl_yr_axis, yfl_yr_axis] = latt2img_coord(xform_params,...
                                     ones(size(yr_axis)) * (obj.Cx_AzAvg + obj.XStarts_ROI(1) - 1),...
                                                  yr_axis + obj.Cy_AzAvg + obj.YStarts_ROI(1) - 1);

                [Xfl, Yfl] = meshgrid(1:size(FlPic,2),1:size(FlPic,1));   
                
                [Xr_from_fl, Yr_from_fl, ~] = img2latt_coord(xform_params, Xfl, Yfl);

                [Xfrec, Yfrec] = meshgrid(obj.XStarts_ROI(1):obj.XStarts_ROI(1) + obj.ImgCropSize - 1,...
                    obj.YStarts_ROI(1):obj.YStarts_ROI(1) + obj.ImgCropSize - 1);
%                 [Xfrec,Yfrec] = meshgrid(1:size(obj.Occs_ImgAvg,2), 1:size(obj.Occs_ImgAvg,1));
                FullRecInterp = @(Xl,Yl) interp2(Xfrec, Yfrec, obj.Occs_ImgAvg, Xl, Yl);
                Rec_FlSpace = FullRecInterp(Xr_from_fl, Yr_from_fl);

                subplot(nrows, ncols, 2);
                imagesc(Rec_FlSpace, [0, 1]);
                axis equal;
                axis image;
                hold on;
                scatter(x_fl_center, y_fl_center, 'rx');
                scatter(xfl_xr_axis, yfl_xr_axis, 'r.');
                scatter(xfl_yr_axis, yfl_yr_axis, 'g.');
                xlim([x_fl_center - 150, x_fl_center + 150]);
                ylim([y_fl_center - 150, y_fl_center + 150]);
%                 scatter(xfl_xr_axis - nhalf_fl_size, yfl_xr_axis - nhalf_fl_size, 'g.');
%                 scatter(xfl_yr_axis - nhalf_fl_size, yfl_yr_axis - nhalf_fl_size, 'r.');
                title('xform to fl sp');
                
            catch err
                disp(err.message);
                lambda1 = 0;
                lambda2 = 0;
                theta1 = 0;
                theta2 = 0;
                offset1 = 0;
                offset2 = 0;
                beam_angle = 0;
            end
               
            % Display bins used for azimuthal average.
            subplot(nrows, ncols, 7) 
            imagesc(obj.showMasks(obj.MaskStack)); 
            axis equal; 
            axis image;
            title(sprintf('Az Avg Bins, Avg Type = %s \n Cx AzAvg = %0.2f, Cy AzAvg = %0.2f',...
                obj.AzAvgType,obj.Cx_AzAvg,obj.Cy_AzAvg));
            
            % 2D gaussian fit to occupation
            subplot(nrows, ncols, 8)
            imagesc(obj.GaussFit_ImgAvg, [0,1]); 
            axis equal; 
            axis image;
            title(sprintf('Gaussian Fit  \nCx = %0.1f, Cy = %0.1f, Sx = %0.1f, Sy = %0.1f, \n Theta = %0.0f deg',...
                obj.GaussFitParams(1:4), mod(obj.GaussFitParams(6), 2*pi) * 180/pi));
            
            % text info
            ax = subplot(nrows, ncols, [5,6]);
            ax.Visible = 'off';
            text_str = sprintf('ROI center style = %s\n Cropped Pic Start Coords\n X = %0.2f\n Y = %0.2f \n lambda1 = %0.3f pix\n lambda2 = %0.3f pix \n angles from fluorescence y-axis \n reconstruction y-axis = theta1 = %0.1f deg \n reconstruction x-axis = theta2 = %0.1f deg \n angle between lattice axes = %0.2f deg \n offset1 = %d \n offset2 = %d \n angle between lattice beams = 2atan(ax/ay) = %0.3f deg', ...
                obj.CenterStyle, obj.CroppedPicStartCoords(1), obj.CroppedPicStartCoords(2),...
                lambda1, lambda2,...
                theta1 * 180/pi, theta2 * 180/pi, (theta2 - theta1) * 180/pi,...
                offset1, offset2,...
                180 - beam_angle * 180/pi);
            t = text(ax, 0, 0.5, text_str);
            t.Interpreter = 'tex';
            t.FontSize = 10;
            t.FontWeight = 'bold';
            
            % azimuthal average versus bin(radial) position
            subplot(nrows,ncols, [13 14])
            errorbar(obj.BinAvg, obj.Occs_AzAvg, obj.Occs_AzAvgUnc);
            hold on;
            grid on;
            for kk = 1:length(obj.BinEdges)
                plot([obj.BinEdges(kk),obj.BinEdges(kk)],[0,1],'black--');
            end
            hold off;
            ylim([0,1]);
            xlabel('bin grid dist');
            ylabel('Occupation');
            title('Az Avg Density vs. bin grid dist');
            
            subplot(nrows, ncols, [15, 16])
            errorbar(obj.RadialPos,obj.Occs_AzAvg,obj.Occs_AzAvgUnc);
            grid on;
            ylim([0, 1]);
            xlabel('dist from 2d center');
            ylabel('Occupation');
            title('Az Avg Density vs. mean dist from fit center');
            
            
            %Atom number vs shot
            subplot(nrows,ncols, [9 10]);
            plot(obj.PictureNumberInFolder, obj.AtomNumbers, 'bo');
            
            hold on;
            grid on;
            
            upper_ylim = max([0.1, 1.2 * max(obj.AtomNumbers)]);
            ylim([0, upper_ylim]);
            upper_xlim = max(max(obj.PictureNumberInFolder),2);
            xlim([1, upper_xlim]);
            
            % show one standard deviation below and above...
            plot([1, upper_xlim], [obj.MeanAtomNum, obj.MeanAtomNum],'b');
            p1 = plot([1, upper_xlim], [obj.MeanAtomNum + obj.AtomNumSD, obj.MeanAtomNum + obj.AtomNumSD], 'b-.');
            plot([1, upper_xlim], [obj.MeanAtomNum - obj.AtomNumSD, obj.MeanAtomNum - obj.AtomNumSD],'b-.');
            
            % show shot noise limits
            p2 = plot([1, upper_xlim], [obj.MeanAtomNum + sqrt(obj.MeanAtomNum), obj.MeanAtomNum + sqrt(obj.MeanAtomNum)], 'r--');
            plot([1, upper_xlim], [obj.MeanAtomNum - sqrt(obj.MeanAtomNum), obj.MeanAtomNum - sqrt(obj.MeanAtomNum)], 'r--');
            legend([p1, p2], {'stdev', 'shot noise'});
            hold off;
            xlabel('Pic # in folder');
            ylabel('Atom Number');
            title(sprintf('NAtoms = %0.1f +/- %0.1f ',...
                obj.MeanAtomNum, obj.AtomNumSD));
            
            % recnostruction phase
            subplot(nrows, ncols, [11 12]);
            try
                [times, phi1s, phi2s, ~] = obj.getRecPhases();
                
                plot(times, phi1s, '.');
                hold on;
                plot(times, phi2s, '.');
                grid on;
                xlabel('Time (s)');
                ylabel('Phis (Latt Sites)');
                legend({'\phi_1', '\phi_2'});
                title('Reconstruction Phase');
               
            catch err
                disp(err.message);
            end
            
            subplot(nrows,ncols, [3, 4]);
            %Shifts vs. picture
            plot(obj.PictureNumberInFolder, obj.XStarts_ROI, 'bo');
            hold on;
            plot(obj.PictureNumberInFolder, obj.X_CoM_fullpic, 'bx');
            plot(obj.PictureNumberInFolder, obj.YStarts_ROI, 'ro');
            plot(obj.PictureNumberInFolder, obj.Y_CoM_fullpic, 'rx');
            grid on;
            ylabel('Shift (Lattice sites)');
            xlabel('Pic # in folder');
            upper_xlim = max(max(obj.PictureNumberInFolder),2);
            xlim([1,upper_xlim]);
            legend({sprintf('XROI Start, SD = %0.2f', std(obj.XStarts_ROI)),...
                sprintf('XCoM = %0.2f, SD = %0.2f', mean(obj.X_CoM_fullpic), std(obj.X_CoM_fullpic)),...
                sprintf('YROI Start, SD = %0.2f', std(obj.YStarts_ROI)),...
                sprintf('YCoM = %0.2f, SD = %0.2f', mean(obj.Y_CoM_fullpic), std(obj.Y_CoM_fullpic))})
            title(sprintf('Cloud CoM Shift'));
            
            try
                subplot(nrows,ncols, [17 18]);
                plot(obj.PictureNumberInFolder, obj.RecErrRates(:, 1), 'bo');
                hold on;
                plot(obj.PictureNumberInFolder, obj.RecErrRates(:, 2), 'bx');
                grid on;
                ylim([0, 1.5e-2]);
                xlabel('Pic # in folder');
                ylabel('Rec Errs');
                upper_xlim = max(max(obj.PictureNumberInFolder),2);
                xlim([1, upper_xlim]);
                grid on;
                title(sprintf('Rec Errs = %0.2f +/- %0.2f%s and %0.2f +/- %0.2f%s',...
                    mean(obj.RecErrRates(:, 1)) * 1e2, std(obj.RecErrRates(:, 1)) * 1e2, '%',...
                    mean(obj.RecErrRates(:, 2)) * 1e2, std(obj.RecErrRates(:, 2)) * 1e2, '%'));
                legend({'', ''});
            catch
            end
            
            suptitle(obj.identifier);
            
            if save_results
                fpath = fullfile(save_dir, file_name);
                savefig(fh, fpath);
                fprintf('Saved to %s\n', fpath);
            end
        end
        
        function [fighandle, dens, dens_std_devs] = ...
                showFlatness(obj, img, dist_grid, display_extra)
            %test density flatness within a region selected by
            %cutoff_radius.
            
            if ~exist('img', 'var') || isempty(img)
                img = obj.Occs_ImgAvg;
            end
            
            if ~exist('dist_grid', 'var') || isempty(dist_grid)
                dist_grid = obj.DistGrid;
            end
            
            if ~exist('display_extra', 'var') || isempty(display_extra)
                display_extra = 0;
            end
            
            cutoff_rads = [1:10,15:5:35];
            bin_sizes = 2:5;
            n_binning_options = length(bin_sizes) + 1;
            
            dens = zeros([length(cutoff_rads), n_binning_options]);
            % standard deviation of binned density
            dens_std_devs = zeros([length(cutoff_rads), n_binning_options]);
            % expected statistical uncertainty for given amount of binning
            % and number of images
            dens_statistical_unc = zeros([length(cutoff_rads), n_binning_options]);
            % portion of the density standard deviation estimated to be due
            % to disorder
            dens_disorder_std = zeros([length(cutoff_rads), n_binning_options]);
            imgs = cell([length(cutoff_rads), n_binning_options]);
            
            n_imgs = size(obj.Occs_Stack, 3);
            
            for ii = 1:length(cutoff_rads)
                for jj = 1:n_binning_options
           
                    cutoff_radius = cutoff_rads(ii);
                    mask = ones(size(img));
                    mask(dist_grid > cutoff_radius) = nan;
                    img_masked = img .* mask; 
                    
                    if jj > 1
                        nbin = bin_sizes(jj - 1);
                        bin_img = binImg(img_masked, nbin, nbin, 'Normal');
                    else
                        nbin = 1;
                        bin_img = img_masked;
                    end
                    imgs{ii,jj} = bin_img;
                    dens(ii, jj) = mean(bin_img(~isnan(bin_img)));
                    dens_statistical_unc(ii, jj) = sqrt(dens(ii, jj) - dens(ii, jj)^2) / sqrt(n_imgs) / nbin;
                    dens_std_devs(ii, jj) = std(bin_img(~isnan(bin_img)));
                    dens_disorder_std(ii, jj) = sqrt(dens_std_devs(ii, jj)^2 - dens_statistical_unc(ii, jj)^2);
                end
            end
            
            
            % plot azimuthal average results in figure
            fighandle = figure;
            nrows = 2;
            ncols = 2 * n_binning_options; %4;
            
            % create legend
            % first entries are the bin sizes
            leg_sizes = [1, bin_sizes];
            leg = cell(size(leg_sizes));
            for ii = 1:length(leg_sizes)
                leg{ii} = sprintf('bin size = %d', leg_sizes(ii));
            end
            % next entries are the statistical uncertainties for bins of a
            % given size
            for ii = 1:length(leg_sizes)
                leg{ii + length(leg_sizes)} = sprintf('unc for single bin = %d', leg_sizes(ii));
            end
            
            subplot(nrows, ncols, [1:n_binning_options]);
            plot(cutoff_rads, dens);
            grid on;
            ax = gca;
            ax.XLim(1) = 0;
            ax.YLim(1) = 0;
            xlabel('cutoff radius (latt sites)');
            ylabel('mean dens');
            title('Mean density vs. cutoff radius');
            legend(leg);
            
            subplot(nrows, ncols, [2*n_binning_options + 1 : 3*n_binning_options]);
            plot(cutoff_rads, dens_std_devs);
            hold on;
            ax = gca;
            ax.ColorOrderIndex = 1;
            plot([min(cutoff_rads), max(cutoff_rads)], [dens_statistical_unc(nbin, :); dens_statistical_unc(nbin, :)], '--');
            grid on;
            ax = gca;
            ax.XLim(1) = 0;
            xlabel('cutoff radius (latt sites)');
            ylabel('sd');
            title('Standard deviation vs. cutoff radius');
            legend(leg);
            
            subplot(nrows, ncols, [3*n_binning_options + 1 : 4*n_binning_options]);
            plot(cutoff_rads, dens_disorder_std);
            grid on;
            ax = gca;
            ax.XLim(1) = 0;
            xlabel('cutoff radius (latt sites)');
            ylabel('density disorder');
            title('Density disorder vs. cutoff radius');
            legend(leg);
                
            for kk = 1:n_binning_options
                subplot(nrows, ncols, n_binning_options + kk);
                imagesc(imgs{end, kk}, [0,1]);
                axis equal;
                axis image;
                ax = gca;
                ax.XTickLabel = '';
                ax.YTickLabel = '';
                if kk == n_binning_options
                    colorbar;
                end
    %             title('nbin = 1');
            end
            
            explanation_str = 'statistical unc n = sqrt(n - n^2) / sqrt(nimgs) / sqrt(nbin * nbin)';
            suptitle(sprintf('%s\n%s', obj.identifier, explanation_str));
              
            if display_extra
                for jj = 1 :length(bin_sizes) + 1
                    fh1 = figure;
                    nrows = 2;
                    ncols = length(cutoff_rads);

                    for ii = 1:length(cutoff_rads)
                        m = imgs{ii, jj};

                        subplot(nrows, ncols,ii)
                        imagesc(m, [0, 1]);
                        axis equal; axis image;
                        title(sprintf('Cutoff radius = %0.2f', cutoff_rads(ii)));

                        subplot(nrows, ncols, ncols + ii)
                        histogram(m(m > 0));
                        xlim([0,1]);
                        ylabel('counts');
                        xlabel('<n>');
                        title(sprintf('Dens = %0.2f, SD = %0.2f',...
                            dens(ii, jj), dens_std_devs(ii, jj)));
                    end
                    if jj == 1
                        bin_size = 1;
                    else
                        bin_size = bin_sizes(jj - 1);
                    end
                    suptitle(sprintf('Binning = %d', bin_size));
    %                 suptitle(sprintf('Flatness: %s', obj.identifier));

                end
            end
           
        end
        
        function [fig_handle_center] = showPositionStability(obj, show_fits)
            
            if ~exist('show_fits', 'var')
                show_fits = 0;
            end
            
            % organize center of mass data
            x_com = mean(obj.X_CoM_fullpic);
            x_com_std = std(obj.X_CoM_fullpic);
            y_com = mean(obj.Y_CoM_fullpic);
            y_com_std = std(obj.Y_CoM_fullpic);

            [x_com_autocorr, lag_pics] = xcorr(obj.X_CoM_fullpic - x_com, 'unbiased');
            [y_com_autocorr, ~] = xcorr(obj.Y_CoM_fullpic - y_com, 'unbiased');
            
            % also compute higher moments
%             [mx_one, my_one, ~] = get_moment(obj.Occs_Stack, 1, [], []);
           
%             [mx_four, my_four, ~] = get_moment(obj.Occs_Stack, 4, [], []);
%             [mx_five, my_five, ~] = get_moment(obj.Occs_Stack, 5, [], []);
            
            % fit individual pictures
            n_imgs = size(obj.Occs_Stack, 3);

            init_params = [50, 50, 20, 20, 0.8, 0, 0];
            fixed_params = [0, 0, 0, 0, 0, 0, 0];
            gauss_fit_params = zeros(n_imgs, 7);
            gauss_fit_std_errs = zeros(n_imgs, 7);

            [xx, yy] = meshgrid(1:100, 1:100);

            for ii = 1:n_imgs
                [fit_params, ~, ffh, std_errs] = fit2D([], [], obj.Occs_Stack(:, :, ii),...
                    [], {'gaussian2D'}, init_params, fixed_params);

                gauss_fit_params(ii, :) = fit_params;
                gauss_fit_std_errs(ii, :) = std_errs;

                if show_fits
                    figure;
                    subplot(2, 1, 1);
                    imagesc(obj.Occs_Stack(:, :, ii));
                    hold on;
                    scatter(fit_params(1), fit_params(2), 'gx');
                    axis equal;
                    axis image;

                    subplot(2, 1, 2);
                    imagesc(ffh(xx, yy));
                    hold on;
                    scatter(fit_params(1), fit_params(2), 'gx');
                    axis equal;
                    axis image;

                    suptitle(sprintf('%s picture %d \n cx = %0.2f +/- %0.2f, cy = %0.2f +/- %0.2f',...
                        obj.DatasetString, ii, fit_params(1), std_errs(1), fit_params(2), std_errs(2)));
                end
            end

            x_fit_mean = mean(gauss_fit_params(:, 1));
            x_fit_std = std(gauss_fit_params(:, 1));
            y_fit_mean = mean(gauss_fit_params(:, 2));
            y_fit_std = std(gauss_fit_params(:, 2));

            fig_handle_center = figure;
            nrows = 2;
            ncols = 3;
            
            subplot(nrows, ncols, 1);
            plot(1:n_imgs, obj.X_CoM_fullpic - x_com);
            hold on;
            plot(1:n_imgs, gauss_fit_params(:, 1) - x_fit_mean);
            xlabel('picture');
            ylabel('x center');
            % ylim(-3, 3]);
            legend({'com', 'fit center'});
            title(sprintf('x center stability\n com std = %0.2f, fit std = %0.2f', x_com_std, x_fit_std));

            subplot(nrows, ncols, 4)
            hx = histogram(gauss_fit_params(:, 1) - x_fit_mean);
            hold on;
            hx_bin_centers = 0.5 * (hx.BinEdges(2:end) + hx.BinEdges(1:end-1));
            init_params_x = [mean(hx_bin_centers), std(hx_bin_centers), max(hx.Values), 0];
            fixed_params_x = [0,0,0,1];
            % fit histogram to gaussian
            [fit_params_x, ~, ffh_x, std_errs_x, chi_sqr_x] = fit1D(hx_bin_centers, hx.Values,...
                [], {'gaussian1D'}, init_params_x, fixed_params_x);

            x_interp = linspace(min(hx.BinEdges), max(hx.BinEdges), 300);
            plot(x_interp, ffh_x(x_interp));
            title(sprintf('histogram x\n gaussian fit sigma = %0.2f', fit_params_x(2)));


            subplot(nrows, ncols, 2);
            plot(1:n_imgs, obj.Y_CoM_fullpic - y_com);
            hold on;
            plot(1:n_imgs, gauss_fit_params(:, 2) - y_fit_mean);
            xlabel('picture');
            ylabel('y center');
            % ylim([min(gauss_fit_params(:, 2)) - 1, max(gauss_fit_params(:, 2)) + 1]);
            title(sprintf('y center stability\n com std = %0.2f, fit std = %0.2f', y_com_std, y_fit_std));

            subplot(nrows, ncols, 5)
            hy = histogram(gauss_fit_params(:, 2) - y_fit_mean);
            hold on;
            hy_bin_centers = 0.5 * (hy.BinEdges(2:end) + hy.BinEdges(1:end-1));
            init_params_y = [mean(hy_bin_centers), std(hy_bin_centers), max(hy.Values), 0];
            fixed_params_y = [0,0,0,1];
            % fit histogram to gaussian
            [fit_params_y, ~, ffh_y, std_errs_y, chi_sqr_y] = fit1D(hy_bin_centers, hy.Values,...
                [], {'gaussian1D'}, init_params_y, fixed_params_y);

            y_interp = linspace(min(hy.BinEdges), max(hy.BinEdges), 300);
            plot(y_interp, ffh_y(y_interp));
            title(sprintf('histogram y\n gaussian fit sigma = %0.2f', fit_params_y(2)));

            subplot(nrows, ncols, 3)
            plot(lag_pics, x_com_autocorr, 'o-');
            xlabel('picture separation');
            ylabel('auto-correlation');
            title('x-position autocorrelation');
            grid on;
            
            subplot(nrows, ncols, 6)
            plot(lag_pics, y_com_autocorr, 'o-');
            xlabel('picture separation');
            ylabel('auto-correlation');
            title('y-position autocorrelation');
            grid on;
            
            suptitle(sprintf('%s\n cloud position distribution', obj.DatasetString)); 
             
        end
        
        function [fig_handle_skewness, skewness_x, skewness_y, mx_two, my_two] = showSkew(obj)
            n_imgs = size(obj.Occs_Stack, 3);
            
            [mx_two, my_two, ~] = get_moment(obj.Occs_Stack, 2, [], []);
            [mx_three, my_three, ~] = get_moment(obj.Occs_Stack, 3, [], []);
            skewness_x = mx_three ./ sqrt(mx_two) .^3;
            skewness_y = my_three ./ sqrt(my_two) .^3;
            skewness_x_mean = mean(skewness_x);
            skewness_x_std = std(skewness_x);
            skewness_y_mean = mean(skewness_y);
            skewness_y_std = std(skewness_y);
            
             % skewness plot
            fig_handle_skewness = figure;
            nrows = 2;
            ncols = 2;
            
            subplot(nrows, ncols, 1)
            plot(1:n_imgs, skewness_x);
            xlabel('picture');
            ylabel('skewness, x');
            % ylim(-3, 3]);
            title(sprintf('x skewness \n mean = %0.2f, std = %0.2f', skewness_x_mean, skewness_x_std));
            
            subplot(nrows, ncols, 3)
            hx_skew = histogram(skewness_x);
            hold on;
            hx_skew_bin_centers = 0.5 * (hx_skew.BinEdges(2:end) + hx_skew.BinEdges(1:end-1));
            init_params_x_skew = [mean(hx_skew_bin_centers), std(hx_skew_bin_centers), max(hx_skew.Values), 0];
            fixed_params_x_skew = [0,0,0,1];
            % fit histogram to gaussian
            [fit_params_x_skew, ~, ffh_x_skew, std_errs_x_skew, chi_sqr_x_skew] = ...
                            fit1D(hx_skew_bin_centers, hx_skew.Values,...
                            [], {'gaussian1D'}, init_params_x_skew, fixed_params_x_skew);
            x_skew_interp = linspace(min(hx_skew.BinEdges), max(hx_skew.BinEdges), 300);
            plot(x_skew_interp, ffh_x_skew(x_skew_interp));
            title(sprintf('histogram x skewnesss\n gaussian fit mean = %0.2f, sigma = %0.2f',...
                fit_params_x_skew(1), fit_params_x_skew(2)));
            
                        
            subplot(nrows, ncols, 2)
            plot(1:n_imgs, skewness_y);
            xlabel('picture');
            ylabel('skewness, y');
            % ylim(-3, 3]);
            title(sprintf('y skewness \n mean = %0.2f, std = %0.2f', skewness_y_mean, skewness_y_std));
            
            subplot(nrows, ncols, 4)
            hy_skew = histogram(skewness_y);
            hold on;
            hy_skew_bin_centers = 0.5 * (hy_skew.BinEdges(2:end) + hy_skew.BinEdges(1:end-1));
            init_params_y_skew = [mean(hy_skew_bin_centers), std(hy_skew_bin_centers), max(hy_skew.Values), 0];
            fixed_params_y_skew = [0,0,0,1];
            % fit histogram to gaussian
            [fit_params_y_skew, ~, ffh_y_skew, std_errs_y_skew, chi_sqr_y_skew] = ...
                            fit1D(hy_skew_bin_centers, hy_skew.Values,...
                            [], {'gaussian1D'}, init_params_y_skew, fixed_params_y_skew);
            y_skew_interp = linspace(min(hy_skew.BinEdges), max(hy_skew.BinEdges), 300);
            plot(y_skew_interp, ffh_y_skew(y_skew_interp));
            title(sprintf('histogram y skewnesss\n gaussian fit mean = %0.2f, sigma = %0.2f',...
                fit_params_y_skew(1), fit_params_y_skew(2)));
            
            suptitle(sprintf('%s \n cloud skewness distribution', obj.DatasetString));
        end
        
        function [fig_handle] = showNNCorr(obj, Indices, save_results, save_dir, file_name)
            %showNNCorr Display correlator average azimuthally and over images.   
            %
            %   showNNCorr(obj,Indices,Axes) Displays the density
            %   correlator C(Indices(1),Indices(2)) versus bin. If the
            %   optional argument Axes as a handle to an axes object,
            %   displays the results on these axes.
            %   
            %   TODO: In next version of DataFolder, rename function to
            %   reflect the fact that this can display any correlator, not
            %   only the nearest neighbor.
            
            % TODO: remove saving options ...
            if ~exist('Indices','var')
                Indices = [0,1];
            end
            
            if ~exist('save_results','var')
                save_results = 0;
            end
            
            if ~exist('save_dir','var') || isempty(save_dir)
                save_dir = pwd;
            end
            
            if ~exist('file_name','var')
                file_name = '';
            end
            

            fig_handle = figure('name','NNCorr');
            d4indices = get_d4_indices(Indices);
            LegendEntries = {'Symm Corr'};
            colormap(jet(size(d4indices, 1)+1));
            cmap = colormap;
            %symmetric correlator
            subplot(2,2,1)
            errorbar(obj.BinAvg,...
                squeeze(obj.Density_Corr_AzAvg(obj.CenterIndex_CorrMatrix+Indices(1),...
                                               obj.CenterIndex_CorrMatrix+Indices(2),:)),...
                squeeze(obj.Density_Corr_AzAvgUnc(obj.CenterIndex_CorrMatrix+Indices(1),...
                                                  obj.CenterIndex_CorrMatrix+Indices(2),:)),...
                'color',cmap(1,:));
            hold on;
            %other correlators
            for ii = 1:size(d4indices,1)
                CurrentIndices = d4indices(ii,:);
                CorrVert = squeeze(obj.Density_Corr_AzAvg_NotSymmetrized(obj.CenterIndex_CorrMatrix + CurrentIndices(1), obj.CenterIndex_CorrMatrix+CurrentIndices(2),:));
                CorrVertUnc = squeeze(obj.Density_Corr_AzAvgUnc_NotSymmetrized(obj.CenterIndex_CorrMatrix+CurrentIndices(1),obj.CenterIndex_CorrMatrix+CurrentIndices(2),:));
%                 errorbar(Axes,obj.BinAvg,CorrVert,CorrVertUnc,'color',cmap(ii+1,:));
                errorbar(obj.BinAvg,CorrVert,CorrVertUnc,'color',cmap(ii+1,:));
                LegendEntries = cat(2,LegendEntries,sprintf('C(%d,%d)',CurrentIndices(1),CurrentIndices(2)));
            end
            hold off;
            grid on;
            xlabel('Position (Latt Sites');
            ylabel('Correlator');
            legend(LegendEntries);
            title('correlator vs. bin position');
            
            ax = gca;
            ax.Color = [0.5, 0.5, 0.5];
            
            subplot(2,2,2)
            errorbar(obj.RadialPos,...
                squeeze(obj.Density_Corr_AzAvg(obj.CenterIndex_CorrMatrix+Indices(1),...
                                               obj.CenterIndex_CorrMatrix+Indices(2),:)),...
                squeeze(obj.Density_Corr_AzAvgUnc(obj.CenterIndex_CorrMatrix+Indices(1),...
                                                  obj.CenterIndex_CorrMatrix+Indices(2),:)),...
                'color',cmap(1,:));
            hold on;
             for ii = 1:size(d4indices,1)
                CurrentIndices = d4indices(ii,:);
                CorrVert = squeeze(obj.Density_Corr_AzAvg_NotSymmetrized(obj.CenterIndex_CorrMatrix + CurrentIndices(1), obj.CenterIndex_CorrMatrix+CurrentIndices(2),:));
                CorrVertUnc = squeeze(obj.Density_Corr_AzAvgUnc_NotSymmetrized(obj.CenterIndex_CorrMatrix+CurrentIndices(1),obj.CenterIndex_CorrMatrix+CurrentIndices(2),:));
%                 errorbar(Axes,obj.BinAvg,CorrVert,CorrVertUnc,'color',cmap(ii+1,:));
                errorbar(obj.RadialPos, CorrVert, CorrVertUnc, 'color', cmap(ii+1,:));
                LegendEntries = cat(2, LegendEntries, sprintf('C(%d,%d)', CurrentIndices(1),CurrentIndices(2)));
            end
            hold off;
            grid on;
            xlabel('Position (Latt Sites');
            ylabel('Correlator');
            legend(LegendEntries);
            title('correlator vs. spatial position');
            
            ax = gca;
            ax.Color = [0.5, 0.5, 0.5];
            
            subplot(2,2,3)
            imagesc(squeeze(obj.nnc2D(:, :, obj.CenterIndex_CorrMatrix + d4indices(1,1),...
                obj.CenterIndex_CorrMatrix + d4indices(1,2))),[-0.05,0.05]);
            if isempty(obj.ColorMap)
                obj.getColorMap();
            end
            colormap(obj.ColorMap);
            axis equal; axis image;
            colorbar;
            title(sprintf('C(%d,%d)',d4indices(1,1),d4indices(1,2)));
            
            subplot(2,2,4)
            imagesc(squeeze(obj.nnc2D(:,:,obj.CenterIndex_CorrMatrix + d4indices(2,1),obj.CenterIndex_CorrMatrix + d4indices(2,2))),[-0.05,0.05]);
            colormap(obj.ColorMap);
            axis equal; axis image;
            colorbar;
            title(sprintf('C(%d,%d)',d4indices(2,1),d4indices(2,2)));
            
            ttl = sprintf('%s\n  NPics = %d',obj.identifier,size(obj.Occs_Stack,3));
            suptitle(ttl);
            
            if save_results
                if isempty(file_name)
                    file_name = sprintf('%s_correlator_dx=%d_dy=%d.fig',obj.identifier,d4indices(1,1),d4indices(1,2));
                end
                fpath = fullfile(save_dir,file_name);
                savefig(fig_handle,fpath);
                fprintf('Saved to %s\n',fpath);
            end
        end
        
        function showNNCorrImbalance(obj, Axes)
            %showNNCorrImbalance
            %
            %    showNNCorrImbalance(obj,Axes)
            %   Display (particle-particle) correlator
            %   C(Indices(1),Indices(2)) vs. azimuthal average bin.
            
            if ~exist('Indices','var')
                Indices = [0,1];
            end
            NAtomsList = sum(sum(obj.Occs_Stack,1),2);
            FigName = sprintf('<dd>_c, NPics = %d, NAtoms = %0.1f +/- %0.1f',size(obj.Occs_Stack,3),mean(NAtomsList),std(NAtomsList));
            
            if ~exist('Axes','var')
                h = figure('name',FigName);
                Axes = axes;
            end
            CenterIndex = obj.NumNeighbors + 1;

            %symmetric correlator
            CorrV1 = squeeze(obj.Density_Corr_AzAvg_NotSymmetrized(CenterIndex+1,CenterIndex+0,:));
            CorrV1Unc = squeeze(obj.Density_Corr_AzAvgUnc_NotSymmetrized(CenterIndex+1,CenterIndex+0,:));
            CorrV2 = squeeze(obj.Density_Corr_AzAvg_NotSymmetrized(CenterIndex-1,CenterIndex,:));
            CorrV2Unc = squeeze(obj.Density_Corr_AzAvgUnc_NotSymmetrized(CenterIndex-1,CenterIndex,:));
            CorrVmean = 0.5*(CorrV1+CorrV2);
            CorrVmeanUnc = 0.5*sqrt(CorrV1Unc.^2+CorrV2Unc.^2);
            
            CorrH1 = squeeze(obj.Density_Corr_AzAvg_NotSymmetrized(CenterIndex+0,CenterIndex+1,:));
            CorrH1Unc = squeeze(obj.Density_Corr_AzAvgUnc_NotSymmetrized(CenterIndex+0,CenterIndex+1,:));
            CorrH2 = squeeze(obj.Density_Corr_AzAvg_NotSymmetrized(CenterIndex+0,CenterIndex-1,:));
            CorrH2Unc = squeeze(obj.Density_Corr_AzAvgUnc_NotSymmetrized(CenterIndex+0,CenterIndex-1,:));
            CorrHmean = 0.5*(CorrH1+CorrH2);
            CorrHmeanUnc = 0.5*sqrt(CorrH1Unc.^2+CorrH2Unc.^2);
            
            CorrSymm = squeeze(obj.Density_Corr_AzAvg(CenterIndex+0,CenterIndex+1,:));
            CorrSymmUnc = squeeze(obj.Density_Corr_AzAvgUnc(CenterIndex+0,CenterIndex+1,:));
            
            subplot(3,1,1)        
            errorbar(obj.BinAvg,CorrSymm,CorrSymmUnc,'b');
            hold on;
            errorbar(obj.BinAvg,CorrVmean,CorrVmeanUnc);
            errorbar(obj.BinAvg,CorrHmean,CorrHmeanUnc);
            legend({'Symmetrized','V','H'});
            grid on;
            
            subplot(3,1,2)
            plot(obj.BinAvg,CorrVmean-CorrHmean,'bo-');
            ylabel('Corr Diff')
            grid on;
            
            subplot(3,1,3)
            plot(obj.BinAvg,(CorrVmean-CorrHmean)./abs(CorrSymm),'bo-');
            ylabel('Corr Diff/Mean')
            grid on;
            ylim([-1,1])
            hold off;
            grid on;

            suptitle(obj.identifier);
        end
        
        function [AllCorrLen, AllCorrLenUnc] = ...
                showCorrelationRange(obj, BinIndices, CorrMatStack, CorrMatUncStack)           
            %showCorrelationRange Show range of density correlator.
            %   
            %   showCorrelationRange(obj,BinIndex)
            %   Plot correlators versus separation for a given azimuthal
            %   average bin specifie by BinIndex.
            
            if ~exist('BinIndices','var')
                BinIndices = 1;
            end
            
            if ~exist('CorrMatStack','var')
                CorrMatStack = obj.Density_Corr_AzAvg;
                CorrMatUncStack = obj.Density_Corr_AzAvgUnc;
            end
            
            if ~exist('CorrMatUncStack','var')
                CorrMatUncStack = obj.Density_Corr_AzAvgUnc;
            end
            
            NNeighbors = (size(CorrMatStack,1) - 1) / 2;
            CenterIndex = NNeighbors + 1;
            QuadrantCorrMatrix = CorrMatStack(CenterIndex:end, CenterIndex:end, BinIndices);
            QuadrantCorrMatrixUnc = CorrMatUncStack(CenterIndex:end, CenterIndex:end, BinIndices);
            [X, Y] = meshgrid(0:NNeighbors, 0:NNeighbors);
            R = sqrt(X.^2 + Y.^2);
            
            if ~exist('Axes','var')
                FigName = sprintf('Correlation Range');
                figure('name',FigName);
                Axes = axes;
            end
            
            AllCorrLen = zeros([1,length(BinIndices)]);
            AllCorrLenUnc = zeros([1,length(BinIndices)]);
            leg = cell([1,length(BinIndices)]);
            plothandles = [];
            
            %fit range
            InitP = [0, 0.05, 1, 0];
            FixedP = [1, 0, 0, 1];
            RsNoZero = R(R ~= 0);
            InterpRs = linspace(min(RsNoZero), max(RsNoZero), 300);
            
            colormap(jet(length(BinIndices)));
            cmap = colormap;
            for ii = 1:length(BinIndices)
                CurrentQuadCorrMat = QuadrantCorrMatrix(:,:,ii);
                CurrentQuadCorrMatUnc = QuadrantCorrMatrixUnc(:,:,ii);
                AbsCorrNoZero = abs(CurrentQuadCorrMat(R~=0));
                [Fp,~,FFH,SE] = fit1D(RsNoZero,AbsCorrNoZero,[],{'exp1D'},InitP,FixedP);
                CorrLen = Fp(3);
                Err = SE(3);
                
                AllCorrLen(ii) = CorrLen;
                AllCorrLenUnc(ii) = Err;
                leg{ii} = sprintf('Bin %0.1f-%0.1f, l = %0.2f+/-%0.2f',...
                    obj.BinEdges(BinIndices(ii)),obj.BinEdges(BinIndices(ii)+1),CorrLen,Err);

                ph = errorbar(Axes,R(:),abs(CurrentQuadCorrMat(:)),CurrentQuadCorrMatUnc(:),'o','color',cmap(ii,:));
                plothandles = cat(1,plothandles,ph);
                hold on;
                plot(Axes,InterpRs,FFH(InterpRs),'color',cmap(ii,:));
            end
            
            grid on;
            xlabel('Correlator Distance (Lattice Sites)')
            ylabel('|C(r)|')
            legend(plothandles,leg);
            title(obj.identifier);
            
            ax = gca;
            ax.Color = [0.5, 0.5, 0.5];
        end
        
        function [fig_handle] = showCorrelationMatrix(obj, BinIndex)
            %showCorrelationMatrix
            %
            %   showCorrelationMatrix(obj,BinIndex)
            %   Display correlation matrix for a given bin. Also display
            %   structure factor.
            %
            %   TODO: Unify showStructureFactor with this function
            
            obj.getColorMap;
            Lims = [-0.05,0.05];
            if ~exist('BinIndex','var')
                BinIndex = 1;
            end
            
            fig_handle = figure('name','Correlation Matrix');
            NRows = 2;
            NCols = 3;
            
            % unsymmetrized correlator
            subplot(NRows,NCols,1)
            imagesc(obj.Density_Corr_AzAvg_NotSymmetrized(:, :, BinIndex), Lims)
            title('Not symmetrized')
            colorbar;
            axis equal; axis image;
            
            % unsymmetrized correlator / unc
            subplot(NRows, NCols, NCols + 1)
            corr_over_unc_notsymm =...
                abs(obj.Density_Corr_AzAvg_NotSymmetrized(:, :, BinIndex)) ./ ...
                abs(obj.Density_Corr_AzAvgUnc_NotSymmetrized(:, :, BinIndex));      
            imagesc(corr_over_unc_notsymm, [-5,5]);
            title('Corr/Unc Not Symmetrized');
            colorbar;
            axis equal; 
            axis image;
            
            % symmetrized correlator
            subplot(NRows, NCols, 2)
            imagesc(obj.Density_Corr_AzAvg(:, :, BinIndex), Lims);
            title('Symmetrized')
            colorbar;
            axis equal; 
            axis image;
            
            % symmetrized correlator /unc
            subplot(NRows, NCols, NCols + 2)
            corr_over_unc_symm = ...
            abs(obj.Density_Corr_AzAvg(:, :, BinIndex)) ./ ...
            (obj.Density_Corr_AzAvgUnc(:, :, BinIndex));
            imagesc(corr_over_unc_symm, [-5, 5]);
            title('Corr/Unc Symmetrized');
            colorbar;
            axis equal; 
            axis image;           
            
            %%%Show structure factor.
            StructFactStack = obj.Static_Structure_Factor;
            StructFactUncStack = obj.Static_Structure_FactorUnc;
            subplot(NRows,NCols,3);
            imagesc(real(StructFactStack(:,:,BinIndex)), [-1,1]);
            ax = gca;
            ax.XTick = 1:2:size(obj.Qxs,2);
            ax.YTick = 1:2:size(obj.Qys,1);
            ax.XTickLabel = round(obj.Qxs(1,1:2:end), 1);
            ax.YTickLabel = round(obj.Qys(1:2:end,1), 1);
            xlabel('qx')
            ylabel('qy')
            title('Structure Factor');
            colorbar;
            axis equal; 
            axis image;
            
            % plot structure factor along cuts in k-space
            subplot(NRows, NCols, NCols + 3);
            N = (size(StructFactStack,1) - 1) / 2;
            %cut from (0,0)->(pi,0)
            Cut1 = abs(transpose(StructFactStack(N + 1 ,N + 1:end, BinIndex)));
            Cut1_Unc = abs(transpose(StructFactUncStack(N + 1, N + 1:end, BinIndex)));
            %(pi,0)->(pi,pi)
            Cut2 = abs(StructFactStack(N + 2:end, end, BinIndex));
            Cut2_Unc = abs(StructFactUncStack(N + 2:end, end, BinIndex));
            %(pi,pi)->(0,0)
            Cut3 = flip(abs(diag(StructFactStack(N + 1:end - 1, N + 1:end - 1, BinIndex))));
            Cut3_Unc = flip(abs(diag(StructFactUncStack(N + 1:end - 1, N + 1:end - 1, BinIndex))));
            errorbar([Cut1; Cut2; Cut3], [Cut1_Unc; Cut2_Unc; Cut3_Unc]);
            
            try
               hold on;
               StrFac_UncBootStr = obj.Static_Structure_Factor_BootstrapUnc;
               Cut1_UncBootStr = abs(transpose(StrFac_UncBootStr(N + 1, N + 1:end, BinIndex)));
               Cut2_UncBootStr = abs(StrFac_UncBootStr(N+2:end,end,BinIndex));
               Cut3_UncBootStr = flip(abs(diag(StrFac_UncBootStr(N + 1:end - 1, N + 1:end - 1, BinIndex))));
               errorbar([Cut1; Cut2; Cut3], [Cut1_UncBootStr; Cut2_UncBootStr; Cut3_UncBootStr]); 
            catch
            end
                
            grid on;            
            ax = gca;
            ax.XTick = [1, N + 1, 2 * N + 1, 3 * N + 1];
            ax.XTickLabel = {'(0,0)', '(0,pi)', '(pi,pi)', '(0,0)'};
            ax.YLim(1) = min([0, ax.YLim(1)]);
            title('Structure Factor');            
            
            colormap(obj.ColorMap)
            ttl = sprintf('%s\n Bin %d, Limits = %0.1f-%0.1f, NPics = %d',...
                obj.identifier, BinIndex,...
                obj.BinEdges(BinIndex), obj.BinEdges(BinIndex+1),...
                size(obj.Occs_Stack, 3));
            suptitle(ttl);
        end
        
        function [fig_handle] = show2DCorrelator(obj, indices)
            %show2DCorrelator
            %
            %   show2DCorrelator(obj,Indices) Given a pair of indices, show
            %   all correlators C(i,j) with equivalent coordinates. Display
            %   2D correlator average over images.
            %   
            %   e.g. if Indices = [0,1] this correspondings to the
            %   correlator C(0,1) which should be identical to C(0,-1) and
            %   is related by symmetry to C(1,0) = C(-1,0).
            %
            %   Each correlator will either have 4 symmetric correlators,
            %   as in the example above, or 8. The correlator C(2,1) is
            %   related to 8 correlators. C(2,-1), C(-2,1), C(-2,-1) and
            %   C(1,2), C(-1,2), C(1,-2), C(-1,-2).
            
             if ~exist('indices', 'var')
                indices = [0,1];
             end
            
            CenterIndex = obj.NumNeighbors + 1;

            D4Indices = get_d4_indices(indices);
           
            %other correlators
            if size(D4Indices,1)==4
                NRows = 2; 
                NCols = 2;
            else
                NRows = 3; 
                NCols = 3;
            end
            
            fig_handle = figure;
            for ii = 1:size(D4Indices,1)
                subplot(NRows,NCols,ii)
                Corr = squeeze(obj.nnc2D(:,:,CenterIndex+D4Indices(ii,1),CenterIndex+D4Indices(ii,2)));
                imagesc(Corr,[-0.1,0.1]);
                axis equal; axis image;
                title(sprintf('C(%d,%d)',D4Indices(ii,1),D4Indices(ii,2)));
            end
            if isempty(obj.ColorMap)
                obj.getColorMap();
            end
            colormap(obj.ColorMap);

            ttl = sprintf('%s\n  NPics = %d',obj.identifier,size(obj.Occs_Stack,3));
            suptitle(ttl);
        end
        
        function [fig_handle] = showG2(obj, indices)
            %showG2
            %
            %   showG2(obj,Indices) Display g2 correlator.
            %
            %   TODO: Add showing all D4 symmetric correlators.
            %   TODO: Add showing 2D g2 correlator.
            
            if ~exist('indices','var')
                indices = [0, 1];
            end
            
            fig_name = sprintf('g2');
            fig_handle = figure('name', fig_name);
            
            CenterIndex = obj.NumNeighbors + 1;
            g2 = squeeze(obj.g2_FullCorr_AzAvg(CenterIndex+indices(1), CenterIndex+indices(2), :));
            g2unc = squeeze(obj.g2_FullCorr_AzAvgUnc(CenterIndex + indices(1), CenterIndex + indices(2), :));
            errorbar(obj.BinAvg, g2, g2unc);
            grid on;
            xlabel('Position (Lattice Sites)');
            ylabel(sprintf('g_2(%d,%d)', indices(1), indices(2)));
            ylim([0, 1.5]);
            
            ttl = sprintf('%s\n NPics = %d, NAtoms = %0.1f +/- %0.1f\n g_2(i,j) = <n_i n_j> / (<n_i><n_j>)',...
                obj.identifier, size(obj.Occs_Stack,3), obj.MeanAtomNum, obj.AtomNumSD);
            suptitle(ttl);
        end
        
        function [fig_handle_array] = showAllPics(obj, picture_numbers)
            %showAllPics
            %
            %   showAllPics(obj) Display all image in image stack used in
            %   data analysis.
            if ~exist('picture_numbers', 'var')
                picture_numbers = 1 : size(obj.Occs_Stack, 3);
            end
            
            indices = 1 : size(obj.Occs_Stack, 3);
            indices = indices( ismember(obj.PictureNumberInFolder, picture_numbers) );
            
            if isempty(indices)
                error('None of the picture number arguments passed to showAllPics were present in this instance.');
            end
            
            fig_handle_array = [];
            for ii = indices
                ttl = sprintf('Pic %d, %s \n ROI Start X=%d, Y=%d, Xcom = %0.2f, Ycom = %0.2f',...
                    obj.PictureNumberInFolder(ii), obj.identifier,...
                    obj.XStarts_ROI(ii), obj.YStarts_ROI(ii),...
                    obj.X_CoM_fullpic(ii), obj.Y_CoM_fullpic(ii));
                fig_handle = figure;
%                 subplot(1, 2, 1)
                imagesc(obj.Occs_Stack(:, :, ii), [0, 1]);
                hold on;
                % center of 2d gaussian fit
                scatter(obj.GaussFitParams(1), obj.GaussFitParams(2), 'ro');
                % azimuthal average center
                scatter(obj.Cx_AzAvg, obj.Cy_AzAvg, 'gx');
                % show center of mass in cropped picture
                scatter(obj.X_CoM_fullpic(ii) - obj.XStarts_ROI(ii) + 1,...
                    obj.Y_CoM_fullpic(ii) - obj.YStarts_ROI(ii) + 1, 200, 'r.');
                legend({'gaussian fit to avg', 'az avg center', 'pic com'});
                axis equal; 
                axis image;
                xlim([1, obj.ImgCropSize]);
                ylim([1, obj.ImgCropSize]);
                
%                 subplot(1, 2, 2);
                % TODO: show fluorescence image
                suptitle(ttl);
                
                fig_handle_array = horzcat(fig_handle_array, fig_handle);
            end
        end
        
        function [ind_var, dens, dens_unc, fig_handle] = showBinVsIndVar(obj, bin_index)
            %showBinVsIndVar(obj,BinIndex)
            %   this function intended for the situation where the folder is a
            %   collection of files with parameters that are changing. That is
            %   a very different assumption than most of the functions here,
            %   which treat a folder as if all the files in it use exactly the
            %   same parameters.
            %   TODO: get rid of this function. Funcionality subsumed by
            %   showAllVsIndVar
            if ~exist('BinIndex','var')
                bin_index = 1;
            end
            [ind_var, dens, dens_unc, ~, ~, fig_handle] = obj.showAllVsIndVar(bin_index);
        end
        
        function [ind_var_unique, all_dens_unique, all_dens_unc_unique,...
                 num_points, std_dev, FigHandle] = ...
                showAllVsIndVar(obj, bin_indices, AtomNumFlag)
            %showAllVsIndVar
            %
            %   showAllVsIndVar Displays the azimuthally averaged density
            %   in a selection of bins versus the independant variable
            %   defined by obj.IndependentVariable
            
            % input checking
            if ~exist('bin_indices', 'var') || isempty(bin_indices)
                bin_indices = 1 : obj.NBins;
            end
            
            if ~exist('AtomNumFlag', 'var')
                AtomNumFlag = 0;
            end
            
            % now implented through the more general get_versus_ind_var
            % function
            AllDens = transpose( obj.Occs_AzAvgStack(bin_indices, :) );
            AllDensUnc = transpose( obj.Occs_AzAvgStackUnc(bin_indices, :) );
            if AtomNumFlag
                AllDensUnc = zeros(size(AllDensUnc));
            end
            [ind_var_unique, all_dens_unique, all_dens_unc_unique,...
                 num_points, std_dev, repetition_num] = obj.get_versus_ind_var(AllDens, AllDensUnc );

             ind_var = obj.IndependentVariable;
             
             if size(AllDens, 1) ~= length(ind_var)
                 warning('AllDens and ind_var are different lengths');
                 if size(AllDens, 1) > length(ind_var)
                     AllDens = AllDens(1 : length(ind_var), :);
                 elseif length(ind_var) > size(AllDens, 1)
                     ind_var = ind_var( 1 : size(AllDens, 1) );
                 end
             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot results
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Leg = cell([1,length(bin_indices)]);
            for ii = 1:length(Leg)
                Leg{ii} = sprintf('Bin %0.1f-%0.1f',...
                    obj.BinEdges(bin_indices(ii)), obj.BinEdges(bin_indices(ii) + 1));
            end
            
            ttl = sprintf('%s \n All Bins Vs. Ind Var', obj.identifier);
            FigHandle = figure('name', 'Ind Var');                
                
            max_repetitions = max(num_points(:));
            colors = jet( max_repetitions );

            subplot(3, 1, 1);
            hold on;
            
            for ii = 1:length(bin_indices)
                
                % if only one bin, plot different repetitions of the
                % independent variable as different colors
                if size(AllDens, 2) == 1
                    for jj = 1 : max_repetitions
                        plot(ind_var( repetition_num(:, ii) == jj ) , AllDens( repetition_num(:, ii) == jj , ii),...
                             '-x', 'Color', colors(jj, :) );
                        hold on;
                    end
                % if multiple bins, plot different bins as different colors
                else
                    plot(ind_var, AllDens(:, ii), 'x');
                end
            end
            
            grid on;
            legend(Leg);
            ylabel('Filling');
            ax = gca;
            ax.YLim(1) = 0;
%                 ylim([0, 1]);
            title('All points');

            subplot(3, 1, 2);
            hold on;
            for ii = 1:length(bin_indices)
                errorbar(ind_var_unique, all_dens_unique(:,ii), all_dens_unc_unique(:,ii), '.-');
            end
            grid on;
            legend(Leg);
            ylabel('Filling');
            xlabel('Ind Var');
            ax = gca;
            ax.YLim(1) = 0;
            title('Averaged points');

            subplot(3, 1, 3);
            hold on;
            plot(ind_var_unique, num_points)
            grid on;
            ylim([0, max(num_points(:)) + 1]);
            ylabel('Num Points');
            xlabel('Ind Var');
            title('number of points averaged');

            suptitle(ttl);
            
        end
        
        function [ind_var_unique, array_unique, array_unc_unique, num_points,...
                  std_dev, repetition_num] = get_versus_ind_var(obj, array, array_unc)
            % replacement for showAllVsIndVar which allows you to
            % use any quantity, instead of only density. The only
            % assumption is that the first dimension of array is the same
            % size as the number of pictures in the folder.
            %
            % array: an array of arbitrary number of dimensions. The first
            % dimension must be the same size as obj.nimgs
            %
            % array_unc: an array the same size as array containing the
            % uncertainty in each element of array.
            ind_var = obj.IndependentVariable;
            if ~exist('array_unc', 'var') || isempty(array_unc)
                array_unc = zeros(size(array));
            end
            
            if ~isequal( size(array), size(array_unc))
                error('array and array_unc must be the same size');
            end
                     
            if isempty(ind_var)
                ind_var = obj.PictureNumberInFolder;
            end
            
            if length(ind_var) > obj.nimgs
                  ind_var = ind_var(obj.PictureNumberInFolder(obj.PictureNumberInFolder <= length(ind_var)));
            end
            
            if isempty(ind_var)
                ind_var = 1 : obj.nimgs;
            end
            
            % sort and round to avoid numerical precision issues with
            % duplicates
            [ind_var, I] = sort(ind_var);
            ind_var = round(ind_var, 10);
            
            % if more images than independent variables, ignore images
            % after last independent variables
            % TODO: crop array of arbitrary dimension along first dimension
            array_size = size(array);
            array = array(1:length(ind_var), :);
            array = reshape(array, [length(ind_var), array_size(2:end)]);
            
            array_unc = array_unc(1:length(ind_var), :);
            array_unc = reshape(array_unc, [length(ind_var), array_size(2:end)]);
            
            % initialize arrays for storing results
            % sort array 
            array(:) = array(I, :);
            array_unc(:) = array_unc(I, :);
            % repetition number is not bin-dependent
            repetition_num = zeros( size(array, 1), 1 );
            
            % determine duplicates
            [ind_var_unique, iu, ia] = unique(ind_var);
            % the unique function works such that
            % ind_var(iu) = ind_var_unique
            % ind_var_unique(ia) = ind_var
            
            % determine repetition number for each shot. This is useful if
            % we are taking many repetitions of the same points.
            for jj = 1 : length(ind_var_unique)
                repetition_num(ind_var == ind_var_unique(jj), 1) = ...
                      1 : sum( ind_var == ind_var_unique(jj) ); 
            end
            
            % average duplicates
            % ia is the index into ind_var_unique, so if two elements of
            % ind_var are the same they will have an identical index in ia
            % number of points at independent variable value
            num_points = accumarray(ia, 1);
            
            % compute average and standard deviation of those points
            remaining_size = size( array );
            remaining_size = remaining_size(2:end);
            accumulator_size = horzcat( length(ind_var_unique), remaining_size);
            
            std_dev = zeros( accumulator_size );
            array_unique = zeros( accumulator_size );
            array_unc_unique = zeros( accumulator_size );
            % loop over bins
            for ii = 1 : numel( array(1, :) )
                array_unique(:, ii) = accumarray(ia, array(:, ii), [], @mean);
                std_dev(:, ii) = accumarray(ia, array(:, ii), [], @std);
                sum_sqr_unc = accumarray(ia, array_unc(:, ii).^2, [], @mean);

                % other option would be to compute std dev of mean from values...
                array_unc_unique(:, ii) = sqrt(sum_sqr_unc) ./ sqrt(num_points);
                if sum(sum_sqr_unc) == 0
                    array_unc_unique(:, ii) = std_dev(:,ii) ./ sqrt(num_points);
                end
            end
 
        end
        
        function [ind_var_all, dens_avg, dens_sdm, num_pts] = ...
                getCombinedSetsVsIndVar(obj, data_sets_stack, bin_indices, display)
            %getCombinedSetsVsIndVar
            %
            %   [IndVarAll,DensAvg,DensSdm,NumPts] = 
            %   getCombinedSetsVsIndVar(obj,DataSets,BinIndices,Display)
            %   Given a stack of DataFolder class objects, DataSets, this
            %   function combines the density in a selection of bins and
            %   the independent variables in the folders. These are
            %   returned in IndVarAll, DensAvg, and DenSdm. NumPts gives
            %   the number of data points averaged for each value of
            %   IndVarAll. This function assumes that the same azimuthal
            %   average bins are being used for all DataFolder instances.
            %
            %   TODO add support for variables besides density...e.g.
            %   correlators, etc.
            %   
            %   TODO rewrite this in a way that uses accumarray as
            %   showAllVsIndVar
            
            if ~exist('bin_indices', 'var') || isempty(bin_indices)
                bin_indices = 1 : obj.NBins;
            end
            
            if ~exist('display','var') || isempty(display)
                display = 0;
            end
            
            ind_var_all = [];
            DensAll = [];
            DensUncAll = [];
            
            ExpandDim = 3;
            
            data_sets_stack = cat(2, obj, data_sets_stack);
            for ii = 1:length(data_sets_stack)
                dfTemp = data_sets_stack(ii);
                IndVarTemp = round(dfTemp.IndependentVariable, 10);
                
                if isempty(IndVarTemp)
                    warning('%s IndependentVariables was empty\n', dfTemp.DatasetString);
                end
                
                if length(IndVarTemp)>size(dfTemp.Occs_AzAvgStack, 2)
                      IndVarTemp = IndVarTemp(dfTemp.PictureNumberInFolder);
                end

                [IndVarTemp,I] = sort(IndVarTemp);
                DensTemp  = transpose(dfTemp.Occs_AzAvgStack(bin_indices, I));
                DensUncTemp = transpose(dfTemp.Occs_AzAvgStackUnc(bin_indices, I));
                
                if isempty(ind_var_all) && isempty(DensAll) && isempty(DensUncAll)
                    ind_var_all = IndVarTemp;
                    DensAll = DensTemp;
                    DensUncAll = DensUncTemp;
                else
                    % does this average duplicates within one set?
                
                    %Now need to add these to our sets...but first need to make
                    %sure they both have the same independent variable points,
                    %or if they don't, to expand them until they do.
                    [IndVar_TempNotFull, I] = setdiff(IndVarTemp, ind_var_all); 
                    [IndVar_FullNotTemp, J] = setdiff(ind_var_all, IndVarTemp);

                    %first expand the full variables, if necessary
                    ind_var_all = cat(2, ind_var_all, IndVar_TempNotFull);
                    NanExtender = nan([length(IndVar_TempNotFull), size(DensAll, 2), size(DensAll, 3)]);
                    DensAll = cat(1, DensAll, NanExtender);
                    DensUncAll = cat(1, DensUncAll, NanExtender);

                    [ind_var_all, IAll] = sort(ind_var_all);
                    DensAll = DensAll(IAll, :, :);
                    DensUncAll = DensUncAll(IAll, :, :);

                    %now expand the temporary variables
                    IndVarTemp = cat(2, IndVarTemp, IndVar_FullNotTemp);
                    NanExtender = nan([ length(IndVar_FullNotTemp), size(DensTemp, 2), size(DensTemp, 3) ]);
                    DensTemp = cat(1, DensTemp, NanExtender);
                    DensUncTemp = cat(1, DensUncTemp, NanExtender);

                    [IndVarTemp, ITemp] = sort(IndVarTemp);
                    DensTemp = DensTemp(ITemp, :);
                    DensUncTemp = DensUncTemp(ITemp, :);

                    %add temporary vars as next layer of full vars
                    DensAll = cat(ExpandDim, DensAll, DensTemp);
                    DensUncAll = cat(ExpandDim, DensUncAll, DensUncTemp);
                end 
            end
            num_pts = sum(~isnan(DensAll), 3);
            dens_avg = nansum(DensAll, 3) ./ num_pts;
            DensSqAvg = nansum(DensAll.^2, 3) ./ num_pts;
           
            NumPtsMinusOne = num_pts - 1;
            NumPtsMinusOne(NumPtsMinusOne == 0) = 1;
            dens_sdm = sqrt(num_pts ./ NumPtsMinusOne) .* sqrt(DensSqAvg-dens_avg.^2) ./ sqrt(num_pts);
            num_pts = num_pts(:, 1);
            
            if display
                Leg = {};
                figure;
                hold on;
                for ii = 1:size(dens_avg,2)
                    errorbar(ind_var_all, dens_avg(:,ii), dens_sdm(:,ii), '.-');
                    Leg{ii} = sprintf('Bin %0.1f-%0.1f',...
                        obj.BinEdges(bin_indices(ii)), obj.BinEdges(bin_indices(ii) + 1));
                end
                grid on;
                legend(Leg);
                title(sprintf('%s\n', data_sets_stack.identifier));
            end
        end
        
        function [fig_handle] = compareSets(obj, obj_stack, corr_indices)
            %compareSets Visual comparison of a stack of DataFolder class
            %objects.
            %
            %   FigHandle = compareSets(obj,objStack,CorrIndices) Accepts a
            %   stack of DataFolder class instances, objStack, and displays
            %   density and correlators for all of the sets. CorrIndices is
            %   a 1x2 vector which specifices the correlator as
            %   C(CorrIndices(1),CorrIndices(2)). Returns a handle to the
            %   figure, FigHandle. This function does not require that the
            %   various DataFolder class objects use the same azimuthal
            %   average binning.
            if ~exist('obj_stack', 'var') || isempty(obj_stack)
                % if you have a stack of objects, simply call
                % obj_stack.compareSets()
                obj_stack = obj;
            else
                % preserve old syntax, which was obj.compareSets(obj_stack)
                % or if you have a stack, it was the awkward idiom
                % obj_stack(1).compareSets(obj_stack(2:end));
                obj_stack = vertcat(obj, obj_stack);
            end
            
            if ~exist('corr_indices', 'var') || isempty(corr_indices)
                corr_indices = [0,1];
            end
            
            n_sets = length(obj_stack);
            cmap = jet(n_sets);
%             cmap = get_color_map(10, [1, 0, 0], [0, 0, 1], [1, 1, 1]);
            
            % create title and legend information
            ttl = sprintf('');
            Leg = cell([1, n_sets]); 
            for ii = 1 : n_sets
                obj2 = obj_stack(ii);
                ttl = sprintf('%s%s\n', ttl, obj2.identifier);
                Leg{ii} = obj2.identifier;
            end
            
            fig_handle = figure('name', 'Compare Sets');
            nrows = 2;
            ncols = 2;
            
            % plot azimuthal average
            subplot(nrows, ncols, 1)
            
            hold on;
            for ii = 1 : n_sets
                obj2 = obj_stack(ii);
                errorbar(obj2.RadialPos, obj2.Occs_AzAvg, obj2.Occs_AzAvgUnc,...
                    '-o', 'color', cmap(ii, :));
            end

            ylim([0,1])
            grid on;
            xlabel('Lattice Sites')
            ylabel('Filling')
            title('Density Vs. Avg Rad Pos');
            ax = gca;
            ax.Color = [0.5, 0.5, 0.5];
                
            % plot azimuthal average, slightly different method of
            % determining radius position
            subplot(nrows, ncols, 2)
            hold on;
            for ii = 1 : n_sets
                obj2 = obj_stack(ii);
                errorbar(obj2.BinAvg, obj2.Occs_AzAvg,...
                    obj2.Occs_AzAvgUnc, '-o', 'color', cmap(ii, :));
            end
            ylim([0,1])
            grid on;
            xlabel('Lattice Sites')
            ylabel('Filling')
            title('Density');
            title('Density Vs. Avg Bin Pos')
            ax = gca;
            ax.Color = [0.5, 0.5, 0.5];
            
            % plot correlators
            subplot(nrows, ncols, 3)
               
            hold on;
            for ii = 1 : n_sets
                obj2 = obj_stack(ii);
                NNCorr = squeeze(obj2.Density_Corr_AzAvg(obj2.CenterIndex_CorrMatrix + corr_indices(1),...
                    obj2.CenterIndex_CorrMatrix + corr_indices(2), :));
                NNCorrUnc = squeeze(obj2.Density_Corr_AzAvgUnc(obj2.CenterIndex_CorrMatrix + corr_indices(1),...
                    obj2.CenterIndex_CorrMatrix + corr_indices(2), :));

                errorbar(obj2.RadialPos, NNCorr, NNCorrUnc, '-o', 'color' ,cmap(ii, :));
            end
            grid on;
            legend(Leg);
            xlabel('Lattice Sites')
            ylabel(sprintf('Correlator (%d,%d)', corr_indices(1), corr_indices(2)))
            title('Correlators');
            ax = gca;
            ax.Color = [0.5, 0.5, 0.5];
                
            % plot atom numbers
            subplot(nrows, ncols, 4)
            hold on;
            for ii = 1 : n_sets
                obj2 = obj_stack(ii);
                plot(obj2.PictureNumberInFolder, obj2.AtomNumbers, 'o', 'color', cmap(ii, :));
            end
            grid on;
            xlabel('Shot')
            ylabel('Atom Number')
            ax = gca;
            ax.XLim(1) = 1;
            ax.YLim(1) = 0;
            ax = gca;
            ax.Color = [0.5, 0.5, 0.5];
            
        end
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Helper functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %philosophy with these is that they shouldn't reference class
        %member directly. Take them as arguments and return values. So in
        %main class initialization function would call these with class
        %members as arguments and save the returned values as other class
        %members. These could in principle be stored as separate functions.
        %Just keep them in the class so I can save a copy of this one file
        %and be (mostly) able to run. So e.g. whenever I run analysis, I
        %should save a copy of the current version of this class with the
        %analysis so it can be re-rerun later without worrying about later
        %modifications.
        
        function [DisplayMask_Img] = showMasks(obj, mask_stack)
            % showMasks Given a stack of disjoint binary arrays, returns an
            % array of numbers specifying which mask a given point is in.
            %
            %   [DisplayMask_Img] = showMasks(MaskStack) Given MaskStack,
            %   an Ny x Nx x NImgs array where M(i,j) = 1 for at most one
            %   of these arrays. Return a 2D array DisplayMask_Img where
            %   DisplayMask_Img(i,j) gives the index along the third
            %   dimension of MaskStack such that MaskStack(i,j,k) = 1.
            
            
            DisplayMask_Img = zeros(size(mask_stack,1),size(mask_stack,2));

            NMasks = size(mask_stack,3);
            for ii = 1:NMasks
                CurrMask = mask_stack(:,:,ii);
                CurrMask(CurrMask==1) = ii;
                DisplayMask_Img = DisplayMask_Img + CurrMask;
            end
        end
        
        function dataset = get_FullDataset(obj, fileIndex)
            %
            %   dataset = get_FullDataset(obj,fileIndex)
            %   
            %   Arguments:
            %
            %   fileIndex: = {year, month, day, folder#} gives the indices
            %   for a certain dataset. Adding extra entries beyond the
            %   first four will not affect the functionality.
            %
            %   Returns:
            %
            %   dataset: = is a nfiles x 6 cell array, where each entry
            %   dataset{ii,:} = {year, month, day, folder#, file#,
            %   picture#}
           
            if ~iscell(fileIndex{1})
                fileIndex={fileIndex};
            end

            dataset = {};
            for k = 1:size(fileIndex, 2)

                PathID = cell2mat(fileIndex{k}(1:5));
                AnExampleFile = findFileGeneral(PathID, '*.mat');
                [Path, ~, ~] = fileparts(AnExampleFile);
                FullSearchString = fullfile(Path, obj.RecFileFormat);
                Files = dir(FullSearchString);
                NumFiles = length(Files);

                for i = 1 : NumFiles
                    dataset(end + 1, :) =  {PathID(1) PathID(2) PathID(3) PathID(4) i 2};
                end
            end
        end
        
        function [filenames] = getReconstructedFilenames(obj, dataset)
            % example input:
            % dataset = {2016 7 21 4 12 3:5;
            %            2016 7 21 4 13 3};
               filenames = cell(size(dataset, 1), 1);
               for i=1:size(dataset,1)
                  [filename] = findFileGeneral(cell2mat(dataset(i,1:5)), '*reconstr.mat');
                  filenames{i} = filename;
               end
        end
        
        function [SymmCorrMat, SymmCorrErrMat] = getSymmetricCorrMat(obj, corr_mat, corr_mat_err)
            %[SymmCorrMat,SymmCorrErrMat] = getSymmetricCorrMat(CorrMat,CorrErrMat)
            %average over D4 symmetric matrix...    

            % TODO: replaced with d4_average function and test that gives
            % same results (with better errorbar estimation)
            %should capture full D4 symmetry...
            NDims = ndims(corr_mat);
            Permutation = [2, 1, 3 : NDims];

            Rotations = corr_mat + rot90(corr_mat, 1) + rot90(corr_mat,2) + rot90(corr_mat,3);
            CorrTranspose = permute(corr_mat, Permutation);
            TransposeRotations = CorrTranspose + rot90(CorrTranspose,1) ...
                               + rot90(CorrTranspose,2) + rot90(CorrTranspose,3);
            SymmCorrMat = (1/8) * (Rotations + TransposeRotations);

            %current feeling is that the rotations represent different things...but the
            %transpose the same thing...
            NDims = ndims(corr_mat_err);
            Permutation = [2, 1, 3 : NDims];
            CorrErrTranspose = permute(corr_mat_err, Permutation);
            %use all errors to get a better idea of the error
            ErrType1 = 0.25 * (corr_mat_err + rot90(corr_mat_err,2) + ...
                               CorrErrTranspose + rot90(CorrErrTranspose,2));
            ErrType2 = 0.25 * (rot90(corr_mat_err,1) + rot90(corr_mat_err,3) ...
                             + rot90(CorrErrTranspose,1) + rot90(CorrErrTranspose,3));
            %but only fair to reduce the error by non-equivalent directions...
            SymmCorrErrMat = 0.5 * sqrt(ErrType1 .^ 2 + ErrType2 .^ 2);
        end
        
        function CharFn_Prod = getCorrPtsInSameBin(obj, MaskStack, NumNeighbors)
            %getCorrPtsInSameBin Creates a matrix describing if two
            %coordinates are in the same bin or not.
            %
            %   SameBinMat = getCorrPtsInSameBin(obj,MaskStack,NumNeighbors)
            %   Mask stack is an Ny x Nx x NImgs stack of binary arrays,
            %   where each slice represents a bin. So each slice has the
            %   value 1 at points in the bin, and 0 otherwise. The various
            %   slices are disjoint. Given an integer NumNeighbors number
            %   of neighbors to consider, SameBinMat is Ny x Nx x NImgs x
            %   2*NumNeighbors+1 x 2*NumNeighbors+1.
            %   $$SameBinMat(y,x,m,i_y,i_x) =
            %   \sum_b Char_b(y,x)*Char_b(y-(i_y - NumNeighbors-1),
            %   x-(i_x-NumNeighbors-1)),$$ where $Char_b(y,x)$ is
            %   the characteristic function of bin b.
            %   
            %   TODO: think about integrating this with function that create
            %   correlation matrix?
            
            % TODO: make this a function in library instead of in class
            
            %CorrObj is Ny x Nx x NCorr x NCorr
            %Mask Stack is Ny x Nx x NBins
            %Shift mask sum and find points which are equal...
            BinNums = obj.showMasks(MaskStack);
            [ShiftIndexMatrix_X,ShiftIndexMatrix_Y] = meshgrid(-NumNeighbors : 1 : NumNeighbors);
            
            FullSize = 2*NumNeighbors + 1;
            %CharFn_Prod(y,x,)
            
            %can do this
            CharFn_Prod = zeros(size(BinNums, 1), size(BinNums, 2), FullSize, FullSize);
            for jj = 1:FullSize
                for kk = 1:FullSize
                    ShiftedBinNums = getShiftedMat(BinNums, ShiftIndexMatrix_X(jj, kk), ShiftIndexMatrix_Y(jj, kk), 0);
%                     ShiftedBinNums = obj.getShiftedMat(BinNums,ShiftIndexMatrix_X(jj,kk),ShiftIndexMatrix_Y(jj,kk),0);
%                     ShiftedBinNums = circshift(BinNums,[ShiftIndexMatrix_Y(jj,kk),ShiftIndexMatrix_X(jj,kk)]);
                    CharFn_Prod(:, :, jj, kk) = (BinNums == ShiftedBinNums);          
                end
            end    
            
            %or maybe this is nicer...
%             ShiftedBinNums = obj.getAllShiftedMats(BinNums,NumNeighbors);
%             ExpandedBinNums = repmat(BinNums,[1,1,FullSize,FullSize]);
%             CharFn_Prod = (ShiftedBinNums == ExpandedBinNums);
        end
        
        function getColorMap(obj)
            %create colormap to use for displaying correlators. Typically
            %use red for negative, white for zero, blue for positive.
            StartColor = [0,0,1];
            EndColor = [1,0,0];
            CenterColor = [1,1,1];
            obj.ColorMap = get_color_map(200, StartColor, CenterColor, EndColor);
        end
        
        function [cloud_density_fractions] = get_density_fraction(obj, n_bin_edges)
            % determine the fraction of our cloud between given densities.
            
            % remove zero density regions
            n_azavg = obj.Occs_AzAvg;
            r_edges = obj.BinEdges;
            
            r_starts = r_edges(1:end - 1);
            r_ends = r_edges(2:end);
            
            % if we have a bunch of zero density points in our cloud, we
            % want to ignore these.
            index_first_zero = find( n_azavg == 0);
            if ~isempty(index_first_zero)
                index_first_zero = index_first_zero(1);

                r_starts = r_starts( 1 : index_first_zero - 1 );
                r_ends = r_ends( 1 : index_first_zero - 1 );
                n_azavg = n_azavg( 1 : index_first_zero - 1 );
            end
            
            % fraction of cloud at each area
            area_fraction = transpose( r_ends .^2 - r_starts .^2) / r_ends(end).^2;

            [~, ~, area_fraction_binned, ~, ~, npts, ~] = azAvg_General(...
                area_fraction, [], n_azavg, n_bin_edges);
                                                   
            cloud_density_fractions = area_fraction_binned .* npts;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % coordinate transformations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [img_rec_in_fl_space, Xfl, Yfl] = get_avg_reconstructed_img_fl_space(obj)
            
            nimgs = obj.nimgs; %size(obj.Occs_Stack, 3);
            
            % fluorescence picture coordinates
            n_border_pix = 5;
            % TODO: where to store these???
            [Xfl,Yfl] = meshgrid(1:801, 1:801);
            rec_in_fl_space_stack = zeros(size(Xfl, 1), size(Xfl, 2), nimgs);

            for ii = 1 : nimgs
                %transform info
                TrformStrct = obj.ReconstrInfoStack(ii).coordTrafo;
                GridParams = TrformStrct.gridParameters;

                %Some thing where cutoff border? Possibly why my transform not
                %working?

                theta1 = GridParams.theta1;
                theta2 = GridParams.theta2;
                lambda1 = GridParams.lambda(1);
                lambda2 = GridParams.lambda(2);
                phi1 = GridParams.phi1;
                phi2 = GridParams.phi2;
                offset1 = TrformStrct.indexOffset(1);
                offset2 = TrformStrct.indexOffset(2);
                borderPx = TrformStrct.borderPx;

                %cropped reconstructed picture
                RecPic = obj.Occs_Stack(:,:,ii); 
                [x_crop_rec, y_crop_rec] = obj.get_uncropped_rec_coords(ii);

                %transformed images.
                xform_params = [theta1, theta2, phi1, phi2, lambda1, lambda2,...
                                offset1, offset2, n_border_pix, n_border_pix];

                % reconstructed image in fluorescence space
                [Xr_from_fl, Yr_from_fl, ~] = img2latt_coord(xform_params, Xfl, Yfl);
                % Define a transformation T(Xl, Yl) = (Xi, Yi)
                % Define our function f_l(Xl, Yl) = M[floor(Yl), floor(Xl)], with M the
                % reconstructon matrix.
                % Then in image space, we have f_i(Xi, Yi) = f_l(T^{-1)(Xi, yi) = ...
                x_start = obj.XStarts_ROI(ii);
                y_start = obj.YStarts_ROI(ii);

                % get indices in cropped picture
                y_indices = round(Yr_from_fl) - y_start + 1;   
                x_indices = round(Xr_from_fl) - x_start + 1;

                mask = ones(size(y_indices));
                mask(y_indices <= 0) = 0;
                mask(y_indices > obj.ImgCropSize) = 0;
                mask(x_indices <= 0) = 0;
                mask(x_indices > obj.ImgCropSize) = 0;

                % for the moment, pick an index that is likely to be zero
                y_indices(~mask) = 1;
                x_indices(~mask) = 1;

                xform_pic = RecPic(sub2ind(size(RecPic), y_indices, x_indices));
                xform_pic( ~mask ) = 0;
                
                rec_in_fl_space_stack(:, :, ii) = xform_pic;
                
            %     figure;
            %     imagesc(rec_in_fl_space_stack(:, :, ii));
            %     axis equal;
            %     axis image;
            end

            img_rec_in_fl_space = mean(rec_in_fl_space_stack, 3);
%             fig_handle = figure;
%             imagesc(img_rec_in_fl_space);
%             axis equal;
%             axis image;
        end
        
        function fighandle = TestTransforms(obj, Index)
            %Test transformation between lattice coordinates and real space
            %coordinates. This appears to be working after correcting for
            %border cropping of fluorescence image by reconstruction.
            if ~exist('Index','var')
                Index = 1;
            end
            
            %border pixels
            %reconstruction crops the fluoresence image before processing.
            %We must account for this cropping in our transformation.
            n_border_pix = 5;
            
            %transform info
            TrformStrct = obj.ReconstrInfoStack(Index).coordTrafo;
            GridParams = TrformStrct.gridParameters;
            
            %Some thing where cutoff border? Possibly why my transform not
            %working?
            
            theta1 = GridParams.theta1;
            theta2 = GridParams.theta2;
            lambda1 = GridParams.lambda(1);
            lambda2 = GridParams.lambda(2);
            phi1 = GridParams.phi1; %think needs a 2*pi???
            phi2 = GridParams.phi2;
            sintheta1 = GridParams.sintheta1;
            costheta1 = GridParams.costheta1;
            sintheta2 = GridParams.sintheta2;
            costheta2 = GridParams.costheta2;
            transformFactor = GridParams.transformFactor;
            offsetxy = GridParams.offsetxy;
            
            offset1 = TrformStrct.indexOffset(1);
            offset2 = TrformStrct.indexOffset(2);
            borderPx = TrformStrct.borderPx;
            
            %full reconstructed picture
            a = load(obj.ReconFPaths{Index});
            FullRecPic = double(a.occupationsRounded>0);
            [Xfrec,Yfrec] = meshgrid(1:size(FullRecPic,2), 1:size(FullRecPic,1));
            FullRecInterp = @(Xl,Yl) interp2(Xfrec, Yfrec, FullRecPic, Xl, Yl);
            
            %cropped reconstructed picture
            RecPic = obj.Occs_Stack(:,:,Index); 
            [Xcrec,Ycrec] = meshgrid(1:size(RecPic, 2), 1:size(RecPic, 1));
            CropRecInterp = @(Xl,Yl) interp2(Xcrec, Ycrec, RecPic, Xl, Yl);
                 
            %Fluorescence picture
            FlPic = readimg(char(obj.FlFPaths{Index}));
            FlPic = FlPic(:,:,2);
%             FlPic = FlPic(n_border_pix + 1:end-n_border_pix, n_border_pix + 1:end-n_border_pix);
            [Xfl,Yfl] = meshgrid(1:size(FlPic,2),1:size(FlPic,1));
            FlInterp = @(Xi,Yi) interp2(Xfl,Yfl,FlPic,Xi,Yi);
                        
            %transformed images.
            xform_params = [theta1, theta2, phi1, phi2, lambda1, lambda2, offset1, offset2, n_border_pix, n_border_pix];
            %Fluorescence image in reconstruction space
            [Xfl_from_rec, Yfl_from_rec, ~] = latt2img_coord(xform_params, Xfrec, Yfrec);
            FlImg_RecSpace = FlInterp(Xfl_from_rec,Yfl_from_rec);
            %Fluorescence image in cropped reconstruction space
            Xfrec_from_crec = Xcrec + obj.XStarts_ROI(Index) - 1;
            Yfrec_from_crec = Ycrec + obj.YStarts_ROI(Index) - 1;
            [Xfl_from_cropped_rec, Yfl_from_cropped_rec, ~] = ....
                latt2img_coord(xform_params, Xfrec_from_crec, Yfrec_from_crec);
            FlImg_CroppedRecSpace = FlInterp(Xfl_from_cropped_rec, Yfl_from_cropped_rec);
            
            %reconstructed image in fluorescence space
            [Xr_from_fl, Yr_from_fl, ~] = img2latt_coord(xform_params, Xfl, Yfl);
            Rec_FlSpace = FullRecInterp(Xr_from_fl, Yr_from_fl);
            
            fighandle = figure;
            %Fluorescence image space
            subplot(3, 3, 1)
            imagesc(FlPic, [120, 300]);
            axis equal; 
            axis image;
            title('Full Fl Imag');
            
            subplot(3, 3, 4)
            imagesc(Rec_FlSpace, [0, 1]);
            axis equal; 
            axis image;
            title('Rec Img, Fl Space');
            
            %Reconstruction space
            subplot(3, 3, 2)
            imagesc(FullRecPic, [0, 1]);
            axis equal; 
            axis image;
            title('Full Rec Img')
            
            subplot(3, 3, 5)
            imagesc(FlImg_RecSpace, [120, 300]);
            axis equal; 
            axis image;
            title('Fl Img, Rec Space');
            
            %Reconstruction space after cropping
            subplot(3,3,3)
            imagesc(RecPic,[0,1]);
            axis equal; axis image
            title('Cropped Rec Img');
            
            subplot(3, 3, 6)
            XStart = obj.XStarts_ROI(Index);
            YStart = obj.YStarts_ROI(Index);
%             FullImgROI = obj.getROI(FullRecPic,XStart,obj.ImgCropSize,YStart,obj.ImgCropSize);
            FullImgROI = getROI(FullRecPic, XStart, obj.ImgCropSize, YStart, obj.ImgCropSize);
            imagesc(RecPic - FullImgROI, [0, 1]);
            axis equal; 
            axis image;
            title(sprintf('CropRecImg - Full \n Equal = %d',...
                isequal(FullImgROI,RecPic)));
            
            subplot(3, 3, 9)
            imagesc(FlImg_CroppedRecSpace, [120, 300]);
            axis equal; 
            axis image;
            title('Fl Img, Crop Rec Space')
            
            suptitle(sprintf('%s\n image %03d', obj.identifier, Index));
            
        end
        
        function [x_coords, y_coords] = get_uncropped_rec_coords(obj, index)
            
            x_start = obj.XStarts_ROI(index);
            y_start = obj.YStarts_ROI(index);
            [x_coords, y_coords] = meshgrid(x_start : x_start + obj.ImgCropSize - 1,...
                                              y_start : y_start + obj.ImgCropSize - 1);

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % comparison of data set objects
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [sets_equal] = compare_dsets(obj, other_dset)
            %compare two DataFolder instances to see if they have the same
            %properties. Only consider properties that exist in the first
            %set. One case where this is useful: you want to compare a
            %loaded set with an instantiated but un-analyzed set. The
            %un-analyzed set only has properties related to analysis
            %settings, but no data. So you will only compare analysis
            %settings. If these are all the same, then you do not need to
            %re-analyze the data.
            fields1 = fieldnames(obj);
            sets_equal = 1;
            ii = 1;
            while sets_equal && ii <= length(fields1)
                compare_1 = obj.(fields1{ii});
                if ~isempty(compare_1)
                    compare_2 = other_dset.(fields1{ii});
                    if isnumeric(compare_1)
                        %remove any nans, because isequal(nan,nan) = 0
                        compare_1 = compare_1(~isnan(compare_1)); 
                        compare_2 = compare_2(~isnan(compare_2));
                    end
                    if ~isequal(compare_1,compare_2)
                        sets_equal = 0;
                    end
                end
                ii = ii + 1;
            end
                    
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % saving and loading routines
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
        function [DataFolderAsStruct] = getStruct(obj)
            % get object as struct. Useful for classes containing instances
            % of DataFolder.
            DataFolderAsStruct = struct();
            fields = fieldnames(obj);
            for ii = 1:length(fields)
                DataFolderAsStruct.(fields{ii}) = obj.(fields{ii});
            end
        end
        
        function saveStruct(obj, save_dir, fname)
            % saveStruct(obj, save_dir, fname)
            % Save class to a structure. This has the advantage that it can
            % be loaded even if you don't have the class code at hand.
            %
            % save_dir
            %
            % fname
            
            if ~exist('save_dir', 'var')
                save_dir = pwd;
            end
            
            if ~exist('fname', 'var') || isempty(fname)
                if ~isempty(obj.identifier)
                    name = obj.identifier;
                else
                    name = 'test';
                end
                fname = fullfile([name,'.mat']);
            end
            
            % make sure file type is correct
            [~,~,ext] = fileparts(fname);
            if isempty(ext)
                fname = fullfile([fname,'.mat']);
            elseif ~strcmp(ext,'.mat')
                error('Extension must be .mat');
            end
            
            if ~exist('save_dir', 'var') || isempty(save_dir)
                save_dir = pwd;
            end
      
            DataFolderAsStruct = obj.getStruct();
            
            fpath = fullfile( save_dir, fname);
            save(fpath,'DataFolderAsStruct');  
        end
        
        function loadStruct(obj, StructToLoad)
            % load a struct and turn it into a class object.
            
            if ischar(StructToLoad)
                a = load(StructToLoad);
                StructToLoad = a.DataFolderAsStruct;
            end

            if ~isstruct(StructToLoad)
                error('StructToLoad was not a structure');
            end
            
            %clear out the fields of our obj first.
            ObjFields = fieldnames(obj);
            for ii = 1:length(ObjFields)
                obj.(ObjFields{ii}) = [];
            end
            
            % loop over structure fields and load those which are fields of
            % DataFolder class.
            StructFields = fieldnames(StructToLoad);
            for ii = 1:length(StructFields)
                try
                    obj.(StructFields{ii}) = StructToLoad.(StructFields{ii});
                catch
                    fprintf('field %s was present in loaded struct, but it is not a field of DataFolder class. Skipped this field. \n',...
                        StructFields{ii});
                end
            end
        end
                
    end
end

