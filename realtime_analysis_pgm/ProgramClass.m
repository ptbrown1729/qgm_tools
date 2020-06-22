classdef ProgramClass < handle
    
    properties
        % top level directory. or imaging files
        RootImageDirectory = '';
        
        % directory storing analysis code, which is the directory of this file.
        % Do not hard code because we want to use local copies with svn...
        RootAnalysisPath = fileparts(mfilename('fullpath')); 
		
        % other paths with required functions
        lib_relative_path = '../lib';
        code_resource_dirs = {};
        
        % image file extensions
        Extensions = {'.aia', '.fits'};
        
        % name of settings files
        % TODO: could change this a pattern
        SettingsFileName = 'settings.m';
        
        % name of log files
        LogFileName = 'log.txt';
        
        % class storing fundamental constants
        Constants
        
        % folder data needed for processing.
        RootFolder = ''; %read subfolders of this folder, if ReadSubfolders = 1.
        % TODO: think it might be better to create a separate class which
        % would store all of these things. Then ProgramClass would simply
        % store a stack of these. This should be easier to manage than a
        % set of many lists which need to be kept synchronized.
        
        % list of all folders to process, as relative paths to obj.RootImageDirectory
        FolderNames
        % list of all folders data is being saved in...basically FolderNames, but fully qualified.
        SavePaths = {};
        % list of Settings class objects to be used in processing respective folders.
        SettingsList = [];
        % list of paths image files are moved to and processed in as our function reads them.
        ProcessPaths = {};
        % list of log file paths to save data from processing folders.
        LogPaths = {}; 
        % list of function handles to 'Experiment' functions.
        ExperimentFnHandlers = {};
        % Last data set processed for each folder.
        MarkerDataSetList = [];
        % Index of folder currently being processed.
        CurrentFolderIndex = 1;
        
        %resources for readfolder
        continue_proc_loop = 1;
        ReadSubFolders = 0;
        FigHandle
        FigName = 'Real Time Analysis';
    end
    methods
        
        function obj = ProgramClass(FolderNames, ReadSubFolders, Constants, RootImageDirectory)
            % Function to initialize a program class object.
            %
            % Arguments
            % ----------------------
            % FolderNames: a cell array of folders to be processed,
            % relative to obj.RootImageDirectory
            %
            % ReadSubFolders: a boolean parameter. If FALSE, only the
            %listed Folders are processed. If TRUE, any subfolders of the
            %given folder are processed.
            %
            % Constants: a ConstantsClass object containing physical data
            % about the given atom. If no argument is given, default to
            % lithium-6
            %
            % RootImageDirectory: Root directory for relative file paths.
            % If no argument is given defaults to
            % //128.112.86.75/lithium/IMAGING DATA
            
            %Currently unclear what happens if you give a cell array of
            %folder names and set ReadSubFolders=TRUE...TODO...figure that
            %out...
            
            if ~exist('RootImageDirectory', 'var') || isempty(RootImageDirectory)
                RootImageDirectory = fullfile('//', '128.112.86.75', 'lithium', 'Imaging Data'); 
            end
            obj.RootImageDirectory = RootImageDirectory;
            
            if ~exist('Constants', 'var') || isempty(Constants)
                Constants = ConstantsClass();
            end
            obj.Constants = Constants;
            
            if ~exist('FolderNames', 'var') || isempty(FolderNames)
                %default behavior. Run todays datapath...
                %maybe hard coding this isn't the best way?
                FolderNames = obj.get_today_folder_path();
            end
            
            if ~exist('ReadSubFolders', 'var') || isempty(ReadSubFolders)
                ReadSubFolders = 1;
            end
            
            obj.run_program(FolderNames, ReadSubFolders)
        end
        
        function run_program(obj, folder_names, read_subfolders)
            % Function for starting analysis loop
            %
            % Arguments:
            % -------------------
            % FolderNames: a cell array of folders to be processed,
            % relative to obj.RootImageDirectory
            %
            % ReadSubFolders: a boolean parameter. If FALSE, only the
            % listed Folders are processed. If TRUE, any subfolders of the
            % given folder are processed.

            obj.checkVersion;
            obj.add_dependencies_to_path;
            
            % ensure FolderNames is a cell array
            if ischar(folder_names)
                obj.FolderNames = {folder_names};
            else
                obj.FolderNames = folder_names;
            end
          
            obj.ReadSubFolders = read_subfolders;
     
            if obj.ReadSubFolders
                obj.RootFolder = folder_names;
                obj.FolderNames = obj.getAllSubFolders(folder_names);
            end
      
            obj.generateFolderData;
            obj.getImageMarkers;
            obj.readfolder; %main function
        end
        
        function folder_name = get_today_folder_path(obj)
            % Get `today' data directory relative to root directory
            
            folder_name = fullfile(datestr(today, 'yyyy'), datestr(today, 'mm'), datestr(today, 'dd'));
        end
        
        function HandleOut = createFig(obj, FigName)
            % create a figure that will hold all the analysis data and
            % include controls for stopping execution, etc.
            %
            
            fhandle = findobj('type', 'figure', 'name', FigName);
            if ~isempty(fhandle)
                HandleOut = fhandle(1);
                if length(fhandle > 1)
                    close(fhandle(2:end));
                end
            else
                HandleOut = figure;
                HandleOut.Name = FigName;
            end
            
            %add some useful controls
            
            % button to stop analysis
            stop_callback = @(hObject, callbackdata) obj.setContinueOff;
            c1 = uicontrol('Style', 'pushbutton', 'String', 'Stop',...
                'Position', [20, 20, 60, 20], 'Callback', stop_callback);
            
            % button to reprocess last image
            reprocess_previous_callback = @(hobject, callbackdata) obj.reprocessLastImage;
            c2 = uicontrol('Style', 'pushbutton', 'String',...
                'Reprocess Previous', 'Position', [100, 20, 120, 20],...
                'Callback', reprocess_previous_callback);
            
            % button to reprocess folder
            reprocess_folder_callback = @(hobject, callbackdata) obj.reprocessCurrentFolder;
            c3 = uicontrol('Style', 'pushbutton', 'String',...
                'Reprocess Folder', 'Position', [240, 20, 100, 20],...
                'Callback', reprocess_folder_callback);
            
            % button to reprocess everything
            reprocess_all_callback = @(hobject, callbackdata) obj.reprocessAll;
            c4 = uicontrol('Style', 'pushbutton', 'String',...
                'Reprocess All', 'Position', [360, 20, 100, 20],...
                'Callback', reprocess_all_callback);
        end
        
        function setContinueOff(obj)
            obj.continue_proc_loop = 0;
        end
        
        function setContinueOn(obj)
            obj.continue_proc_loop = 1;
        end
        
        function checkVersion(obj)
            %checks matlab version. Raises warning if your version is too early.
            versionDate = version('-date');
            if datenum(versionDate) < datenum(2015, 07, 1)
                warning('Parts of this program were written using matlab 2015b. Some functions do not run on matlab 2014 and earlier.');
            end
        end
        
        function add_dependencies_to_path(obj)
		% add subdirectories necessary for this program to work correctly.
            
            lib_path = fullfile(obj.RootAnalysisPath, obj.lib_relative_path);
            obj.code_resource_dirs = strsplit(genpath(lib_path), ';');   

			% Folders on path.
			path_cell = regexp(path, pathsep, 'split');

			%Check if directory is already on path, because adding to path takes a long
			%time.
			for ii = 1 : length(obj.code_resource_dirs)
% 				current_path = fullfile(obj.RootAnalysisPath, obj.sub_dirs{ii});
                current_path = obj.code_resource_dirs{ii};
                
				if ~any(strcmpi(current_path, path_cell))
					addpath(current_path);
				end
			end
        end
        
        function SubFolders = getAllSubFolders(obj, FolderName)
            %SubFolders = getAllSubFolders(obj, FolderName)
            %Get a cell array of strings listing all sub folders of a given
            %folder that appear to be data sets. Only returns folders that
            %have settings files and aren't hidden (i.e. don't begin with
            %dots)
            %TO DO...get subfolder for multiple FolderNames???
            
            if iscell(FolderName)
                FolderName = FolderName{1};
            end
            
            SubFolders = {};
            Folders = dir(fullfile(obj.RootImageDirectory, FolderName));
            for ii = 1:length(Folders)
                SettingsPath = fullfile(obj.RootImageDirectory, FolderName, Folders(ii).name, obj.SettingsFileName);
                
                %exclude files we don't want.
                if Folders(ii).isdir && ...
                   ~strcmp(Folders(ii).name(1),'.') && ...
                   exist(SettingsPath, 'file') && ...
                   length(Folders(ii).name)>2
                
                   %implement check for settings file. Otherwise don't add that file.
                   CurrentFilePath = fullfile(FolderName, Folders(ii).name);
                   SubFolders = {SubFolders{:}, CurrentFilePath};
                end
            end
            
        end
        
        function addFolder(obj, new_folder_path)
            % addFolder(obj, new_folder_path)
            %
            % add a new folder and associated data (i.e. paths, settings, markers, etc.)
            % to the class instance.
            
            if ischar(obj.FolderNames)
                obj.FolderNames = {obj.FolderNames};
            end
            
            try
                new_save_path = fullfile(obj.RootImageDirectory, new_folder_path);
                new_settings_path = fullfile(new_save_path, obj.SettingsFileName);
                new_settings = SettingsClass(new_settings_path);
                new_process_path = fullfile(obj.RootImageDirectory, new_folder_path);
                
                if ~isdir(new_process_path)
                    mkdir(new_process_path);
                end
                
                NewLogPath = fullfile(new_process_path, obj.LogFileName);
                NewExpFn = str2func(new_settings.Experiment);
                
                %marker
                if ~exist(NewLogPath, 'file')
                    Paths = '';
                else
                    [~, ~, Paths] = parseLog(NewLogPath);
                end
                
                if isempty(Paths)
                    DataSet = DatasetClass;
                else
                    MarkerPath = char(Paths(end));
                    DataSet = path2dataset(MarkerPath);
                    DataSet.Second = DataSet.Second + 1;
                end
                Success = 1;
            catch
                Success = 0;
            end
            
            if Success
                obj.FolderNames = cat(2, obj.FolderNames, new_folder_path);
                obj.SavePaths = cat(2, obj.SavePaths, new_save_path);
                obj.ProcessPaths = cat(2, obj.ProcessPaths, new_process_path);
                obj.SettingsList = cat(2, obj.SettingsList, new_settings);
                obj.LogPaths = cat(2, obj.LogPaths, NewLogPath);
                obj.ExperimentFnHandlers = {obj.ExperimentFnHandlers{:}, NewExpFn};
                obj.MarkerDataSetList = cat(2, obj.MarkerDataSetList, DataSet);
            end
            
        end
        
        function addNewSubFolders(obj)
            % addNewSubFolders(obj)
            %
            % add any subfolders of the RootFolder that aren't already
            % listed on obj.FolderNames. (i.e. subfolders that have been
            % created since the last time this function was called).
            
            if ~obj.ReadSubFolders
                return;
            end
            
            sub_folders = obj.getAllSubFolders(obj.RootFolder);
            
            for ff = 1:length(sub_folders)
                if 0 == sum(strncmp(sub_folders{ff}, obj.FolderNames, length(sub_folders{ff})))
                    obj.addFolder(sub_folders{ff});
                end
            end
            
        end
        
        function generateFolderData(obj)
            % generateFolderData(obj)
            %
            % From folder names already stored as properties of your class instance
            % , get a collection of paths, settings files,
            % etc. that we will need for processing the folders.
            
            obj.SavePaths = {};
            obj.SettingsList = [];
            obj.ProcessPaths = {};
            obj.LogPaths = {};
            obj.ExperimentFnHandlers = {};
            
            %make sure obj.FolderNames is cell, not char.
            if ischar(obj.FolderNames)
                obj.FolderNames = {obj.FolderNames};
            end
            
            for ff = 1:length(obj.FolderNames)
                
                % location image files are being saved.
                current_fname = obj.FolderNames{ff};
                current_save_path = fullfile(obj.RootImageDirectory,...
                    current_fname);
                obj.SavePaths = cat(2, obj.SavePaths, current_save_path);
                
                % settings files
                current_settings_path = fullfile(current_save_path, obj.SettingsFileName);
                current_settings = SettingsClass(current_settings_path); 
                obj.SettingsList = cat(2, obj.SettingsList, current_settings);
                
                % Create folder and log file for processed data.
                current_process_path = fullfile(obj.RootImageDirectory, current_fname);
                obj.ProcessPaths = cat(2, obj.ProcessPaths, current_process_path);
                if ~isdir(current_process_path)
                    mkdir(current_process_path);
                end
                
                current_log_path = fullfile(current_process_path, obj.LogFileName);
                obj.LogPaths = cat(2, obj.LogPaths, current_log_path);
                
                % Get function handle of analysis function. experiment variable in settings
                % file directs to the appropriate function.
                obj.ExperimentFnHandlers = {obj.ExperimentFnHandlers{:},...
                    str2func(current_settings.Experiment)};
            end
            
        end
        
        function getImageMarkers(obj)
            %getImageMarkers(obj)
            %
            % Get "marker" images from log files. The timestamp of the
            % marker image lets us decide which images to process. We only
            % process files that are later than the marker image.
            % possible change this to save an m file?
            
            % loop over folders and aggregate marker image data set objects
            obj.MarkerDataSetList = [];
            for ff = 1:length(obj.FolderNames)
                
                if ~exist(obj.LogPaths{ff}, 'file')
                    Paths = '';
                else
                    [~, ~, Paths] = parseLog(obj.LogPaths{ff});
                end
                
                if isempty(Paths)
                    %obj.MarkerImageList(ff) = 1;
                    DummyDataset = DatasetClass;
                    obj.MarkerDataSetList = cat(2, obj.MarkerDataSetList, DummyDataset);
                else
                    MarkerPath = char(Paths{end});
                    DataSet = path2dataset(MarkerPath);
                    DataSet.Second = DataSet.Second + 1;
                    obj.MarkerDataSetList = cat(2, obj.MarkerDataSetList, DataSet);
                end
            end
            
        end
        
        function reprocessLastImage(obj)
            % reprocessLastImage(obj)
            %
            % reprocess the most recently processed image. Continue the
            % processing loop after this.
            
            obj.setContinueOn;
            CurrentSecond = obj.MarkerDataSetList(obj.CurrentFolderIndex).Second;
            obj.MarkerDataSetList(obj.CurrentFolderIndex).Second = CurrentSecond - 1;
            obj.readfolder; 
        end
        
        function reprocessCurrentFolder(obj)
            % reprocessCurrentFolder(obj)
            %
            % reprocess the current folder. Continue processing loop after this.
            
            obj.setContinueOn;
            obj.MarkerDataSetList(obj.CurrentFolderIndex) = DatasetClass;
            obj.readfolder;
        end
        
        function reprocessAll(obj)
            % reprocessAll(obj)
            %
            % Reprocess everything, then continue processing loop.
            
            obj.setContinueOn;
            for ii = 1:length(obj.MarkerDataSetList)
                obj.MarkerDataSetList(ii) = DatasetClass;
            end
            obj.readfolder;
        end
          
		function [dest] = mv_file_to_proc_folder(obj, save_path, proc_path, fname)
			%move file from the path it was saved in (SavePaths) to processing path (ProcessPaths)
            
			max_attempts = 5;
            move_attempts = 0;
            Moved = 0;
			while ~Moved && move_attempts < max_attempts
				try
                    move_attempts = move_attempts + 1;
					%Move file to processed folder.
					src = fullfile(save_path, fname);
					dest = fullfile(proc_path, fname);
					if ~strcmp(src, dest)
						movefile(src, dest); %movefile fails if paths are the same.
					end
					
					%Rename fits files using date/time metadata
					%stored in the file.
					[~, ~, ext] = fileparts(fname);
					% if strcmp(FileName(end-4:end),'.fits')
					if strcmp(ext, '.fits')
						dest = timeStampFitsFile(dest);
						[~, fname, ext] = fileparts(dest);
						fname = fullfile([fname, ext]);
					end
					
					Moved = 1;
					
				catch Exception
					MsgString = getReport(Exception);
					fprintf('Caught exception in readfolder trying to move file %s\n', fname);
					disp(MsgString);
					pause(2);
					continue
				end
			end
        end
		
        function displayData(obj, vals, keywords)
            % display data in terminal
            
            display = cat(2, num2cell(vals)', keywords');
            disp(display);
        end
        
        function readfolder(obj)
            % readfolder(obj)
            %
            % the main processing function. Continuously loop through all
            % folders looking for images, and process these as we see them.
            
            fprintf('analysing folder %s \n', obj.RootFolder);
            obj.generateFolderData;
            obj.FigHandle = obj.createFig(obj.FigName);
            obj.setContinueOn;
            
            % main processing loop.
            while obj.continue_proc_loop
                tic;
                obj.addNewSubFolders(); %check for new subfolders
                
                % loop over folders
                for nfolder = 1:length(obj.FolderNames)
                    obj.CurrentFolderIndex = nfolder;
                    
                    if ~exist(obj.SavePaths{nfolder}, 'dir')
                        fprintf('%s does not exist\n', obj.SavePaths{nfolder});
                        continue;
                    end
                    
                    %read files with allowed extensions
                    files = [];
                    for ii = 1:length(obj.Extensions)
                        files_ext = dir(fullfile(obj.SavePaths{nfolder}, ['*', obj.Extensions{ii}]));
                        files = horzcat(files, files_ext);
                    end
                    
                    %If there are image files, process them.
                    if isempty(files)
                        continue;
                    end
                    
                    % loop over files
                    for ii = 1:length(files)
                        if ~obj.continue_proc_loop
                            continue
                        end
                        
                        fname = files(ii).name;
                        dest = obj.mv_file_to_proc_folder(obj.SavePaths{nfolder},...
                            obj.ProcessPaths{nfolder}, fname); 

                        current_data_set = path2dataset(dest);

                        if current_data_set < obj.MarkerDataSetList(nfolder)
                            %only process files after the marker image.
                            continue
                        end
                        
                        try
                            disp(fname);
                            %MAIN ANALYSIS PART OF FUNCTION
                            
                            %read image data
                            [images, ~, time_stamp, vals, keys] = readimg(dest);
                            
                            % store image data in class
                            CurrentDataClass = ImageDataClass;
                            CurrentDataClass.Vals = vals;
                            CurrentDataClass.Keys = keys;
                            CurrentDataClass.OriginalImages = images;
                            CurrentDataClass.TimeStamp = time_stamp;
                            CurrentDataClass.Dataset = path2dataset(dest);
                            CurrentDataClass.Settings = obj.SettingsList(nfolder);
                           
                            %analysis tasks
                            % point to the function which will do the image
                            % processing and analysis
                            CurrentExperimentFnHandler = obj.ExperimentFnHandlers{nfolder};
                            
                            % do image processing and analysis
                            [CurrentDataClass] = CurrentExperimentFnHandler(CurrentDataClass,...
                                CurrentDataClass.Settings, obj.Constants, []);
                            
                            %display text data
                            obj.displayData(CurrentDataClass.Vals, CurrentDataClass.Keys);
                            
                            % display image data
                            if ~isvalid(obj.FigHandle)
                                %make sure figure handle is valid
                                obj.FigHandle = obj.createFig(obj.FigName);
                            end
                       
                            showPlotsUpdate(CurrentDataClass,...
                                CurrentDataClass.Settings, obj.FigHandle);
                            drawnow;
                            
                            %file processed
                            fprintf('%s processed\n', fullfile(obj.FolderNames{nfolder}, fname));
                            
                            %write information to log file
                            writeLog(obj.LogPaths{nfolder},...
                                CurrentDataClass.Vals,...
                                CurrentDataClass.Keys,...
                                CurrentDataClass.TimeStamp);
                            
                        catch Exception
                            MsgString = getReport(Exception);
                            disp(MsgString);
                            %write information to log file
                            obj.LogPaths{nfolder} = obj.LogPaths{nfolder};
                            writeLog(obj.LogPaths{nfolder}, [0], {'error'}, time_stamp);
                        end
                        
                        obj.MarkerDataSetList(nfolder) = current_data_set; %increment image counter.
                        obj.MarkerDataSetList(nfolder).Second =...
                        obj.MarkerDataSetList(nfolder).Second + 1;
                    end
                    
                    pause(0.5);  
                end
                loopT = toc;
            end
        end
        
    end
    
end