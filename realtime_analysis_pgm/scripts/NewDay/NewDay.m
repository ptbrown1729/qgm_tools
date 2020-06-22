%Prepare system for new day, by creating imaging data folders, other data folder, and
%copying settings files into some of these.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define relevant folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parent directory for everything
% fullfile and other functions are nice because they ensure this can run on
% windows or linux machines
root_dir = fullfile('//', '128.112.86.75', 'lithium');
% paranet directory for all imaging data
image_dir = 'Imaging Data';
% parent directory for other processing data
other_data_dir = 'Other Data';
% data diagnostics files are stored in
diagnostics_files_dir = fullfile(root_dir, 'Diagnostics', 'Potential Positions');
% other files to copy with paths relative to root_dir
other_file_paths = {};

% Get path of this script. 
new_day_dir = fileparts(mfilename('fullpath'));
% settings file dir
settings_files_dir = fullfile(new_day_dir, 'settings_files');
% miscellaneous script files to copy
script_files_dir = fullfile(new_day_dir, 'utility_scripts');
% other files to copy with paths relative to new_day_path
other_script_file_paths = {fullfile('..', '..', '..', 'lattice', 'Fit_Lattice_Depth.m'),...
                           fullfile('..', '..', '..', 'lattice', 'LatticeBandData_Interpolation.txt')};

% define startup folders. Settings files will be automatically detected if
% they are stored in settings_files_dir
startup_folder_list = {'000_MOT_And_CMOT', '001_ODT', '002_PotentialPositions','003_PCake_Axis_Imaging',...
    '004_Hopping', '005_PCake_TrpFrq', '006_LattMod', '007_Fluorescence_AngleCal'};

                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create startup folders and corresponding settings files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Path for today's imaging data
Today = date;
% Today = datenum(2019, 04, 19);
today_img_subdir = fullfile(root_dir, image_dir,...
    datestr(Today,'yyyy'), datestr(Today, 'mm'), datestr(Today, 'dd'));

% Create individual startup folders (i.e. 000_CMOT..., etc.) and copy
% settings files
for ii = 1 : length(startup_folder_list)
    
    try
        % create folders
        current_dir = fullfile(today_img_subdir, startup_folder_list{ii});
        if ~isdir(current_dir)
            mkdir(current_dir)
            fprintf('Successfully created directory %s\n', current_dir);
        end

        % get folder number. Use regular expressions to match 3-digit number
        re = '(\d{3})*';
        out = regexp(startup_folder_list{ii}, re, 'match');
        folder_num_str = out{1};

        % copy settings files
        settings_files = dir(fullfile(settings_files_dir, sprintf('%s*_settings.m', folder_num_str)));
        if isempty(settings_files)
            continue;
        end
        
        if length(settings_files) > 1
            fprintf('warning, more than one settings file found with the same number as folder %s\n', startup_folder_list{ii});
        end

        src_path = fullfile(settings_files.folder, settings_files.name);
        dest_path = fullfile(current_dir, 'settings.m');

        if ~exist(dest_path, 'file')
            copyfile(src_path, dest_path);
            fprintf('Copied %s to %s\n', src_path, dest_path);
        end

    catch err 
        disp(err);
    end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create path for other data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
today_data_dir = fullfile(root_dir, other_data_dir,...
    datestr(Today, 'yyyy'),...
    sprintf('%s-%s', datestr(Today, 'mm'), datestr(Today, 'mmm')),...
    datestr(Today, 'dd'));

if ~isdir(today_data_dir)
    mkdir(today_data_dir)
    fprintf('Successful created other-data directory %s\n', today_data_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% copy utility script files for trapping frequencies, etc. to folder
% In matlab 2016a and before folder is not a field of file object?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
files = dir(fullfile(script_files_dir, '*.m'));
text_files = dir(fullfile(script_files_dir, '*.txt'));
files = vertcat(text_files, files);

for ii = 1:length(files)

    try
    src_path = fullfile(script_files_dir, files(ii).name);
    dest_path = fullfile(today_data_dir, files(ii).name);
    if ~exist(dest_path, 'file')
        copyfile(src_path, dest_path);
        fprintf('Copied file %s\n', files(ii).name);
    end
    catch err
        disp(err);
    end
end

% copy diagnostic sequences into other data for a continuous record
diagnostic_dir = fullfile(diagnostics_files_dir, '*.seq');
diag_files = dir(diagnostic_dir);

for ii = 1:length(diag_files)

    try
        src_path = fullfile(diagnostics_files_dir, diag_files(ii).name);
        dest_path = fullfile(today_data_dir, diag_files(ii).name);
        if ~exist(dest_path, 'file')
            copyfile(src_path, dest_path);
            fprintf('Copied file %s\n', diag_files(ii).name);
        end
    catch err
        disp(err);
    end
end

% copy other files with paths relative to root_dir
for ii = 1 : length( other_file_paths )
    try
        src_path = fullfile(root_path, other_file_paths{ii});
        [~, fname, ext] = fileparts(src_path);
        dest_path = fullfile(today_data_dir, sprintf('%s%s', fname, ext));
        if ~exist(dest_path, 'file')
            copyfile(src_path, dest_path);
            fprintf('Coped file %s\n', fname);
        end
    catch err
        disp(err);
    end  
end

% copy other files with paths relative to new_day_dir
for ii = 1 : length( other_script_file_paths )
    try
        src_path = fullfile(new_day_dir, other_script_file_paths{ii});
        [~, fname, ext] = fileparts(src_path);
        dest_path = fullfile(today_data_dir, sprintf('%s%s', fname, ext));
        if ~exist(dest_path, 'file')
            copyfile(src_path, dest_path);
            fprintf('Coped file %s\n', fname);
        end
    catch err
        disp(err);
    end  
end
