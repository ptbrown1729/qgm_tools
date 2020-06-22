function [duplicates] = find_duplicate_fnames(dir_paths)
% Find duplicate file names contained within a directory and any
% subdirectories. This is useful for searching for duplicate file names in
% a project.
% 
% duplicates is a cell array. Each cell is another cell array, giving
% information about duplicates. This second level of cell arrays has the
% form {file_name, {directory_1, directory_2, ..., directory_n} }. So 
% duplicates = { {file_name_1, {directory_11, directory_12, ...} }, 
% {file_name_2, {directory_21, directory_22,...} }, ...};

    if ischar(dir_paths)
        dir_paths = {dir_paths};
    end

    if ~exist('dir_paths', 'var') || isempty(dir_paths)
        dir_paths = {pwd};
    end
    
    sub_dirs = {};
    for ii = 1:length(dir_paths)
        sub_dirs_temp = strsplit( genpath(dir_paths{ii}), ';');
        sub_dirs = horzcat(sub_dirs, sub_dirs_temp);
    end

%     sub_dirs = strsplit( genpath(dir_paths{1}), ';' );
        
    % ignore hidden files
    % '.*' = any string
    % sprintf('%s.', filesep) will match \. on a windows machine, which is
    % how a hidden file will start. With this expression, we will exclude
    % folders that are contained in hidden folders also.
    hidden_rexp = sprintf('.*%s..*', filesep); 
    indices_hidden_files = [];
    if 1
        for ii = 1:length(sub_dirs)
            % hijack file parts to get non-empty btm_dir if and only if
            % directory is hidden.
            matches = regexp(sub_dirs{ii}, hidden_rexp, 'match');
            if ~isempty(matches)
                indices_hidden_files = cat(2, indices_hidden_files, ii);
            end
        end
        sub_dirs(indices_hidden_files) = [];
    end
    
    % get list of all files
    files = [];
    for ii = 1:length(sub_dirs)
        files_temp = dir(sub_dirs{ii});
        files_temp = files_temp(~[files_temp.isdir]);
        files = vertcat(files, files_temp);
    end
    
    fnames = {files.name};
    duplicates = {};
    duplicates_indices = [];
    for ii = 1:length(files)
        if ~any(ii == duplicates_indices)
        
            fn = @(name) strcmp(name, fnames{ii});
            duplicate_log_array = cellfun(fn, fnames);
            if sum(duplicate_log_array) < 2
            else
                indices = find(duplicate_log_array);
                duplicates_indices = cat(2, duplicates_indices, indices);
                dup_files = files(indices);
                dup_cell = {dup_files(1).name, {dup_files.folder}};
                duplicates = cat(2, duplicates, {dup_cell});
            end
        end
    end
    
end