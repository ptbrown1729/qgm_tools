function [values, keys, stamps] = parseLog(log_path)
%[Vals,Keys,Stamps]=parseLog(LogPath)
%Vals is an NxM array, where N is the number of log entries, M is the
%number of Keywords.
%Keys is a 1xM cell array of key words
%In log file, every image must have the same keywords.

% TODO: I wrote this function before I was familiar with higher level
% matlab functions like dlmwrite. Might be worth reformatting how the
% log files work, i.e. using headings instead of including the key name
% next to the value on each line of the file. I did this initially because
% I couldn't guarantee that the lines would be the same between different
% .aia files. But if that is the case, the log file isn't going to be
% useful anyway, so it might be better to return an error.

%for parsing log files created with readfolder.m
keys=[];
values=[];
stamps = {};

if ~exist(log_path, 'file')
    warning('%s is not a valid file path.', log_path);
    return
end

try
    fid = fopen(log_path, 'r');
    
    %parse log file and return array giving information of all files stored
    if fid ~= -1
        first_line = fgetl(fid);
        first_line_split = strsplit(first_line, ',');
        n_entries = length(first_line_split);
        indices_keys = 2:2:n_entries;
        indices_vals = 1:2:n_entries - 1;
        values = cat(1, values, first_line_split(indices_vals));
        keys = cat(1, keys, first_line_split(indices_keys));

        % read each line of the file, and split each line at commas
        while ~feof(fid)
            %get current line
            line = fgetl(fid);
            line_split = strsplit(line,',');
            if length(line_split) == n_entries
                % if our line matches others
                values = cat(1, values, line_split(indices_vals));
                keys = cat(1, keys, line_split(indices_keys));
            else
                % if there was an error on this line, return nans
                nan_strings = cell(1, length(indices_vals));
                nan_strings(:) = {'nan'};
                values = cat(1, values, nan_strings);
                empty_strings = cell(1, length(indices_keys));
                empty_strings(:) = {''};
                keys = cat(1, keys, empty_strings);
            end
        end
    else
        warning('parseLog.m Failure when tried to open file %s.', log_path)
        fclose('all');
    end
    
    fclose(fid);
    
catch Exception
    MsgString = getReport(Exception);
    disp(MsgString);
    
	try
		fclose(fid);
    catch Exception
        MsgString = getReport(Exception);
        disp(MsgString);
	end

    stamps = '';
    values = [];
    keys = {};
end

% if the file was read successfully, get file stamps, only return
% the keys list once, and convert everything to the correct types.
if ~isempty(values)
    stamps = values(:, 1);
    values=str2double(values(:, 2:end));
    keys=keys(1, 2:end);
else
    warning('Did not obtain any data from file %s.', log_path);
end

end
