function [fileName, directoryName] = findFileGeneral(seqId, pattern)
%
% Arguments:
%
% seqId: [year month day folder_index file_index]
%
% pattern: a string (usually including a wild card character "*")
% specifying the pattern of files to find in the given folder.
%
% Returns:
% 
% fileName: file path as a string
%
% directoryName: as a string

rootFolder = fullfile('//','128.112.86.75','lithium','IMAGING DATA');
%rootFolder = 'Z:/IMAGING DATA/'; %only works if you map file server to Z: drive. 
% changed on 09/05/16 by PB

%   Import the picture referenced in the line "index" of "fileContent"
%   and returns it in the structure "data"
% pattern can be like *abc*def* or *.txt or *_123_*.sif
%   -----------------------------------------------------------------------

if (~isnumeric(seqId) || length(seqId)~=5)
    fileName = [];
    directoryName = [];
    fprintf('error seqId not valid\n');
	seqId
    return;
end

%yearFolder = sprintf('%04d/',seqId(1));
yearFolder = sprintf('%04d', seqId(1));
monthStr = sprintf('%02d', seqId(2));
dayStr = sprintf('%02d', seqId(3));
datasetNoStr = sprintf('%03d', seqId(4));


%monthDirName = [rootFolder yearFolder monthStr '*'];
monthDirName = fullfile(rootFolder, yearFolder, sprintf('%s*', monthStr));
rootContent = dir(monthDirName); % adapted to new folder structure
rootContent = rootContent(cellfun(@(x)x,{rootContent(:).isdir}));  % kick out non-folder entries
if length(rootContent) ~= 1
    fprintf('error trying to access %s\n',monthDirName);
    fileName = [];
    directoryName = [];
    return;
end

monthStr = rootContent(1).name;

dayDirName = fullfile(rootFolder, yearFolder, monthStr, sprintf('%s*', dayStr));
rootContent = dir(dayDirName); % adapted to new folder structure
rootContent = rootContent(cellfun(@(x)x,{rootContent(:).isdir}));  % kick out non-folder entries

if length(rootContent) ~= 1
    fprintf('error trying to access %s\n',dayDirName);
    fileName = [];
    directoryName = [];
    return;
end

dirName = fullfile(rootFolder, yearFolder, monthStr, rootContent(1).name);

if (~rootContent(1).isdir)
    fprintf('error trying to access %s\n',dirName);
    fileName = [];
    directoryName = [];
    return;
end

dirContentSearchString = fullfile(dirName, sprintf('%s*',datasetNoStr));
dirContent = dir(dirContentSearchString);

dirContent = dirContent(cellfun(@(x)x,{dirContent(:).isdir}));  % kick out non-folder entries
if (isempty(dirContent))
    fprintf('folder not found %s\n',dirContentSearchString); 
    fileName = [];
    directoryName = [];
    return;
elseif (length(dirContent) > 1)
    fprintf('folder ambigious %s\n',dirContentSearchString); 
    fileName = [];
    directoryName = [];
    return;
end
        
% get name of the image file
directoryName = fullfile(dirName, dirContent(1).name);
fileNamePattern = fullfile(directoryName, pattern);

files = dir(fileNamePattern);
files = files(cellfun(@(x)~x,{files(:).isdir}));  % kick out non-file entries
S = [files(:).datenum].'; 
[~,idx] = sort(S);
files = files(idx); % Cell array of names in order by datenum. 

if (isempty(files))
    fileName = [];
    fprintf('no files found with %s\n',fileNamePattern);
    return;
end

if (seqId(5) <= numel(files))
    fileName = fullfile(directoryName, files(seqId(5)).name);
else
    fileName = [];
%     fprintf('file idx larger than number of files in folder\n');
end

end

