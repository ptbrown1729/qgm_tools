function [Paths] = getDataSetPath(Date, DataSetNums, BoolProcessed)
%[Paths] = getDataSetPath(Date, DataSetNums)
%Returns full path for data set. To get a date, Date = datenum(2016,2,1); 
% e.g. for yyyy,mm,dd.
%
% Arguments:
%
% DataSetNums: is a 1D array. If it is a single number, then Paths is a char.
%If it is longer, Paths is a cell array of paths.
%
% BoolProcessed: If is 1, this function returns the path to the processed
%subfolder of the data set. Otherwise, it simply returns the data set
%folder.
%

RootImageDirectory = fullfile('//', '128.112.86.75', 'lithium', 'IMAGING DATA');
DateStr = fullfile(datestr(Date, 'yyyy'), datestr(Date, 'mm'), datestr(Date, 'dd'));

if ~exist('BoolProcessed', 'var')
    BoolProcessed = 0;
end

Paths = {};
for ii = 1:length(DataSetNums)
    DataSetNum = DataSetNums(ii);
    DatePath = fullfile(RootImageDirectory, DateStr);

    if isinf(DataSetNum) || floor(DataSetNum) ~= DataSetNum || ...
       DataSetNum < 0 || DataSetNum > 1000
        error('getDataSetPath DataSetNum was not an integer, or was greater than one thousand or less than zero');
    end
    
    folder_name = sprintf('%03d_*', DataSetNum);
    DirPath = fullfile(DatePath, folder_name);
    DataFile = dir(DirPath);
    
    if isempty(DataFile)
        disp('No Such Dataset')
        Path = '';
    elseif length(DataFile)>1
        disp('More than one such directory')
        Path = '';
    else
        if BoolProcessed
            Path = fullfile(DatePath,DataFile.name,'processed');
        else
            Path = fullfile(DatePath,DataFile.name);
        end
    end
    Paths = cat(2,Paths,Path);
end

if length(Paths)<=1
    Paths = char(Paths);
end

end

