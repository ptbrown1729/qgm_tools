function [Dataset, PathAndStampConsistent, AllFilesCorrectlyNamed] = path2dataset(Path, Warnings)
%[Dataset,PathAndStampConsistent,AllFilesCorrectlyNamed] = path2dataset(Path)
%Dataset is a DatasetClass object.
%PathAndStampConsistent is a boolean describing whether or not date from
%path and time stamp agree.
%AllFilesCorrectlyNamed is a boolean describing whether or not all files
%have the correct time stamps. If not, Dataset.FileIndex will be set to 0.
%Dataset.PictureIndex is always set to 1.

if ~exist('Warnings', 'var')
    Warnings = 0;
end

    [PathStr,Name,Ext] = fileparts(Path); %fileparts treats the last colon as a delimeter if no slashes are present.
    FName = strcat(Name,Ext);

%get information from path
%assuming we have */yyyy/mm/dd/nnn_*/*.fits, match this with a regulatr
%expression.
try
    %     Exp = ['.*\(?<year>\d+)\(?<month>\d+)\(?<day>\d+)\(?<DataSetIndex>\d{3})_.*\.*.fits'];
    %     [Tokens,Matches] = regexp(Path,Exp,'tokens','match');
    %     PathYear = str2double(Tokens{1}{1});
    %     PathMonth = str2double(Tokens{1}{2});
    %     PathDay = str2double(Tokens{1}{3});
    %     PathDataSetIndex = str2double(Tokens{1}{4});
    %     PathFailed = 0;
    Folders = strsplit(PathStr,'\');
    PathYear = str2double(Folders(end-3));
    PathMonth = str2double(Folders(end-2));
    PathDay = str2double(Folders(end-1));
    DataFolder = Folders(end);
    DataFolderChar = char(DataFolder);
    PathDataSetIndex = str2double(DataFolderChar(1:3));
    PathFailed = 0;
    
    if isnan(PathYear)||isnan(PathMonth)||isnan(PathDay)||isnan(PathDataSetIndex)
        PathFailed = 1;
    end
catch
    PathFailed = 1;
    PathYear = 0; PathMonth = 0; PathDay = 0; PathDataSetIndex = 0;
end

%this same info should come from timestamp
try
    %use regexp...pretty gory, but just means yyyy-mm-dd-hh;mm;ss with
    %optional file extension
    StampExp = ['(?<year>\d+)-(?<month>\d+)-(?<day>\d+)-(?<hour>\d+);(?<minute>\d+);(?<second>\d+)(.fits.|aia)?'];
    [Tokens, Matches] = regexp(FName,StampExp,'tokens','match');
    StampYear = str2double(Tokens{1}{1});
    StampMonth = str2double(Tokens{1}{2});
    StampDay = str2double(Tokens{1}{3});
    StampHour = str2double(Tokens{1}{4});
    StampMinute = str2double(Tokens{1}{5});
    StampSecond = str2double(Tokens{1}{6});
    StampFailed = 0;
catch
    StampFailed = 1;
    StampYear = 0; StampMonth = 0; StampDay = 0; StampHour = 0; StampMinute = 0; StampSecond = 0;
end

%reconcile path and stamp timestep information
if PathFailed&&StampFailed
    PathAndStampConsistent = 0;
    Year = 0; Month = 0; Day = 0; Hour = 0; Minute = 0; Second = 0;
    DataSetIndex = 0;
elseif PathFailed
    PathAndStampConsistent = 0;
    Year = StampYear; Month = StampMonth; Day = StampDay;
    Minute = StampMinute; Hour = StampHour; Second = StampSecond;
    DataSetIndex = 0;
elseif StampFailed
    PathAndStampConsistent = 0;
    Year = PathYear; Month = PathMonth; Day = PathDay; DataSetIndex = PathDataSetIndex;
    Minute = 0; Hour = 0; Second = 0;
else
    PathAndStampConsistent = (PathYear==StampYear&&PathMonth==StampMonth&&PathDay==StampDay);
    if ~PathAndStampConsistent && ~PathFailed && ~StampFailed
        warning('Path and time stamp are inconsistent');
    end
    
    DataSetIndex = PathDataSetIndex;
    Year = StampYear; Month = StampMonth; Day = StampDay;
    Minute = StampMinute; Hour = StampHour; Second = StampSecond;
end


%get shot number
try
    Files = [dir(fullfile(PathStr,'*.fits')),dir(fullfile(PathStr,'*.aia'))];
    NameArray = {Files.name};
    %check if all names are in the correct format.
    NumberMatching = sum(~cellfun(@isempty,regexp(NameArray,StampExp)));
    if NumberMatching == length(Files) && ~isempty(Files)
        %if all files have the correct format, we know we are getting the
        %correct index here.
        [~,FileIndex] = find(strcmp(NameArray,FName));
        AllFilesCorrectlyNamed = 1;
    elseif NumberMatching == 0 || isempty(Files)
        FileIndex = 0;
        if Warnings
            warning('File not found in directory. File index set to zero.');
        end
    else
        %if they don't, still return index, but raise warning.
        [~,FileIndex] = find(strcmp(NameArray,FName));
        if Warnings
            warning('Not all files in this folder have the correct time stamps. File order may change.');
        end
        AllFilesCorrectlyNamed = 0;
    end
catch
    FileIndex = 0;
end

%asign values to DatasetClass fields.
Dataset = DatasetClass('');
Dataset.Year = Year;
Dataset.Month = Month;
Dataset.Day = Day;
Dataset.Hour = Hour;
Dataset.Minute = Minute;
Dataset.Second = Second;
Dataset.DataSetIndex = DataSetIndex;
Dataset.FileIndex = FileIndex;
Dataset.PictureIndex = 1;
end