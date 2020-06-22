
function readfolder(FolderNames,Mode)
%readfolder(FolderNames,Mode)
%Takes a cell array FolderNames, each relative to the root image directory.
%Mode is 'Continous' or 'Single'. If no second argument is given, will run
%Continously.
%The root data anlysis function. This script should run continuously while
%data is being taken. It processes data in real time using other functions.
%All information from these functions is saved in a text file in the process
%directory.

%Check matlab version
versionDate = version('-date');
if datenum(versionDate)<datenum(2015,07,1)
    warning('Parts of this program were written using matlab 2015b. Some functions including showPlots are known not to run on matlab 2014 and earlier.');
end
%Directory storing images.
RootImageDirectory = '\\128.112.86.75\lithium\Imaging Data';
%Directory storing analysis scripts
RootAnalysisPath = '\\128.112.86.75\lithium\Imaging Data\Image Analysis';
%Add relevant subdirectories to matlab path. Tried genpath, but was very slow.
addAnalysis2Path(RootAnalysisPath)

%class storing fundamental constants
Constants = ConstantsClass;

%List for each folder being read
SavePaths = {};
SettingsList = [];
ProcessPaths = {};
LogPaths = {};
ExperimentFnHandlers = {};

if ischar(FolderNames)
    FolderNames = {FolderNames};
end

%For each folder, generate paths and settings files.
for ff = 1:length(FolderNames)
    %location image files are being saved.
    CurrentFolderName = FolderNames{ff};
    CurrentSavePath = fullfile(RootImageDirectory,CurrentFolderName);
    SavePaths = cat(2,SavePaths,CurrentSavePath);
    
    %settings files
    CurrentSettingsPath=fullfile(CurrentSavePath,'settings.m');
    CurrentSettings = classFromSettingsFile(char(CurrentSettingsPath));
    SettingsList = cat(2,SettingsList,CurrentSettings);
    
    %Create folder and log file for processed data.
    CurrentProcessPath = fullfile(RootImageDirectory,CurrentFolderName,'\processed');
    ProcessPaths=cat(2,ProcessPaths,CurrentProcessPath);
    if ~isdir(CurrentProcessPath)
        mkdir(CurrentProcessPath);
    end
    CurrentLogPath=fullfile(CurrentProcessPath,'log.txt');
    LogPaths = cat(2,LogPaths,CurrentLogPath);
    
    %Get function handle of analysis function. experiment variable in settings
    %file directs to the appropriate function.
    %cd(fullfile(RootAnalysisPath,'Experiment')) %This is added to matlab
    %path by addAnalysis2Path above...
    ExperimentFnHandlers = {ExperimentFnHandlers{:},str2func(CurrentSettings.Experiment)};
end

%loop checking for files in folder
if ~exist('Mode','var')
    Mode = 'Continuous';
end
Continue = 1;

if strcmp(Mode,'Continuous')
    fprintf('starting loop...\n');
end

while Continue == 1
    tic;
    for nfolder = 1:length(FolderNames)
        %fprintf('processing folder #%d\n',nfolder);
        %move to folder and begin
        if exist(SavePaths{nfolder},'dir')
            files = [dir(fullfile(SavePaths{nfolder},'*.aia'));dir(fullfile(SavePaths{nfolder},'*.fits'))];
            %If there are image files, process them.
            if ~isempty(files)
                for ii=1:length(files)
                    FileName=files(ii).name;
                    disp(FileName)
                    Moved = 0;
                    while Moved == 0
                        try
                            %Move file to processed folder.
                            Location=fullfile(SavePaths{nfolder},FileName);
                            Destination = fullfile(ProcessPaths{nfolder},FileName);
                            movefile(Location,Destination);
                            
                            %Rename fits files using date/time metadata
                            %stored in the file.
                            if strcmp(FileName(end-4:end),'.fits')
                                Destination = timeStampFitsFile(Destination);
                                [~,FileName,Ext] = fileparts(Destination);
                                FileName = strcat(FileName,Ext);
                            else
                            end
                            
                            Moved = 1;
                            
                        catch Exception
                            MsgString = getReport(Exception);
                            disp('Caught exception in readfolder trying to process file.');
                            disp(MsgString);
                            pause(2);
                            continue
                        end
                    end
                    
                    try
                        %MAIN ANALYSIS PART OF FUNCTION
                        %run functions to extract data from file
                        [Images,~,Stamp,Vals,Keys] = readimg(Destination);
                        %analysis tasks
                        CurrentSettings = SettingsList(nfolder);
                        CurrentExperimentFnHandler = ExperimentFnHandlers{nfolder};
                        [MoreVals,MoreKeys,ImgCrops,Fits]=CurrentExperimentFnHandler(Images,Vals,Keys,CurrentSettings,Constants);
                        Vals=cat(2,MoreVals,Vals);
                        Keys=cat(2,MoreKeys,Keys);
                        %display data
                        displayData(Vals,Keys);
                        %Show images
                        Continue = showPlots(ImgCrops,Fits,CurrentSettings,Vals,Keys,Stamp);
                        %file processed
                        disp(strcat(FolderNames{nfolder},'/',FileName,' processed'));
                        %write information to log file
                        CurrentLogPath = LogPaths{nfolder};
                        writeLog(CurrentLogPath,Vals,Keys,Stamp);
                    catch Exception
                        MsgString = getReport(Exception);
                        disp(MsgString);
                        %write information to log file
                        CurrentLogPath = LogPaths{nfolder};
                        writeLog(CurrentLogPath,[0],{'error'},Stamp);
                    end
                    
                end
            end
            pause(0.5);
        else
            disp(strcat(SavePaths{nfolder},' does not exist'));
        end
    end
    if strcmp(Mode,'Single')
        Continue = 0;
    end
    loopT = toc;
    %fprintf('readfolder loop T = %0.2f s \n',loopT);
end
