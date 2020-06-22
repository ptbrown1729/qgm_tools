function readSubFolders(FolderName)
%Process all immediate subfolders of FolderName. This function


RootImageDirectory = '\\128.112.86.75\lithium\Imaging Data';
%FolderPath = fullfile(RootImageDirectory,FolderName,'*_*');
FolderPath = fullfile(RootImageDirectory,FolderName);

Continue = 1;
while Continue == 1
%     tic;
    SubFolders = {};
    Folders = dir(FolderPath);
    if ~isempty(Folders)
        %List subdirectories
        
        
        for ii = 1:length(Folders)
            SettingsPath = fullfile(RootImageDirectory,FolderName,Folders(ii).name,'settings.m');
            %exclude files we don't want.
            if Folders(ii).isdir && ~strcmp(Folders(ii).name(1),'.') && exist(SettingsPath,'file') && length(Folders(ii).name)>2
                %implement check for settings file. Otherwise don't add that file.
                CurrentFilePath = fullfile(FolderName,Folders(ii).name);
                SubFolders = {SubFolders{:},CurrentFilePath};
            end
        end
        readfolder(SubFolders,'Single');
        
    else
    end
%     IterationT = toc;
%     fprintf('readSubFolders Iteration took %0.3f Seconds \n',IterationT);
end


end

