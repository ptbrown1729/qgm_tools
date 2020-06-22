SET FullPath=%~dp0 
rem get the path from the batch file location. Change to that drive.
%FullPath:~0,2% 
rem cd to the same folder as the batch file
cd "%FULLPATH%" 
matlab -r "try; p = ProgramClass(); catch; end; exit;"
PAUSE