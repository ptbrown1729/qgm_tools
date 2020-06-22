rem Windows batch file for creating new folders and copying useful scripts
rem into the "OTHER DATA" directory. This file calls the matlab script "NewDay.m"
rem which does the bulk of the work. This batch file can be easily 
rem run every morning using windows task scheduler.

SET FullPath=%~dp0 
rem get the path from the batch file location. Change to that drive.
%FullPath:~0,2% 
rem cd to the same folder as the batch file
cd "%FULLPATH%" 

rem for debugging, uncomment and run the following lines:
rem matlab -nosplash -nodesktop -r "run('NewDay.m');"
rem pause

rem normal operation, run this line. This opens matlab and runs the "NewDay.m" script.
matlab -nosplash -nodesktop -r "try; run('NewDay.m'); catch; end; exit;"
