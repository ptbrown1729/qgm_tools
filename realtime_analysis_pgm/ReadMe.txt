Updated 2017-02-25 by Peter Brown

A Brief Overview
-----------------------------------------
A simple framework for automated analysis of imaging data in realtime along with a number of useful data analysis functions.

I assume that imaging data is stored in a folder hierarchy based on date. The data folders should have the format \yyyy\mm\dd\nnn_anystring.

Each data folder contains image files (.fits or .aia are supported), a settings file "settings.m" which describes parameters for the analysis (e.g. imaging system magnification, initial guesses for fit parameters, etc) and a log file "log.txt"

The main function, called ProgramClass.m, runs continuously, searching for new files in todays date path (by default) and analyzes them according to the rules in the settings file.

My goal is to create a general frame work for analyzing images that rarely needs to be changed. Functions such as readfolder and its dependent function should not be edited, to the extent that is possible.

Using this code
-----------------------------------------
This code is maintained in a subversion (svn) repository on our file server (\lithium\SVN\RealTimeAnalysisSVN). To use it, checkout a copy to your local computer. For more information on how to do this, read about subversion online or in the documentation found in the svn folder.


Main Code -- Critical Files
-----------------------------------------
"ProgramClass.m"
The main function which keeps track of data sets, settings files, experiment analysis functions, writing logs, etc.

"ConstantsClass.m"
Physical constants data

"SettingsClass.m"
Settings files are transformed into SettingsClass instances. This is a more convenient way to transmit this data between functions in ProgramClass.m

"ImageDataClass"
A class for holding data about an image such as path, actual images read out, region of interest, etc. Convenient way to transmit this data between functions in ProgramClass.m

"DatasetClass.m"
Class for data set identification. Contains information about year, month, day, data set number, etc. Useful for checking consistency between different sorts of timestamps, comparing data sets (i.e. did one come later than the other?) etc. 

"\Helper Functions"
A collection of functions for doing things such as timestamp manipulations, file manipulations, log file writing and reading, coordinate transformations, etc. Used extensively by ProgramClass.m

"\experiment"
Contains functions which do the bulk of the data analysis. These are the files that should be modified for new types of experiments. i.e. if you change the number of pictures you are taking, or switch for Gaussian fitting to Thomas-Fermi you can add a new function here, then point to it using your settings file.

Additional Files
-----------------------------------------
"\azimuthalavg"
Data analysis functions supporting azimuthal averaging

"\display"
Functions related to displaying data. Probably the most useful is "watchLog.m" which continuously watches variables in a log file.

"\fileviewer"
A matlab GUI tool for visualizing od data. Not well maintained. Could use a lot of work.

"\Fitting"
Functions for fitting. "fit1D.m" and "fit2D.m" are powerful. They allow you to fix variables in your fit function and perform fitting of functions as long as they are of a standard for.

"\Fitting\Functions\2D"
2D functions, e.g. "gaussian2D.m"

"\Fitting\Functions\1D"
1D functions

"\imageproc"
Image processing functions. Calculate optical density, bin image data, normalize od, get rid of background stripes, get average of all pictures in folder, etc.

"\physdata"
Get physical data, such as atom number or intensity from od image. Phase space density, temperature from time-of-flight.

"\sandbox"
Files still under test. To be removed.

"\scripts\NewDay"
Script files for use with windows Task Scheduler to automatie basic directory creation tasks that need to be performed at the start of every day.

"\statistics"
Perform statistics on log file.

"aia.pdf"
Description of .aia file format which is used by Simplicio software (see ...)

"readimg.m"
main function for reading .aia or .fits files

"ImgAnalysis.py"
Data analysis functions similar to many of the matlab ones (e.g. azimuthal average) written in python (primarily using NumPy)


