# qgm_tools
This repository contains various tools I developed while studying the Fermi-Hubbard model using a quantum gas microscope. Various functions were written by other members of the Bakr lab at Princeton University.

## lib
various useful functions for image processing, fitting, data analysis etc.

## qgm_analysis
code intended to analyze output of quantum gas microscope experiment. The most important file is [DataFolder.m](qgm_analysis/DataFolder.m),
which was used to analyze a single dataset, consisting either of many images taking at the same parameters,
or a series of images taken while one parameter varied. Many of the other files, such as [arpes_data.m](qgm_analysis/arpes_data.m), 
[DiffuseSet.m](qgm_analysis/DiffuseSet.m), and [CorrDataSetAtt.m](qgm_analysis/CorrDataSetAtt.m) are classes which organize many different kinds of data together 
(e.g. different spin states, etc.), and perform analysis that requires integrating information from these various datasets.
The individual data sets are typically instances of the class defined in DataFolder.m. These files are dependent on the 
functions available in [lib](lib), and an effort has been made to factor out most usefu functions and place them there, 
while using these files to apply those functions to certain types of data. Many of the other files in this folder are not
 useful and should be removed.

## realtime_analysis_pgm
Code to analyze absorption images and display in real time. The general idea is: data from a camera is saved in a folder,
and analyzed as it arrives. The specific analysis to be performed is described by a combination of a settings files (settings.m)
located in the same directory the image files are saved, and an analysis script contained in [scripts](realtime_analysis_pgm/scripts)
which is specified by the settings file. The results of the analysis are saved in a log file.

TODO: settings.m is now a reserved name in the most recent versions of Matlab, so this name needs to be
changed.  

The main program is run through [ProgramClass.m](realtime_analysis_pgm/ProgramClass.m). Settings files are parsed and
stored in instances of the class defined in [SettingsClass](realtime_analysis_pgm/SettingsClass.m).
Instances of [ConstantsClass.m](realtime_analysis_pgm/ConstantsClass.m) contain physical data, including
information about the specific atomic species. [DatasetClass.m](realtime_analysis_pgm/DatasetClass.m).
[ImageDataClass.m](realtime_analysis_pgm/ImageDataClass.m). These files may relay on functions in [lib](lib),
and these will need to be added to the Matlab path before running this program.

[watchLog.m](realtime_analysis_pgm/watchLog.m) is a utility to continuously check the log file for new information and plot
certain desired portions of it.  

[scripts/NewDay](realtime_analysis_pgm/scripts/NewDay) contains an example of a collection of settings files and folders
that were automatically created each night in preparation for the next days calibration experiments.

The files in [analysis_old_style](realtime_analysis_pgm/analysis_old_style) are retained mostly to give other examples of
what analysis script files should look like.

## lattice
optical lattice calculations. The most important file is [lattice_single_particle.py](lattice/lattice_single_particle.py). A number of examples
of its usage are found in [lattice/scripts](lattice/scripts). Other useful code for fitting the lattice depth based on estimation of
the band structure (from e.g. intensity modulation spectroscopy) are found in [lattice/lattice_depth_fitting](lattice/lattice_depth_fitting).
[lattice.nb](lattice/lattice.nb) contains some useful plots of the lattice profile in realspace. [lattice_notes.tex](lattice/lattice_notes.tex) 
contains information about the "standard" 2D optical lattice, as well as the D4 lattice that we are more interested in
here. [library.bib](lattice/library.bib) is a bibtex file which contains bibliographic information for useful references for
some of the calculations included here. [fermi_gas.py](lattice/fermi_gas.py) performs useful calculations of correlators
and etc. for a non-interacting Fermi gas in a lattice. Complimentary matlab functions can be found in [lib/non-int-fg](lib/non-int-fg).
[lattice_unittest.py](lattice/lattice_unittest.py) has tests to verify certain functions in lattice_single_particle.py
  
## gpib_devs
code for interacting with GPIB devices using matlab.

## dqmc_tools
code used to automate running QuestQMC on the Princeton Feynman cluster, and parse the results.

## optical-potentials
calculation of trap depths, atomic polarizabilities, etc.

## real-space-hamiltonians
Diagonalization of single-particle Hamiltonians with various potentials. Based on a tutorial given by Erich Mueller at
 the Princeton Center for Theoretical Science meeting on Topological and Strongly Correlated Phases in Cold Atoms in 2015.
 These files have not been touched since 2016, but some of the ideas are used in [lattice/lattice_single_particle.py](lattice/lattice_single_particle.py). 

