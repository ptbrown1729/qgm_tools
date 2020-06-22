# qgm_tools
This repository contains various tools I developed (or modified) while studying the Fermi-Hubbard model using a quantum gas microscope.

## lib
various useful functions for image processing, fitting, data analysis etc.

## qgm_analysis
code intended to analyze output of quantum gas microscope experiment

## realtime_analysis_pgm
code to analyze absorption images and display in real time

## lattice
optical lattice calculations. Right now there are what amounts to 2-3 major versions of this code stuffed in this folder. Some in python, some in matlab. The most recent is in python. Needs to be cleaned up, and previous versions should be removed

## gpib_devs
code for interacting with gpib devices using matlab. Some of these functions are based on code written by Debayan Mitra

## dqmc_tools
code used to automate running QuestQMC on the Princeton Feynman cluster, and parse the results. Several of the Matlab versions were originally written by Peter Schauss (http://ultracold.phys.virginia.edu/public_html/).

## optical-potentials
calculation of trap depths and etc.

## real-space-hamiltonians
Diagonalization of single-particle Hamiltonians with various potentials. Based on a tutorial given by Erich Mueller at the Princeton Center for Theoretical Science meeting on Topological and Strongly Correlated Phases in Cold Atoms in 2015
