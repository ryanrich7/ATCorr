ATCorr

This is a simple class which extracts the kinematic quantities necessary for the AT correction for CREX. 

Ryan Richards, Sept 2021

ca48 tables -- Ca48 Horowitz tables for the cross section and asymmetry for a given beam energy and angle

CREXdata.h -- header file needed by GetKine which reads the Horowitz and computes the asymmetry for each event.

GetKine.C -- Draws the kinematic distributions for a given run and detector (upstream, downstream, ATs) and saves them to a pdf. The central values and the stat errors for each histogram is printing out for a given run. You can store that information in a file. 


LeftMainDet.txt and RightMainDet.txt -- Example files for the upstream main detectors

ATCorr.h -- header 

ATCorr.C -- implementation

An instance of the class takes two files and reads in the data.  

Main Methods
Init -- computes kinematic quantites for detector combinations 

PlotMainDetPlots -- Fits the data for the main detectors (can be upstream or downstream) and their combinations

PlotATDetPlots -- Fits the data for the AT detectors and their combinations
