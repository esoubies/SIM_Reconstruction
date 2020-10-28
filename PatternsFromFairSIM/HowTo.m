%--------------------------------------------------------------------------
% This script shows how to generate patterns from FairSim estimated
% parameters.
%
% Copyright (2020) E. Soubies emmanuel.soubies@irit.fr
%--------------------------------------------------------------------------
clear; clc;
addpath ../Utils/

% The file FairSimParameters.xml contains the patterns parameters for the 
% SIM images in the Data/ folder that have been estimated with the FairSIM 
% ImageJ plugin. 

% Below is the line to execute in order to obtain the SIM patterns saved as
% patterns.tif in the Data/ folder.
[patt,~,~,~]=PatternsFromFairSimXml('FairSimParameters.xml');