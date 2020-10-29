%--------------------------------------------------------------------------
% This script reproduces the Figure 3 of the paper
%
% [1] E. Soubies and M. Unser, "Computational super-sectioning for single-slice structured-illumination microscopy" 
%     IEEE Transactions on Computational Imaging, vol. 5 no. 2, pp. 240-250, 2019. 
%
% Copyright (2018) E. Soubies emmanuel.soubies@irit.fr
%--------------------------------------------------------------------------
clear;close all;

%% Parameters
% -- Paths
basedir=cd;
dataname=[basedir,'/Data/Acqdata.tif'];      % File name data image
psfname=[basedir,'/Data/PSF.tif'];           % File name psf
pattname=[basedir,'/Data/patterns.tif'];     % File name patterns
gtname=[];                    % File name ground truth (if no let empty)
outFolder=basedir;             % Folder to save results
sav=0;                        % Boolean if true save the result

% -- Data
valback=0;       % Background value
downFact=[2 2];  % Downsampling factors

% -- Objective Function
lamb=1e-4;       % Hyperparameter (can be an array to loop)
nbOutSl=2;       % The number of considered out-of-focus slices will be (2*nbOutSl-1)
symPsf=1;        % Boolean true if psf is symmetric (in this case 2*nbOutSl out-of-focus slides are considered on the same side of the psf in z)
Reg=1;           % Choice regul: 1 for TV, 2 for Hessian-Schatten
