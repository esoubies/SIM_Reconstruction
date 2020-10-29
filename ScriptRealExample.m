%--------------------------------------------------------------------------
% This script is an example of use of the SIM reconstruction method
% proposed in [1] on the real SIM data located within the folder Data/
% Patterns have been estimated with FairSIM (see ScriptPatternFromFairSIM.m) 
% and the PSF has been manually tuned (physical parameters not available).
% It is thus possible to further increase the reconstruction results with
% more accurate patterns and PSF.
%
% [1] E. Soubies and M. Unser, "Computational super-sectioning for single-slice structured-illumination microscopy" 
%     IEEE Transactions on Computational Imaging, vol. 5 no. 2, pp. 240-250, 2019. 
%
% Copyright (2018) E. Soubies emmanuel.soubies@irit.fr
%--------------------------------------------------------------------------
clear;close all;clc;

%% Parameters
% -- Paths
basedir=cd;
dataname=[basedir,'/Data/Acqdata.tif'];      % File name data image
psfname=[basedir,'/Data/PSF.tif'];           % File name psf
pattname=[basedir,'/Data/patterns.tif'];     % File name patterns
outFolder=[basedir,'/RealExample/'];         % Folder to save results
gtname=[];                    % File name ground truth (if no let empty)
sav=0;                        % Boolean if true save the result

% -- Data
valback=0;       % Background value
downFact=[2 2];  % Downsampling factors

% -- Objective Function
lamb=1e-4;       % Hyperparameter (can be an array to loop)
nbOutSl=1;       % The number of considered out-of-focus slices will be (2*nbOutSl-1)
symPsf=1;        % Boolean true if psf is symmetric (in this case 2*nbOutSl out-of-focus slides are considered on the same side of the psf in z)
Reg=2;           % Choice regul: 1 for TV, 2 for Hessian-Schatten

% -- SIM reconstruction
alg=1;                    % Algorithm (1: ADMM / 2: Primal-Dual)
maxIt= 100;               % Max iterations
ItUpOut=round(maxIt/10);  % Iterations between to call to OutputOpti
rhoDTNN= 1e-2;            % rho parameter (ADMM) associated to data term
rhoReg= 1e-2;             % rho parameter (ADMM) associated to the regularization (must be greater or equal than rhoDTNN if iterCG=0)
split=2;                  % Splitting strategy for data fidelity:
valId=2;                  % Scaling (>1) of the identity operator for the reformulation in [1] (only for splitting 3)

%% Run script
% - Reconstruction with no out-of-focus plane
run ./SimScript2D

pp=20; % Discard border effects for better visualization
figure;subplot(2,2,1);imagesc(xopt(pp:sz(1)*2-pp,pp:sz(2)*2-pp)); axis image; title('Reconstructed image');
subplot(2,2,2);imagesc(log(1+abs(fftshift(fftn(xopt))))); axis image; title('Reconstructed image FFT'); 
subplot(2,2,3);imagesc(sum(y,3)); axis image; title('Widefield image');
subplot(2,2,4);imagesc(log(1+abs(fftshift(fftn(sum(y,3)))))); axis image; title('Widefield image FFT'); 