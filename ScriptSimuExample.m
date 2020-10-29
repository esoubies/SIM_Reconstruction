%--------------------------------------------------------------------------
% This script is an example of use of the SIM reconstruction method
% proposed in [1] on a small simulated example without out-of-focus signal.
%
% [1] E. Soubies and M. Unser, "Computational super-sectioning for single-slice structured-illumination microscopy" 
%     IEEE Transactions on Computational Imaging, vol. 5 no. 2, pp. 240-250, 2019. 
%
% Copyright (2018) E. Soubies emmanuel.soubies@irit.fr
%--------------------------------------------------------------------------
clear;close all;clc;

%% Simulation
% -- Paths
gtpath='./Data/object.tif';          % file name ground truth 
expFolder='SimulatedExample/';   % experiment folder
sav=1;                             % to save results

% -- PSF
lamb=488;                % Illumination wavelength
res=32;                  % Resolution (nm)
Na=1.4;                  % Objective numerica aperture
nl=1.518;                % Refractive index of the objective medium (glass/oil)
ns=1.333;                % Refractive index of the sample medium (water)

% -- Patterns
orr=[0 pi/3 2*pi/3];   % Patterns orientations (vector)                
ph=linspace(0,2*pi,4); % Patterns lateral phases (vector)
ph=ph(1:end-1);  
a=0.9;                 % Amplitude coefficient 
bet=asin(Na/nl);       % Angle between side beams and the optic axis (e.g. bet asin(Na/nl))

% -- Acquisition
downFact=[2 2];  % Downsmpling factor (e.g. [2 2 2]) 
photBud=500;    % Photon Budget

% -- Run Simulator
disp('######## Generate SIM data ########');
run './SimuSIM2D.m'
disp('###################################');

%% Reconstruction 
% -- Paths
basedir='./SimulatedExample/';
dataname=[basedir,'/AcqData.tif'];      % File name data image
psfname=[basedir,'/PSF'];               % File name psf
pattname=[basedir,'/patterns.tif'];     % File name patterns
gtname=['./Data/object.tif'];           % File name ground truth (if no let empty)
outFolder=[basedir];                    % Folder to save results
sav=1;                                  % Boolean if true save the result

% -- Data
valback=0;       % Background value
downFact=[2 2];  % Downsampling factors

% -- Objective Function
lamb=2e-4;       % Hyperparameter (can be an array to loop)
nbOutSl=0;       % The number of considered out-of-focus slices will be (2*nbOutSl-1)
symPsf=1;        % Boolean true if psf is symmetric (in this case 2*nbOutSl out-of-focus slides are considered on the same side of the psf in z)
Reg=2;           % Choice regul: 1 for TV, 2 for Hessian-Schatten

% -- SIM reconstruction
alg=1;                    % Algorithm (1: ADMM / 2: Primal-Dual)
maxIt= 100;               % Max iterations
ItUpOut=round(maxIt/10);  % Iterations between to call to OutputOpti
rhoDTNN= 1e-3;            % rho parameter (ADMM) associated to data term
rhoReg= 1e-3;             % rho parameter (ADMM) associated to the regularization (must be greater or equal than rhoDTNN if iterCG=0)
split=2;                  % Splitting strategy for data fidelity:
valId=2;                  % Scaling (>1) of the identity operator for the reformulation in [1] (only for splitting 3)

%% Run script
disp('########## Reconstruction #########');
run ./SimScript2D
disp('###################################');
figure;subplot(1,2,1);imagesc(xopt); axis image; title('Reconstructed image');
subplot(1,2,2);imagesc(log(1+abs(fftshift(fftn(xopt))))); axis image; title('Reconstructed image FFT'); 