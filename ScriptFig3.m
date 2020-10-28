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

% -- Initialize Figures
figIters=figure; hold on; set(gca,'Fontsize',14); grid;xlabel('Iterations'); ylabel('f(x)/f(x^0)');
set(gca,'yscale','log')
figtime=figure; hold on; set(gca,'Fontsize',14); grid;xlabel('Elapsed time (min)'); ylabel('f(x)/f(x^0)');
set(gca,'yscale','log')

%% Algorithms
% -- Commun parameters
ItUpOut=5;   
maxIt=100;      % Max iterations

% *** ADMM ***
alg=1;
rhoDTNN= 7e-2;
rhoReg= 7e-2;

% - Full Splitting (iter GC = 1)  
split=1;       
iterCG=1;
run ./SimScript2D
F0=Opt.cost*zeros(H.sizein);
% Plots
figure(figIters)
semilogy([0,Opt.OutOp.iternum],[F0,Opt.OutOp.evolcost]/F0,'LineWidth',2);
drawnow;
figure(figtime)
semilogy(linspace(0,Opt.time/60,length(Opt.OutOp.evolcost)+1),[F0,Opt.OutOp.evolcost]/F0,'LineWidth',2);
drawnow;

% - Full Splitting (iter GC = 3)  
iterCG=3;
run ./SimScript2D
F0=Opt.cost*zeros(H.sizein);
% Plots
figure(figIters)
semilogy([0,Opt.OutOp.iternum],[F0,Opt.OutOp.evolcost]/F0,'LineWidth',2);
drawnow;
figure(figtime)
semilogy(linspace(0,Opt.time/60,length(Opt.OutOp.evolcost)+1),[F0,Opt.OutOp.evolcost]/F0,'LineWidth',2);
drawnow;

% - Splitting proposed in [1]
split=2;    
rhoDTNN= 7e-4;
rhoReg= 7e-4;
valId=2; 
run ./SimScript2D
F0=Opt.cost*zeros(H.sizein);
% Plots
figure(figIters)
semilogy([0,Opt.OutOp.iternum],[F0,Opt.OutOp.evolcost]/F0,'LineWidth',2);
drawnow;
figure(figtime)
semilogy(linspace(0,Opt.time/60,length(Opt.OutOp.evolcost)+1),[F0,Opt.OutOp.evolcost]/F0,'LineWidth',2);
drawnow;

% *** Primal-Dual (Splitting proposed in ref [12] of [1]) ***
alg=2;
rho=1.99;     
tau=2.5;
run ./SimScript2D
F0=Opt.cost*zeros(H.sizein);
figure(figIters)
semilogy([0,Opt.OutOp.iternum],[F0,Opt.OutOp.evolcost]/F0,'LineWidth',2);
drawnow;
figure(figtime)
semilogy(linspace(0,Opt.time/60,length(Opt.OutOp.evolcost)+1),[F0,Opt.OutOp.evolcost]/F0,'LineWidth',2);
drawnow;


figure(figIters)
legend('ADMM FS (1 iter CG)','ADMM FS (3 iters CG)', ...
       'ADMM Inner-loop free',...
       'Primal-Dual');
figure(figtime)
legend('ADMM FS (1 iter CG)','ADMM FS (3 iters CG)', ...
       'ADMM Inner-loop free',...
       'Primal-Dual');


