%----------------------------------------------------
% This script performs 2D SIM Reconstruction using the ADMM of [1] with
%    - Data-term: Least-Squares 
%    - regul: TV or Hessian-Schatten norm
%
%  arg min_x \sum_p ||SHW_px - y_p||_2^2 + lamb*R(Lx) + i_{>=0}(x)
%
% where W_p is a diagonal operator moddeling the p-th structured-illumination
%       H is the convolution operator (with PSF)
%       S is a downsampling operator
%
% References:
% [1] E. Soubies and M. Unser, "Computational super-sectioning for single-slice structured-illumination microscopy" 
%     IEEE Transactions on Computational Imaging, vol. 5 no. 2, pp. 240-250, 2019. 
%
% Before to run it set the following papareters:
% *** Paths *** 
%  outFolder=      -> Folder to save results
%  dataname=       -> File name data image
%  psfname=        -> File name psf
%  gtname=         -> File name ground truth (if no let empty)
%  pattname=       -> File name patterns 
%  sav=            -> Boolean if true save the result
%
% *** Data ***
%  valback=        -> Background value
%  downFact=       -> Downsmpling operator (usually [2 2])
%
%  Both data (y) and patterns (w) should be tif 2D stacks
%    - [w1,w2,w3,...,wP] 
%    - [y1,y2,y3,...,yP]
%  where yp is the 2D data associated to the wp pattern.
%
% *** Objective Function ***
%  lamb=           -> Hyperparameter (can be an array to loop)
%  Reg=            -> Choice regul: 1 for TV, 2 for Hessian-Schatten 
%  nbOutSl=        -> Number of considered out-of-focus slices
%  symPsf=         -> Boolean true if psf is symmetric 
%                     (in this case 2*nbOutSl out-of-focus slides are considered on the same side of the psf in z)
%
% *** SIM reconstruction ***
%  maxIt=          -> Max iterations
%  alg=            -> Algorithm (1: ADMM / 2: Primal-Dual)
%  ItUpOut=        -> Iterations between to call to OutputOpti
%
% -- ADMM
%  rhoDTNN=        -> rho parameter (ADMM) associated to data term
%  rhoReg=         -> rho parameter (ADMM) associated to the regularization
%                     (must be greater or equal than rhoDTNN if iterCG=0)
%  split=          -> Splitting strategy for data fidelity:  
%                       0: no-splitting
%                       1: u_p=SHW_px (full splitting)
%                       2: u_p=W_px   (reformulation proposed in [1] is used)
%  iterCG=         -> Max number of conjugate gradient iterations for splittings (0 and 1)
%  valId=          -> Scaling (>1) of the identity operator for the reformulation in [1] (only for splitting 2)
%
% -- Primal-Dual (Splitting strategy proposed in ref [12] of [1])
%  tau=            -> tau parameter (Primal-Dual)
%  rho=            -> rho parameter (Primal-Dual) in ]0,2[
%
% Copyright (2018) E. Soubies emmanuel.soubies@irit.fr
%----------------------------------------------------
global isGPU
addpath Utils

%% Reading data
% -- PSF
if strcmp(psfname(end-2:end),'tif')
    psf=double(loadtiff(psfname));psf=psf/sum(psf(:));
else
    load(psfname);
end
% -- Patterns
if ~isempty(pattname)
    patt=double(loadtiff(pattname));
    if min(patt(:)) < 0, patt=patt-min(patt(:)); end                       % Ensure that patterns are positive]
    nbPatt=size(patt,3);
end
% -- Data
y=double(loadtiff(dataname))-valback;maxy=max(y(:));y=y/maxy;              % Load and normalize data
% -- Ground Truth
if ~isempty(gtname)                                                        % Load ground truth if available
    gt=double(loadtiff(gtname))-valback;            
else
    gt=[];
end
% -- Some parameters
sz=size(y);
szUp=size(psf);szUp=szUp(1:2); 
if alg==2,split=1;  end

%% Conversion CPU/GPU is necessary
psf=gpuCpuConverter(psf);
y=gpuCpuConverter(y);

%% Common operators and costs
% -- Data term
szPsf=size(psf);
if size(psf,3)==1
    fp=1;
else
    fp=floor(szPsf(3)/2)+1;
end   
if nbOutSl==0
    psf=psf(:,:,fp);
else
    if symPsf
        psf=psf(:,:,[fp+nbOutSl:fp+nbOutSl*2-1,fp:fp+nbOutSl-1]);      
    else
        psf=psf(:,:,fp+nbOutSl:-1:fp-nbOutSl+1);
    end
end
psf=psf/sum(psf(:));
otf=fftn(fftshift(fftshift(psf,1),2));
H=LinOpConv(otf,1); 
if nbOutSl==0  % Downsmpling operator
    S=LinOpDownsample(szUp,downFact);
else
    S=LinOpDownsample([szUp,size(otf,3)],[downFact,size(otf,3)]);
end
Fn={};Hn={};F0=[];
if ~isempty(pattname)
    nn=sum(patt.^2,3);                      % Trick to perform the last ADMM step in Fourier
    patt=patt/(sqrt(max(nn(:))));
    for i=1:size(patt,3)
        L2=CostL2([],y(:,:,i));             % L2 cost function
        switch split
            case 0
                if i==1
                    F0=L2*S*H*LinOpDiag(H.sizein,patt(:,:,i));
                else
                    F0=F0+L2*S*H*LinOpDiag(H.sizein,patt(:,:,i));
                end
            case 1
                Fn{i}=L2;
                Fn{i}.doPrecomputation=true;
                Hn{i}=S*H*LinOpDiag(H.sizein,patt(:,:,i));
            case 2
                Fn{i}=(L2*(S*H));
                Fn{i}.doPrecomputation=true;
                Hn{i}=LinOpDiag(H.sizein,patt(:,:,i));
        end
    end
else
    L2=CostL2([],y);             % L2 cost function
    LS=(L2*(S*H));
    LS.doPrecomputation=true;
end
% -- Regularization
if Reg==1
    Opreg=LinOpGrad(H.sizein,[1 2]);             % TV regularizer: Operator Gradient
    Freg=CostMixNorm21(Opreg.sizeout,4);         % TV regularizer: Mixed Norm 2-1
elseif Reg==2
    Opreg=LinOpHess(H.sizein,'circular',[1 2]);  % Hessian-Shatten: Hessian Operator
    Freg=CostMixNormSchatt1(Opreg.sizeout,1);    % Hessian-Shatten: Mixed Norm 1-Schatten (p=1)
end
% -- Non-Negativity constraint
pos=CostNonNeg(H.sizein);
if alg==1 && (~isempty(pattname) && split==2)
    assert(valId >= 1,'valId must be greater than 1');
    OpD=LinOpDiag(H.sizein,sqrt(valId-sum(patt.^2,3)));
else
    OpD=LinOpIdentity(H.sizein);
end
%% SIM Reconstruction 
for ii=1:length(lamb)
    if alg==1     % ADMM
        FF=[Fn,{pos},{lamb(ii)/(size(otf,3))*Freg}];
        HH=[Hn,{OpD},{Opreg}];
        rho=[ones(size(Fn))*rhoDTNN,rhoDTNN,rhoReg];
        Opt=OptiADMM(F0,FF,HH,rho);
        Opt.OutOp=MyOutputOpti(1,gt,round(maxIt/10),[1:nbPatt,nbPatt+2]);
    elseif alg==2 % Primal-Dual
        FF=[Fn,{lamb(ii)/(size(otf,3))*Freg}];
        HH=[Hn,{Opreg}];
        Opt=OptiPrimalDualCondat(F0,pos,FF,HH);
        Opt.tau=tau;
        Opt.rho=rho;
        T=HH{1}.makeHtH();for ll=2:length(HH), T=T+HH{ll}.makeHtH(); end;
        Opt.sig=1/(tau*T.norm);
        Opt.OutOp=MyOutputOpti(1,gt,round(maxIt/10),2:(2+nbPatt));
    end
    Opt.OutOp.saveXopt=0;
    Opt.CvOp=TestCvgStepRelative(1e-5);
    if alg==1 && (split~=2 && iterCG~=0)
        Opt.CG.maxiter=iterCG;
    end
    Opt.ItUpOut=ItUpOut;              % call OutputOpti update every ItUpOut iterations
    Opt.maxiter=maxIt;                % max number of iterations
    Opt.run(zeros(Hn{1}.sizein));     % run the algorithm zeros(H.sizein)
    if isGPU==1
        xopt=gather(Opt.xopt(:,:,nbOutSl+1));
        fftxopt=gather(log(1+abs(fftshift(fftshift(Sfft(xopt,3),1),2))));   % because Sfft uses zeros_
    else
        xopt=Opt.xopt(:,:,nbOutSl+1);
        fftxopt=log(1+abs(fftshift(fftshift(Sfft(xopt,3),1),2)));
    end   
    if sav
        if ~isempty(pattname)
            saveastiff(single(xopt),[outFolder,'Recons_nbPl',num2str(nbOutSl*2-1),'_lamb',num2str(lamb(ii)),'.tif']);
            saveastiff(single(fftxopt),[outFolder,'Recons_nbPl',num2str(nbOutSl*2-1),'_lamb',num2str(lamb(ii)),'-FFT.tif']);
        else
            saveastiff(single(xopt),[outFolder,'Deconv_nbPl',num2str(nbOutSl*2-1),'_lamb',num2str(lamb(ii)),'.tif']);
            saveastiff(single(fftxopt),[outFolder,'Deconv_nbPl',num2str(nbOutSl*2-1),'_lamb',num2str(lamb(ii)),'-FFT.tif']);
        end
    end
end