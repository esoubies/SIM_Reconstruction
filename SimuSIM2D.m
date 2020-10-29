%----------------------------------------------------
% Script to simulate SIM data. Illumination patterns are generated using:
%
%  W(x,y) = 1 + a cos(2*(kx x + ky y + ph))
%
% where a is the amplitude of the modulation, ph the phase of the
% illumination pattern and k=[kx,ky] is the wave vector given by
%   k=(2 pi ns)/lamb [cos(th), sin(th)]*sin(bet) 
% with th the lateral orientation of the illumination, nl the refractive
% index of the objective medium and bet the angle between side beams and
% the optic axis (of the two interfering beams used to generate the
% pattern) which should be lower than asin(NA/nl).
%
% If the given image is 3D, the script perform a 3D acquisition and returns
% only the focal plane.
%
% Before to run it set the following papareters:
% -- General 
% sav=             -> to save results
% expFolder=       -> experiment folder
% gtpath=          -> file name ground truth 
% 
% -- PSF
% lamb=            -> Illumination wavelength
% res=             -> Resolution
% Na=              -> Objective numerica aperture
% nl=              -> Refractive index of the objective medium (glass/oil)
% ns=              -> Refractive index of the sample medium (water)
%
% -- Patterns
% orr=             -> Patterns orientations (vector)
% ph=              -> Patterns lateral phases (vector)
% a=               -> Amplitude coefficient
% bet=             -> Angle between side beams and the optic axis (e.g. bet asin(Na/nl))
%
% -- Acquisition
% downFact=        -> Downsmpling factor (e.g. [2 2])
% photBud=         -> Photon Budget (average per voxel of ground truth)
%
% Copyright (2018) Emmanuel Soubies, emmanuel.soubies@epfl.ch
%----------------------------------------------------
addpath ./Utils 
javaaddpath ./Utils/PSFGenerator.jar

%% Reading data
fprintf('Reading data .............');
im=double(loadtiff(gtpath)); im=im/max(im(:));
sz=size(im);
fprintf(' done \n');

%% PSF Generation
fprintf('PSF Generation ...........');
fc=2*Na/lamb*res;
if length(sz)==3
    setConfigFilePSF(sz(2),sz(1),sz(3),res,resZ,Na,lamb,'./Utils/configPSF.txt');
    psf = PSFGenerator.compute('./Utils/configPSF.txt');psf=psf/sum(psf(:));
    ll=linspace(-0.5,0,sz(1)/2+1);
    lr=linspace(0,0.5,sz(1)/2);
    [X,Y]=meshgrid([ll,lr(2:end)],[ll,lr(2:end)]);
    [th,rho]=cart2pol(X,Y);
    OTF=fftn(fftshift(psf));OTF=fftshift(fftshift(OTF).*repmat((rho<fc),[1 1 size(OTF,3)]));psf=real(fftshift(ifftn(OTF)));
    figure;subplot(1,2,1);imagesc(log(1+abs(fftshift(OTF(:,:,1))))); axis image; title('OTF'); colormap(fire(200));viscircles(floor(sz(1:2)/2)+1,fc*sz(1));caxis([0 0.01*max(abs(OTF(:)))]);
    subplot(1,2,2);imagesc(psf(:,:,floor(sz(3)/2)+1)); axis image; title('PSF'); caxis([0 0.01*max(psf(:))]);
else   
    ll=linspace(-0.5,0,sz(1)/2+1);
    lr=linspace(0,0.5,sz(1)/2);
    [X,Y]=meshgrid([ll,lr(2:end)],[ll,lr(2:end)]);
    [th,rho]=cart2pol(X,Y);
    OTF=fftshift(1/pi*(2*acos(abs(rho)/fc)-sin(2*acos(abs(rho)/fc))).*(rho<fc));
    figure;subplot(1,2,1);imagesc((fftshift(OTF))); axis image; title('OTF'); colormap(fire(200));viscircles(floor(sz/2)+1,fc*sz(1));
    psf=real(fftshift(ifft2(OTF)));
    subplot(1,2,2);imagesc(psf); axis image; title('PSF'); caxis([0 0.01*max(psf(:))]);
end
fprintf(' done \n');

%% Patterns Generation
fprintf('Patterns Generation ......');
patt=zeros([sz(1:2),length(orr)*length(ph)]);
[X,Y]=meshgrid(0:sz(2)-1,0:sz(1)-1);X=X*res;Y=Y*res;
it=1;
for ii=1:length(orr)
    k=2*pi*ns/lamb*[cos(orr(ii)), sin(orr(ii))]*sin(bet);
    for jj=1:length(ph)
        patt(:,:,it)=1+ a*cos(2*(k(1)*X+k(2)*Y + ph(jj)));
        it=it+1;
    end
end
nbPatt=size(patt,3); % Normalization such that the mean of each pattern is 1/#Patterns
for ii=1:nbPatt
    tmp=patt(:,:,ii);
    patt(:,:,ii)=patt(:,:,ii)/(mean(tmp(:))*nbPatt);
end
figure;subplot(1,2,1);imagesc(patt(:,:,1)); axis image; title('Example pattern');
subplot(1,2,2);imagesc(log(1+abs(fftshift(fftn(patt(:,:,1)))))); axis image; title('Example pattern FFT'); 
viscircles(floor(sz(1:2)/2)+1,fc*sz(1));
fprintf(' done \n');

%% Image Noiseless Acquisition
fprintf('Acquisition simulation ...');
% - LinOp Downsampling and integration over camera pixels
if length(sz)==3
    SS=LinOpSum(sz,3);
else
    SS=LinOpIdentity(sz);
end
S=LinOpDownsample(sz(1:2),downFact);
% - LinOpConv (PSF)
OTF=Sfft(fftshift(fftshift(psf(:,:,end:-1:1),1),2),3);
H=LinOpConv(OTF,1,[1 2]);
% - Acquisition
acqNoNoise=zeros([S.sizeout,size(patt,3)]);
fprintf(' Pattern # ');
for it=1:size(patt,3)
    fprintf([num2str(it),'..']);
    D=LinOpDiag(sz,patt(:,:,it));
    acqNoNoise(:,:,it)=S*SS*H*D*im;
end
fprintf(' done \n');
    
%% Add noise and Save
for ii=1:length(photBud)
    % - Add Noise
    acq=acqNoNoise;
    if photBud>0
        tmp=sum(acqNoNoise,3);
        factor = photBud(ii)./mean(tmp(:)) ;
        acqNoNoise = acqNoNoise.* factor;
        im = im.*factor;
        acq = random('Poisson',acqNoNoise);
    end
    SNR=20*log10(norm(acqNoNoise(:))/norm(acq(:)-acqNoNoise(:)));
    disp(['SNR = ',num2str(SNR),' dB']);
    
    % - Save
    if sav
        saveastiff(single(acqNoNoise),[expFolder,'/AcqDataNoiseless.tif']);
        saveastiff(single(log(1+abs(fftshift(fftshift(Sfft(acqNoNoise,3),1),2)))),[expFolder,'/AcqDataNoiseless-FFT.tif']);
        saveastiff(single(acq),[expFolder,'/AcqData.tif']);
        saveastiff(single(log(1+abs(fftshift(fftshift(Sfft(acq,3),1),2)))),[expFolder,'/AcqData-FFT.tif']);
        saveastiff(single(sum(acq,3)),[expFolder,'/WFData.tif']);
        saveastiff(single(log(1+abs(fftshift(fft2(sum(acq,3)))))),[expFolder,'/WFData-FFT.tif']);
        saveastiff(single(sum(acqNoNoise,3)),[expFolder,'/WFDataNoiseless.tif']);
        saveastiff(single(log(1+abs(fftshift(fft2(sum(acqNoNoise,3)))))),[expFolder,'/WFDataNoiseless-FFT.tif']);
    end
end
save([expFolder,'/PSF'],'psf');
saveastiff(single(patt),[expFolder,'/patterns.tif']);

figure;subplot(1,2,1);imagesc(acq(:,:,1)); axis image; title('Example Acquired image');
subplot(1,2,2);imagesc(log(1+abs(fftshift(fftn(acq(:,:,1)))))); axis image; title('Example Acquired image FFT'); 
figure;subplot(1,2,1);imagesc(sum(acq,3)); axis image; title('Widefield image');
subplot(1,2,2);imagesc(log(1+abs(fftshift(fftn(sum(acq,3)))))); axis image; title('Widefield image FFT'); 



