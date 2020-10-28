function [patt,orr,phOff,k0]=PatternsFromFairSimXml(file,varargin)
% Generate pattern from FairSim estimated parameters.

%% Read file
st = xml2struc(file);

%% Extract parameters
band=str2double(st.fairsim.sim_dash_param.nr_dash_bands.Text);
res_xy=str2double(st.fairsim.sim_dash_param.microns_dash_per_dash_pxl.Text);   % Lateral resolution
sz=str2double(st.fairsim.sim_dash_param.img_dash_size_dash_pxl.Text);          % Image size (square)
nbOrr=str2double(st.fairsim.sim_dash_param.nr_dash_angles.Text);               % Number of orientations
nbPh=str2double(st.fairsim.sim_dash_param.nr_dash_phases.Text);                % Number of phases
nbPatt=nbOrr*nbPh;                                                             % Total number of patterns
phShift=linspace(0,2*pi,nbPh+1); phShift=phShift(1:end-1);                     % Phase shift
k0=0;
for ii=1:nbOrr
    phOff(ii)=str2double(getfield(getfield(getfield(st.fairsim.sim_dash_param,['dir_dash_',num2str(ii-1)]),'phase_dash_offset'),'Text'));
    shift=str2num(getfield(getfield(getfield(st.fairsim.sim_dash_param,['dir_dash_',num2str(ii-1)]),'shift'),'Text'));
    k0(ii)=norm(shift);
    sinA=shift(2)/norm(shift);
    cosA=shift(1)/norm(shift);
    orr(ii)=-sign(sinA)*acos(cosA);
end

if nargin==1
    szFinal=[sz,sz];
else
    szFinal=varargin{1};
end

%% Build patterns
% Optical parameters (not really used here, k0 is used)
lamb=0.488;           % Illumination wavelength 
Na=1.4;               % Objective numerica aperture
ni=1.518;             % refractive index of the objective working medium
ns=1.518;             % refractive index of the imaged medium (sample)
if band==3
    par=[orr',phOff'+phShift(1),zeros(nbOrr,1)];
    for ii=2:nbPh
        par=[par;orr',phOff'+phShift(ii),zeros(nbOrr,1)];
    end
    bet=asin((mean(k0)*lamb)/(sz*res_xy*2*ni));
    patt=squeeze(PatternsGeneration([szFinal(1)*2 szFinal(2)*2 1],par,[1 1 1],bet,ns,lamb,[res_xy/2 res_xy/2 1],0));
elseif band==2
    patt=zeros([szFinal(1)*2,szFinal(2)*2,nbPatt]);
    [X,Y]=meshgrid(0:szFinal(2)*2-1,0:szFinal(1)*2-1);X=X*res_xy/2;Y=Y*res_xy/2;
    it=1;
    for ii=1:length(orr)
        bet=asin((k0(ii)*lamb)/(sz*res_xy*2*ni));
        k=2*pi*ns/lamb*[cos(orr(ii)), sin(orr(ii))]*sin(bet);
        for jj=1:length(phShift)
            patt(:,:,it)=1+ cos(2*(k(1)*X+k(2)*Y + phOff(ii)/2+ phShift(jj)/2));
            it=it+1;
        end
    end
    patt=patt-min(patt(:))+1e-3;
    nbPatt=size(patt,3); % Normalization such that the mean of each pattern is 1/#Patterns
    for ii=1:nbPatt
        tmp=patt(:,:,ii);
        patt(:,:,ii)=patt(:,:,ii)/(mean(tmp(:))*nbPatt);
    end
end


end
