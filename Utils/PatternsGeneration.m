function P=PatternsGeneration(sz,par,a,bet,ns,lamb,res,wf)
%--------------------------------------------------------------------------
% function P=PatternsGeneration(sz,par,a,bet,ns,lamb,res,wf)
%
% Generates 3D SIM illumination patterns [1]:
%        W(x,y,z) = a0 + a1 cos(kx x + ky y + phxy)cos(kz z + phz)  
%                      + a2 cos(2*[kx x + ky y + phxy])
%
% Inputs : sz    -> Size of the patterns
%          par   -> Orientation and phase parameters. Array N x 3
%                       [orr1 phXY1 phZ2
%                             ...
%                        orrN phXYN phZN]
%                   where N is the number of patterns 
%          a     -> Amplitude coefficients
%          bet   -> Angle between side beams and the optic axis
%          ns    -> Refractive index sample medium
%          lamb  -> Excitation wavelength
%          res   -> Resolution res=[resX,resY,resZ]
%          wf    -> Boolean true if widefield (constant) pattern
%
% Output : P     -> Patterns, array of size [sz,size(par,1)]
%
% References:
%    [1] Mats GL Gustafsson et al, Three-dimensional resolution doubling in wide-field
%        fluorescence microscopy by structured illumination. Biophysical journal, 2008.
%
% Emmanuel Soubies (2017) emmanuel.soubies@epfl.ch
%--------------------------------------------------------------------------

% Check
assert(size(par,2)==3,'par must be a nbPatterns x 3 matrix');

% Pattern Generation
P=zeros([sz,size(par,1)]);            % Patterns initialization
k0=2*pi*ns/lamb;                      % Amplitude wave vector
kk=k0*[sin(bet) sin(bet) 1-cos(bet)]; % Wave vector
[X,Y,Z]=meshgrid(0:sz(2)-1,0:sz(1)-1,0:sz(3)-1);
X=X*res(1);Y=Y*res(2);Z=Z*res(3);    
for it=1:size(par,1)
    P(:,:,:,it)=a(1)+ ...
        a(2)*cos(kk(1)*cos(par(it,1))*X + kk(2)*sin(par(it,1))*Y + par(it,2)).*cos(kk(3)*Z + par(it,3))+...
        a(3)*cos(2*(kk(1)*cos(par(it,1))*X + kk(2)*sin(par(it,1))*Y + par(it,2)));
end
P=P-min(P(:));
if wf
    P(:,:,:,end+1)=ones(sz);
end

% Normalization such that the mean of each pattern is 1/#Patterns
nbPatt=size(P,4);
for ii=1:nbPatt
    tmp=P(:,:,:,ii);
    P(:,:,:,ii)=P(:,:,:,ii)/(mean(tmp(:))*nbPatt);
end
end