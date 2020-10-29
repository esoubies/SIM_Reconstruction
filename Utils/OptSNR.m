function snr=OptSNR(xref,x,mask)
% function snr=OptSNR(xref,x)
% Computes an optimized SNR between xref (reference image)
% and x. The returned value maximizes
% 10*log10(||xref||/(||xref - a*x + b||))
% with respect to scalars a and b
if nargin<3 || isempty(mask), mask=ones(size(x))>0; end
    N=numel(x);
    sumref=sum(xref(:));
    sumx=sum(x(:));
    dotprod=dot(xref(:),x(:));
    a=(N*dotprod-sumref*sumx)/(N*norm(x(:))^2-sumx^2);
    b=(a*sumx-sumref)/N;
    xx=(a*x(mask>0)-b);
    rr=xref(mask>0);
    snr=20*log10(norm(rr(:))/norm(rr(:)-xx(:)));
end