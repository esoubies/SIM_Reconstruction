classdef MyOutputOpti < OutputOptiSNR
    % MyOutputOpti : Derived class from :class:`OutputOpti` which implement
    % the SNR by doing a regression between the two images
    %
    % See also :class:`OutputOpti`
    
    %%     Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.    
    
    methods
    	%% Constructor
        function this=MyOutputOpti(computecost,xtrue,iterVerb,costIndex) 
            if nargin <4, costIndex=[]; end
            if nargin <3, costIndex=[]; iterVerb=[]; end
            if nargin <2, costIndex=[]; iterVerb=[]; xtrue=[]; end
            if nargin <1, costIndex=[]; iterVerb=[]; xtrue=[]; computecost=[]; end
            this@OutputOptiSNR(computecost,xtrue,iterVerb,costIndex) ;
        end
        %% Update method
        function update(this,opti)
            % Computes SNR, cost and display evolution.
            str=sprintf('Iter: %5i',opti.niter);
            if this.computecost
                cc=this.computeCost(opti);
                str=sprintf('%s | Cost: %4.4e',str,cc);
                this.evolcost(this.count)=cc;
            end
            
            if ~isempty(this.xtrue)
                snr=this.computeSNR(opti);
                str=sprintf('%s | SNR: %4.4e dB',str,snr);
                this.evolsnr(this.count)=snr;
            end
            
            if this.saveXopt
                this.evolxopt{this.count}=opti.xopt;
            end
            this.iternum(this.count)=opti.niter;
            this.count=this.count+1;
            if opti.verbose && (opti.niter~=0 && (mod(opti.niter,this.iterVerb)==0) || (opti.niter==1 && this.iterVerb~=0))
                disp(str);
            end
        end
        function snr=computeSNR(this,opti)
            % Reimplemented from :class:`OutputOpti`
            snr=OptSNR(this.xtrue,opti.xopt);
        end
    end
end
