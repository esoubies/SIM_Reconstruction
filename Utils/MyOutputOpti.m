classdef MyOutputOpti < OutputOpti
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
            if nargin <4, costIndex=[]; end;
            if nargin <3, costIndex=[]; iterVerb=[]; end;
            if nargin <2, costIndex=[]; iterVerb=[]; xtrue=[]; end;
            if nargin <1, costIndex=[]; iterVerb=[]; xtrue=[]; computecost=[]; end;
            this@OutputOpti(computecost,xtrue,iterVerb,costIndex) ;             
        end
        %% Update method
        function snr=computeSNR(this,opti)
            % Reimplemented from :class:`OutputOpti`
            snr=OptSNR(this.xtrue,opti.xopt);   
        end
    end
end
