classdef OutputOpti2DSIM < OutputOpti
    % OutputOpti2DSIM
    %
    % :param fp: focal plane to extract from xopt to compute the SNR
    %
    % Extract the focal plane of the optimized variable to compute the SNR
    %
    % See also :class:`Opti` :class:`OutputOpti`  :class:`OutputOptiADMM`

    %%    Copyright (C) 2018
    %     E. Soubies emmanuel.soubies@epfl.ch
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
    
    properties (SetAccess = public,GetAccess = public)
        fp;                % Focal plane to use for computing SNR
    end
    
    methods
        %% Constructor
        function this=OutputOpti2DSIM(computecost,iterVerb,costIndex)
             this@OutputOpti(computecost,iterVerb,costIndex);
        end
        function snr=computeSNR(this,opti)
            % Reimplemented from parent class :class:`OutputOptiADMM`
            % Extract the focal plane (attribute fp) of the optimized
            % variable xopt to compute the snr.
            tmp=opti.xopt(:,:,this.fp);
            snr=OptSNR(this.xtrue,tmp);
        end
    end
end