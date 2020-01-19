function scanInfo = spifi_get_scaninfo(idxSess)
% returns info about a certain scan (fmri-run), e.g., nSlices, resolution
%
%   scanInfo = spifi_get_scaninfo(idxSess)
%
% IN
%
% OUT
%
% EXAMPLE
%   spifi_get_scaninfo
%
%   See also
 
% Author:   Lars Kasper
% Created:  2019-06-16
% Copyright (C) 2019 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3. 
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
%

switch idxSess
    case {1,4}
        scanInfo.resolution_mm = [190/232 230/280 1];
    case {2,3, 5, 6}
        scanInfo.resolution_mm = [190/160 230/192 1];
end

switch idxSess
    case {1,2}
        scanInfo.TR = 3.306;
    case {3}
        scanInfo.TR = 3.126;
    case {4}
        scanInfo.TR = 2.506;
end
scanInfo.nSlices = 36;
