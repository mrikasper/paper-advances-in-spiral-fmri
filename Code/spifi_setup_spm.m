function spifi_setup_spm()
% Set up SPM with better defaults for fMRI
%
%   spifi_setup_spm(input)
%
% IN
%
% OUT
%
% EXAMPLE
%   spifi_setup_spm
%
%   See also
%
% Author:   Lars Kasper
% Created:  2018-03-09
% Copyright (C) 2018 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the Zurich fMRI Methods Evaluation Repository, which is released
% under the terms of the GNU General Public License (GPL), version 3. 
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.

spm('defaults','FMRI');
spm_get_defaults('stats.maxmem', 4e9);
spm_get_defaults('stats.resmem', true); % read all images at once
spm_get_defaults('mat.format', '-v7.3'); % enable saving files >4GB

% make masks great again!
spm_get_defaults('mask.thresh', 0.2);

spm_jobman('initcfg');