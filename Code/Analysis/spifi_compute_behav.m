function spifi_compute_behav(fileBehav, fileMultipleConditions, scanInfo)
% Converts behavioral .log file of vis quarterfield paradigm to multiple
% conditions .mat file
%
%    spifi_compute_behav(fileBehav, fileMultipleConditions, scanInfo);
%
% IN
%   scanInfo    to adjust onsets for TR delay of trigger
%
% OUT
%
% EXAMPLE
%   spifi_compute_behav
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

[onsets, durations, names] = ...
    get_multiple_conditions_visual(fileBehav, 'fileMultipleConditions', ...
    fileMultipleConditions, 'isVerbose', true);

% correction for extra TR in the beginning!
onsets = cellfun(@(x) x - scanInfo.TR, onsets, 'UniformOutput', false);
save(fileMultipleConditions, 'onsets', 'names', 'durations');
