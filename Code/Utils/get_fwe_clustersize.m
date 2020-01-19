function nMinVoxels = get_fwe_clustersize(SPM, iContrast, pClusterFormingThreshold)
% returns minimum cluster size that respect cluster level family-wise error
% (FWE) correction for multiple comparison for a given contrast and
% cluster-forming threshold
%
%  nMinVoxels = get_fwe_clustersize(SPM, iContrast, pClusterFormingThreshold)
%
%
% IN
%
% OUT
%   nMinVoxels  min. number of voxels in cluster that survives FWE
%               correction at the cluster level
% EXAMPLE
%   get_fwe_clustersize
%
%   See also

% Author:   Lars Kasper
% Created:  2019-11-14
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
pClusterFWE = 0.05;

if isfile(SPM)
    spm_mat_root = fileparts(SPM);
    SPM = load(SPM, 'SPM');
else
    spm_mat_root = SPM.swd;
end

if nargin < 3
    pClusterFormingThreshold = 0.001;
end

% setup xSPM input struct

xSPM.swd = spm_mat_root;
xSPM.title = '';
xSPM.Ic = iContrast;
xSPM.n = 1;
xSPM.Im = [];
xSPM.pm = [];
xSPM.Ex = [];
xSPM.k = 0;

secondlevelStatsThreshold = 'cluster';
switch secondlevelStatsThreshold
    case 'peak'
        xSPM.u = pClusterFormingThreshold; % in this case, u is the p value threshold - in the
        % output xSPM, u will be the stat threshold (e.g.,
        % T value)
        xSPM.thresDesc = 'FWE';
    case 'cluster'
        xSPM.u = pClusterFormingThreshold; % in this case, u is the p value threshold - in the
        % output xSPM, u will be the stat threshold (e.g.,
        % T value)
        xSPM.thresDesc = 'none';
end

[SPM, xSPM] = spm_getSPM(xSPM);

table = spm_list('Table', xSPM);

idxRowSignificantClusters = find(cell2mat(cellfun(@(x) ~isempty(x) && x<pClusterFWE, ...
    table.dat(:,3), 'UniformOutput', false)));

nVoxelsPerSignificantCluster = cell2mat(table.dat(idxRowSignificantClusters,5));

nMinVoxels = min(nVoxelsPerSignificantCluster);