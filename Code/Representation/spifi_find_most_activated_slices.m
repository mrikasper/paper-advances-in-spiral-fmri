function [idxSlice, centre] = spifi_find_most_activated_slices(statMap, threshold)
% for each ortogonal view (sagittal, coronal, transversal), find slices
% with most activated voxels for given contrast map and threshold
%
%   [idxSlice, centre] = spifi_find_most_activated_slices(fileStatMap, threshold)
%
%
% IN
%   statMap         string, filename (with path) 
                        % OR
%                   MrImage
%                   of statistical map to be
%                   thresholded; For T-maps, absolute values are considered
%   threshold       statistical threshold which is applied to map before
%                   counting active voxels
% OUT
%
% EXAMPLE
%   spifi_find_most_activated_slices
%
%   See also
 
% Author:   Lars Kasper
% Created:  2019-09-06
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
% determine slices (for all 3 orientations) with most significant voxels from combined contrast
    for iDim = 1:3 
        if isa(statMap, 'MrImage')
            X = statMap.copyobj;
        else
            X = MrImage(statMap);
        end
        % permute data to put slice dimensions to 3rd dimension
        switch iDim
            case 1
                X = X.permute([2 3 1]);
            case 2
                X = X.permute([1 3 2]);
        end
        X.dimInfo.dimLabels = {'x','y','z'}; % hack, because roi extraction always along z dim
        M = X.copyobj.abs.binarize(threshold); % note: for t-maps, this takes positive and negative peaks into account
        X.extract_rois(M);
        %X.compute_roi_stats(); % not necessary here, maybe slow
        [nVoxels(iDim), idxSlice(iDim)] = max(X.rois{1}.perSlice.nVoxels);
    end
    centre = X.dimInfo.index2sample(idxSlice,[1:3]);
    