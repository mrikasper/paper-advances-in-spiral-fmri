function [subjectM, mniM] = warp_visual_cortex_mask_to_subject(idxSubj)
% Warps mask of early visual cortex (human OC V1-3 ventral/dorsal) to
% individual subject space using inverse deformation field from
% Segmentation
%
%   [M, warpedM] = warp_visual_cortex_mask_to_subject(idxSubj)
%
% IN
%   
% OUT
%   mniM         mask in MNI space
%   subjectM   warped mask in subject space
%
% EXAMPLE
%   warp_visual_cortex_mask_to_subject
%
%   See also

% Author:   Lars Kasper
% Created:  2021-09-20
% Copyright (C) 2021 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
%
if nargin < 1
    idxSubj = 7;
end

thresholdMask = 0.6; % TODO: make parameter in options.preproc.warp2subject.thresholdMask
nVoxelsDilation = 3; % dilation of subject mask to be more inclusive of voxels; set to 0 or [] to skip this step
thresholdImage = [0 0.8];
selectedSlicesMni = [79:114]-2;
selectedSlicesSubject = 1:36;

options = spifi_get_analysis_options();
details = spifi_get_subject_details(idxSubj, options);

fileRoiArray = options.paths.data.spm_anatomy_toolbox.maps_early_visual;
fileUnderlayMni = details.preproc.warp2mni.out{1};
fileUnderlaySubject = details.preproc.warp2mni.source{1};
fileInverseDeformationField = details.preproc.anat.inverseDeformationField;

doPlot = 1; % for plots, 2 for debug plots
% if true, use GUI export from SPM Anatomy toolbox (recommended!)
% if false, recalculation within this function using new thresholds; note
% that there is a shift between "anatomical MNI" (PMaps are saved in) and SPM's MNI space
% which hasn't been corrected for the re-calculation version yet
doUseExportedMaskFromToolbox = true; 


if doPlot
    Y = MrImage(fileUnderlayMni);
end

if doUseExportedMaskFromToolbox
    mniM = MrImageSpm4D(options.paths.data.spm_anatomy_toolbox.mask_visual_cortex);
else
    % load probability map data from different visual cortex regions
    %, then visualize, threshold and mask
    
    nRois = numel(options.paths.data.spm_anatomy_toolbox.maps_early_visual);
    for iRoi = 1:nRois
        X{iRoi} = MrImageSpm4D(fileRoiArray{iRoi});
        
        if doPlot >= 2
            X{iRoi}.plot();
        end
        
        M{iRoi} = X{iRoi}.binarize(thresholdMask);
    end
    
    % plot overlay with mean func
    if doPlot
        Y.threshold(thresholdImage).plot_overlays(X, 'overlayMode', 'mask', ...
            'selectedSlices', selectedSlicesMni, 'colorBar', 'off', ...
            'overlayAlpha', 0.5)
        set(gcf, 'Name', 'MNI space - Mean Functional with HOC V1-V3 mask')
    end
    
    % Sum of all probabblity maps and mask at threshold
    mniM = X{1}.copyobj;
    for iRoi = 2:nRois
        mniM = mniM + X{iRoi};
    end
    
    mniM = mniM.binarize(thresholdMask);
end

if doPlot >=2
    mniM.plot('z', selectedSlicesMni);
end

mniM.save('fileName', details.preproc.warp2subject.mask_early_visual.source);
if doPlot
    Y.threshold(thresholdImage).plot_overlays(mniM, ...
        'overlayMode', 'mask', ...
        'selectedSlices', selectedSlicesMni, 'colorBar', 'off', ...
        'overlayAlpha', 0.5);
        set(gcf, 'Name', 'MNI space - Mean Functional with Combined mask')
end

%% Warp to subject space
Z = MrImage(fileUnderlaySubject); % important for dimension of mask!


inverseDeformationField = MrImage(fileInverseDeformationField);
subjectM = mniM.apply_deformation_field(inverseDeformationField, ...
    'voxelSize', Z.dimInfo.resolutions);
subjectM = subjectM.reslice(Z);
subjectM = subjectM.binarize(thresholdMask);

if ~isempty(nVoxelsDilation) && nVoxelsDilation > 0
   subjectM = subjectM.imdilate(strel('disk', nVoxelsDilation));
end

subjectM.save('fileName', details.preproc.warp2subject.mask_early_visual.out);

if doPlot
   Z.threshold(thresholdImage).plot_overlays(subjectM, ...
        'overlayMode', 'mask', ...
        'selectedSlices', selectedSlicesSubject, 'colorBar', 'off', ...
        'overlayAlpha', 0.5);
        set(gcf, 'Name', sprintf('Subject %d space - Mean Functional with Combined mask', idxSubj))
end