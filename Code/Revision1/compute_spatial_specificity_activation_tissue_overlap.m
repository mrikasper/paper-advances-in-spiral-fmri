function [fhArray, absY, figOptions] = ...
    compute_spatial_specificity_activation_tissue_overlap(varargin)
% Function compute_spatial_specificity_activation_tissue_overlap
% Following reviewer's suggestions, quantify spatial accuracy via overlap
% between activated voxels and gray matter tissue
%
%   [fhArray, Y, figOptions] = compute_spatial_specificity_contour_congruency(figOptions)
%
% DETAILS
% Using UniQC and SPM, perform
% 1. unified segmentations of multi-echo structural images
% (co-registered to functional data)
% 2. Extract the GW/WM/CSF and border ROIs by thresholding TPMs
% 3. Compute ROI distribution in contrast images
% 4. Provide histograms
%
% OUT
%   fhArray     created figure handles
%   absY        absolute value of statistical map with all tissue ROIs
%   See also main_create_figures

% Author:   Lars Kasper
% Created:  2020-08-13
% Copyright (C) 2020 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaults.doUseSmoothedMaps = false;
defaults.doUseVisualCortexMask = false;
% if true, threshold is used to consider only voxels above a certain tissue probability before edge detection
defaults.thresholdTPM = 0.9;
% individual probability thresholds for boundaries between 2 tissue
% classes, indicating uncertainty about the class by having individual
% probabilities above this threshold
defaults.thresholdTPMBoundary = 0.3;
defaults.thresholdRange = [3.22 6]; %[3.2,8]; % T or F value % T = 3.22 <=> p <= 0.001
defaults.FWEcorrection = 'cluster'; % 'none'; 'cluster', 'peak'; %TODO: peak!
defaults.doRecalcTPMs = false;
defaults.sourceTPMs = 'meanME'; % 'meanME' or 'T1' or echo number (1-6)
defaults.displayRangeStruct = [];
defaults.displayRangeFunc = [0 0.8];
defaults.selectedSlice = 3:4:34; %27
defaults.idxContrast = 1;
defaults.verbosity_level = 0;
defaults.iSubj = 3;
defaults.iSess = 1;
defaults.iRecon = 1;

% merge with input fields
figOptions = propval(varargin, defaults);
strip_fields(figOptions);

% sharedParameters = {'z', selectedSlice, 'rotate90', 1, 'colorBar', 'off'};
sharedParametersOverlays = {'selectedSlices', selectedSlice, 'rotate90', 1, ...
    'colorBar', 'off', 'overlayColorMaps', {'winter','autumn'}};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = spifi_get_analysis_options();
options.glm.doUseSmoothed = doUseSmoothedMaps;
options.representation.fig_specificity_contours.doUseVisualCortexMask = doUseVisualCortexMask;

options.representation.fig_specificity_contours.struct_select = sourceTPMs;

if isempty(displayRangeStruct)
    switch upper(sourceTPMs)
        case 'T1'
            displayRangeStruct = [0 1000];
        otherwise % ME, any echo or mean
            displayRangeStruct = [0 1.2];
    end
end

details = spifi_get_subject_details(figOptions.iSubj, options, ...
    figOptions.iSess, figOptions.iRecon);


figOptionsFromContours = details.representation.fig_specificity_contours;

if doUseVisualCortexMask
    fileMaskVisual = figOptionsFromContours.mask_visual;
    maskVisualCortex = MrImageSpm4D(fileMaskVisual);
end

fileStruct = figOptionsFromContours.struct_biascorrected;
switch idxContrast
    case 1
        fileStatMap = details.glm.tcon;
    otherwise
        error('Check how to load other contrasts');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data and initial plot, UniQC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = MrImageSpm4D(fileStruct);
Y = MrImageSpm4D(fileStatMap);

% resample X to higher res cunctional
rX = X.resize(Y.dimInfo);

fhArray = zeros(0,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocess data: Segment, threshold, edge detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Segment both images independently, retrieve tissue probability maps (TPMs)
for iTissue = 1:3
    structTPMs{iTissue} = MrImage(figOptionsFromContours.structTPMs{iTissue});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% two-sided T-test, use absolute values above minimum threshold
absY = abs(Y);
maskActivation = absY.binarize(thresholdRange(1));

% remove small clusters (multiple comparison correction)
switch FWEcorrection
    case 'cluster'
        % remove small clusters following SPM's FWE-cluster correction
        nMinVoxels = get_fwe_clustersize(details.glm.spm_mat, idxContrast);
        maskActivation = maskActivation.remove_clusters(...
            'nPixelRange', [1 nMinVoxels-1]);
    case 'peak'
        % TODO
    case 'none'
        % nothing to remove
end

% remove all voxels outside visual cortex
if doUseVisualCortexMask
    maskActivation = maskActivation.*maskVisualCortex;
end


% Canny uses local intensity gradient maxima, seems to work better than global (Sobel)
maskGM = structTPMs{1}.binarize(thresholdTPM);
maskWM = structTPMs{2}.binarize(thresholdTPM);
maskCSF = structTPMs{3}.binarize(thresholdTPM);

% boundary defined by min. individual probabilities for two classes,
% only consider voxels that have not been classified before
maskUnidentified = (maskWM.*-1 + 1).*(maskGM.*-1 + 1).*(maskCSF.*-1 + 1);

maskWMSurface = maskUnidentified.*(structTPMs{1}.binarize(thresholdTPMBoundary)).*(structTPMs{2}.binarize(thresholdTPMBoundary));
maskUnidentified = maskUnidentified.*(maskWMSurface.*-1 + 1);

maskPialSurface =  maskUnidentified.*(structTPMs{1}.binarize(thresholdTPMBoundary)).*(structTPMs{3}.binarize(thresholdTPMBoundary));
maskUnidentified =  maskUnidentified.*(maskPialSurface.*-1 + 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Rois as sanity check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbosity_level >=2
    fhArray(end+1,1) = rX.threshold(displayRangeStruct).plot_overlays({Y, Y.*-1}, 'overlayThreshold', thresholdRange, ...
        'overlayMode', 'map', sharedParametersOverlays{:});
    set(fhArray(end), 'Name', 'Structural Image with significantly activated voxels');
end

% FWE-corrected activation on structural overlay
if verbosity_level >=1
    fhArray(end+1,1) = rX.threshold(displayRangeStruct).plot_overlays(...
        {Y.*maskActivation, Y.*-1.*maskActivation}, 'overlayThreshold', thresholdRange, ...
        'overlayMode', 'map', sharedParametersOverlays{:});
    set(fhArray(end), 'Name', 'Structural Image with significantly activated voxels (FWE-corrected)');
end

maskPlotArray = {maskActivation, maskGM, maskWM, maskCSF, maskUnidentified, ...
    maskPialSurface, maskWMSurface};

roiNameArray = {'Activated Voxels', 'GM', 'WM', 'CSF', ...
    'Unidentified/Ambiguous', 'Pial Surface', ...
    'WM/GM Surface'};

for iRoi  = 1:numel(maskPlotArray)
    maskPlotArray{iRoi}.name = roiNameArray{iRoi};
    
    if verbosity_level >= 2
        fhArray(end+1,1) = maskPlotArray{iRoi}.plot_overlays({Y, Y.*-1}, 'overlayThreshold', thresholdRange, ...
            'overlayMode', 'map', sharedParametersOverlays{:});
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract ROI stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Limit Rois to activated voxels only
for iRoi = 1:numel(maskPlotArray)
    maskRois{iRoi} = maskPlotArray{iRoi}.*maskActivation;
end

% plot activation map with all ROIs as edges
if verbosity_level >= 3
    fhArray(end+1,1) = absY.threshold([0, thresholdRange(2)]).plot_overlays(...
        maskPlotArray, 'overlayMode', 'edge', sharedParametersOverlays{:});
    set(fhArray(end), 'Name', 'Activation map with Edges of all tissue ROIs');
end

if verbosity_level >=2
    % plot structural image with all ROIs
    fhArray(end+1,1) = rX.threshold(displayRangeStruct).plot_overlays(maskPlotArray, ...
        'overlayMode', 'mask', sharedParametersOverlays{:}, 'overlayAlpha', 0.5);
    set(fhArray(end), 'Name', 'Structural Image with Overlay of all tissue ROIs');
end

absY.extract_rois(maskRois);
absY.compute_roi_stats();

if verbosity_level >=2
    for iRoi = 1:numel(maskRois)
        fhArray(end+1,1) = absY.rois{iRoi}.plot('axisType', 'relative', 'selectedSlices', selectedSlice);
        % set(gcf, 'Name', sprintf('ROI: %s', roiNameArray{iRoi}));
    end
end

absY.save('fileName', ...
    details.representation.fig_specificity_activation_tissue_overlap.output_spm_with_tissue_rois);

% include smoothing status and roi setting in figure names
for iFig = 1:numel(fhArray)
    set(fhArray(iFig), 'Name', ...
        sprintf('%s - doUseSmoothedMaps%d - doUseVisualCortexMask%d', ...
        get(fhArray(iFig), 'Name'), ...
        figOptions.doUseSmoothedMaps, ...
        figOptions.doUseVisualCortexMask));
end