function [fhArray, distMapX, figOptions] = ...
    compute_spatial_specificity_contour_congruency(figOptions)
% Function compute_spatial_specificity_contour_congruency
% Following reviewer's suggestions, quantify spatial accuracy via contour
% line matches between functional and structural data
%
%   distMapX = compute_spatial_specificity_contour_congruency(figOptions)
%
% DETAILS
% Using UniQC and SPM, perform
% 1. separate unified segmentations of mean
% spiral functional and multi-echo structural images (already
% co-registered)
% 2. Compute contours via edge() command
% 3. Per voxel, compute min. distance of within-slice contour of other
% image
% 4. Provide histograms and "heat maps" of distance
%
% OUT
%   fhArray     created figure handles
%   distMapX    same dimension as structural image, contains tissue edges
%               from structural image-derived TPM, intensity values on edge reflect
%               distance to closest in-plane edge of functional TPM
%           .rois{1}
%               ROI statistics of all edge voxels, summary of distances for
%               all edges (per slice and volume)
%   figOptions  used figure options
%               Note: contains calculated structTPMs and funcTPMs for later
%               reruns of this function
%
%
%   See also main_create_figures create_figure_contour_congruency_summary_all_subjects

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

% if true, threshold is used to consider only voxels above a certain tissue probability before edge detection

defaults.doThresholdTPM = true;
defaults.doPlotInitContours = false;% plot classical sobel-contours as in other paper figures;
defaults.thresholdTPM = 0.9;
defaults.doRecalcTPMs = true;
defaults.doRecalcVisualCortexMask = true;
defaults.displayRangeStruct = [0 1.2];
defaults.displayRangeFunc = [0 0.8];
% Use certain tissue for comparison
% 1=GM, 2=WM, 3=CSF
defaults.idxTissue = 1;
defaults.structTPMs = []; % tissue probability maps (TPMs) of structural image
defaults.funcTPMs = []; % TPM of functional image
defaults.doPlotInitContours = false;
defaults.selectedSlice = 3:4:34;
defaults.iSubj = 7;
defaults.iSess = 1;
defaults.iRecon = 1;
defaults.struct_select = 'meanME';
defaults.verbosity_level = 1; % 1 = plots for figure 5(9), 2 = diagnostic plots
defaults.doUseFiguresForSaving = true;
defaults.mask_visual = '';
defaults.struct = '';
defaults.struct_biascorrected = '';
defaults.func_biascorrected = '';
defaults.output_edge_distance_map = 'EdgeDistanceMap.nii';
defaults.func = '';
defaults.coreg = '';
% merge with input fields
figOptions = propval(figOptions, defaults);
strip_fields(figOptions);

if numel(selectedSlice) <= 4
    sharedParameters = {'z', selectedSlice, 'rotate90', 1, 'colorBar', 'off', ...
        'plotLabels', false, 'nRows', 1};
    sharedParametersOverlays = {'selectedSlices', selectedSlice, 'rotate90', 1, ...
        'colorBar', 'off', 'plotLabels', false, 'nRows', 1};
else
    sharedParameters = {'z', selectedSlice, 'rotate90', 1, 'colorBar', 'off', ...
        'plotLabels', false};
    sharedParametersOverlays = {'selectedSlices', selectedSlice, 'rotate90', 1, ...
        'colorBar', 'off', 'plotLabels', false};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if doRecalcVisualCortexMask
    warp_visual_cortex_mask_to_subject(figOptions.iSubj);
end

% options:
%   structural:     mean of ME
%                   last echo of ME
%                   MPRAGE (T1 w/o monitoring)
%   functional:      mean realigned (without bias correction for segmentation, but with for edge detection?)
%                   ...definitely as underlay!

fileMaskVisual = figOptions.mask_visual;
fileStruct = figOptions.struct;
fileFunc = figOptions.func;
fileStructTPMs = figOptions.structTPMs;
fileFuncTPMs = figOptions.funcTPMs;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocess data: Segment, threshold, edge detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Segment both images independently, retrieve tissue probability maps (TPMs)
if doRecalcTPMs || ~isfile(fileStructTPMs{idxTissue}) || ~isfile(fileFuncTPMs{idxTissue});
    
    % load data
    X = MrImageSpm4D(fileStruct);
    Y = MrImageSpm4D(fileFunc);
    
    % select mean multi-echo or specific one
    switch(figOptions.struct_select)
        case {'mean', 'meanME'}
            X = mean(X,'t');
        case {1,2,3,4,5,6}
            X = X.select('t', figOptions.struct_select);
        case 't1'
            % all good, already 3D
    end
    
    % resize to higher resolution (functional!)
    rX = X.resize(Y.dimInfo);
    
    % coregister structural to functional
    coregX = rX.copyobj();
    coregX.coregister_to(Y);
    coregX.save('fileName', figOptions.coreg.out);
    
    % segment both structural and functional and save!
    [~, structTPMs, ~, xBiasField] = coregX.segment();
    [~, funcTPMs, ~, yBiasField] = Y.segment();
    
    for iTissue = 1:3
        structTPMs{iTissue}.save('fileName', figOptions.structTPMs{iTissue});
        funcTPMs{iTissue}.save('fileName', figOptions.funcTPMs{iTissue});
    end
    
    % bias-correct and save! Use same biasfield from structural for both
    % images
    biasCorrectedX = coregX.*xBiasField;
    biasCorrectedY = Y.*xBiasField;
    
    biasCorrectedX.save('fileName', figOptions.struct_biascorrected);
    biasCorrectedY.save('fileName', figOptions.func_biascorrected);
    
end

inputEdgeReference = MrImageSpm4D(fileStructTPMs{idxTissue});
inputEdgeComparison = MrImageSpm4D(fileFuncTPMs{idxTissue});

%% Compute contours from (thresholded) TPMs

% Canny uses local intensity gradient maxima, seems to work better than global (Sobel)
if doThresholdTPM
    edgeX = inputEdgeReference.threshold(thresholdTPM).edge('Canny'); % consider only voxels that are most likely to be GM
    edgeY = inputEdgeComparison.threshold(thresholdTPM).edge('Canny'); % consider only voxels that are most likely to be GM
else % maybe less arbitrary to work w/o threshold
    edgeX = inputEdgeReference.edge('Canny');
    edgeY = inputEdgeComparison.edge('Canny');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot resulting edges on different overlays...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fhArray = zeros(0,1);

X = MrImage(figOptions.struct_biascorrected);
Y = MrImage(figOptions.func_biascorrected);

% the classical contours we plot
if doPlotInitContours
    Y.threshold(displayRangeFunc).plot_overlays(...
        X, 'overlayMode', 'edge', sharedParametersOverlays{:});
    fhArray(end+1,1) = gcf;
    X.edge('sobel',0.05).plot(sharedParametersOverlays{:});
    fhArray(end+1,1) = gcf;
    Y.edge('sobel',0.06).plot(sharedParametersOverlays{:});
    fhArray(end+1,1) = gcf;
end

% plot input images
if figOptions.verbosity_level > 1
    X.threshold(displayRangeStruct).plot(sharedParameters{:});
    fhArray(end+1,1) = gcf;
    Y.threshold(displayRangeFunc).plot(sharedParameters{:});
    fhArray(end+1,1) = gcf;
end

% overlay image with contour
X.threshold(displayRangeStruct).plot_overlays(edgeX, sharedParametersOverlays{:})
fhArray(end+1,1) = gcf;
if doUseFiguresForSaving
    axis off; title('');
end
Y.threshold(displayRangeFunc).plot_overlays(edgeY, sharedParametersOverlays{:})
fhArray(end+1,1) = gcf;
if doUseFiguresForSaving
    axis off; title('');
end

% overlay both contour sets
if figOptions.verbosity_level > 1
    edgeX.plot_overlays(edgeY, sharedParametersOverlays{:})
    fhArray(end+1,1) = gcf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate edge distance heat map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSlices = edgeX.dimInfo.z.nSamples;

distMapX = edgeX.copyobj;
distMapX.name = 'Edge Distance Map';
distMapX.data = -1*ones(size(distMapX.data));


for iSlice = 1:nSlices
    [xX,yX] = find(edgeX.data(:,:,iSlice));
    [xY,yY] = find(edgeY.data(:,:,iSlice));
    
    
    [nn, distnn] = dsearchn([xY,yY],[xX,yX]);
    
    
    for n = 1:numel(distnn)
        distMapX.data(xX(n),yX(n), iSlice) = distnn(n);
    end
end

%% Plotting
% hot or jet?
distMapX.plot('displayRange', [-1,4], 'colorMap', 'jet', ...
    sharedParameters{:})
if doUseFiguresForSaving
    axis off; title('');
else
    colorbar;
end
fhArray(end+1,1) = gcf;

if figOptions.verbosity_level > 0
    Y.name = ['distmap_' Y.name];
    Y.threshold(displayRangeFunc).plot_overlays(distMapX+1, 'overlayMode', 'map', 'overlayColorMaps', {@(x) jet(7)}, ...
        'overlayThreshold', [-0.1 4], sharedParametersOverlays{:});
    fhArray(end+1,1) = gcf;
    if doUseFiguresForSaving
        axis off; title('');
    else
        colorbar;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Histogram creation for quantification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maskVisualCortex = MrImageSpm4D(fileMaskVisual);

% only evaluate distances for edge voxels, therefore make edges an ROI

% Focus on slices that were selected / remove first/last slice (partial
% coverage)
% also: create same rois, but multiply with visual cortex mask
idxExcludedSliceRoi{1} = [];
idxExcludedSliceRoi{2} = setdiff(1:edgeX.dimInfo.nSamples('z'),figOptions.selectedSlice);
idxExcludedSliceRoi{3} = [1 edgeX.dimInfo.nSamples('z')];
roiNames = {'all slices', 'all but border slices', 'displayed slices'};
nExcludedSliceRois = numel(idxExcludedSliceRoi);

for idxRoi = 1:nExcludedSliceRois
    
    roiX{idxRoi} = edgeX.copyobj;
    for iSlice = idxExcludedSliceRoi{idxRoi}
        roiX{idxRoi}.data(:,:,iSlice) = 0;
    end
    roiX{idxRoi}.name = roiNames{idxRoi};
    roiX{idxRoi+nExcludedSliceRois} = roiX{idxRoi}.*maskVisualCortex;
    roiX{idxRoi+nExcludedSliceRois}.name = [roiX{idxRoi}.name ' - visual cortex'];
    
end

% show visual cortex mask
if figOptions.verbosity_level > 0
    tmpName = Y.name;
    Y.name = 'Visual cortex mask on edges';
    Y.threshold(displayRangeFunc).plot_overlays(roiX([end/2 end]), ...
        sharedParametersOverlays{:})
    if doUseFiguresForSaving
        axis off; title('');
    else
        colorbar;
    end
    set(gcf, 'Name', Y.name);
    Y.name = tmpName;
end

distMapX.extract_rois(roiX);
distMapX.compute_roi_stats();
distMapX.save('fileName', figOptions.output_edge_distance_map);

% plot selected ROIs
for idxRoi = 2 + [0 nExcludedSliceRois]
    distMapX.rois{idxRoi}.plot('dataGrouping', 'perSlice', 'axisType', 'relative', 'selectedSlices', selectedSlice)
    fhArray(end+1,1) = gcf;
    set(gcf, 'Name', sprintf('ROI Histograms (slices): %s',  roiX{idxRoi}.name));
    distMapX.rois{idxRoi}.plot('dataGrouping', 'perVolume', 'axisType', 'relative')
    fhArray(end+1,1) = gcf;
    set(gcf, 'Name', sprintf('ROI Histograms (volume): %s',  roiX{idxRoi}.name));
    xlabel('Distance to functional contour edge in (voxels)');
    xlim([-1 5])
end

set(fhArray, 'WindowStyle', 'docked');
% figure('Name', 'Histogram: Distribution of edge distances')
% hist(distMapX.data(distMapX.data>-1),numel(unique(distMapX.data(distMapX.data>-1))))
