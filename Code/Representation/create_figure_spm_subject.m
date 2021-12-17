function [fh, figOptions] = create_figure_spm_subject(figOptions)
% Creates activation plots of single subject, overlaid on mean functional /
% structural
%
%   fh = create_figure_spm_subject(figOptions)
%
% IN
%   figOptions  See also spifi_get_analysis_options:
%               options.representation.fig_...
% OUT
%   fh          vector of figure handles useful in the corresponding paper
%               figure
% EXAMPLE
%   create_figure_spm_subject
%
%   See also
%
% Author:   Lars Kasper
% Created:  2018-09-20
% Copyright (C) 2018 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich

% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.

options = spifi_get_analysis_options();
if nargin < 1
    figOptions = options.representation.fig_spm_subject;
end

options.glm.doUseSmoothed = figOptions.doUseSmoothedMaps;
details = spifi_get_subject_details(figOptions.iSubj, ...
    options, figOptions.iSess, figOptions.iRecon);

files = {
    details.preproc.func.biascorrected
    details.preproc.anat.biascorrected
    ''
    ''
    };


switch figOptions.contrastType
    case 'differential'
        files{3} = details.glm.tcon_differential1;
        files{4} = details.glm.tcon_differential2;
        funcColorMaps = {
            'winter'
            'autumn'
            };
        overlayAlpha = 1;
    case 'condition'
        files{3} = details.glm.tcon_condition1;
        files{4} = details.glm.tcon_condition2;
        funcColorMaps = {
            @(nColors) get_colormap_gradient(9/10*[0 1 1], [0 1 1], nColors)
            @(nColors) get_colormap_gradient(9/10*[1 1 0], [1 1 0], nColors)
            };
        overlayAlpha = 0.6;
end

nFiles = numel(files);


if figOptions.doPlotMask
    
    % different mask colors for smoothed (red/blue) and unsmoothed
    % (yellow/cyan)
    if figOptions.doUseSmoothedMaps
        funcColorMaps = {
            @(nColors) get_colormap_gradient(9/10*[0 0 1], [0 0 1], nColors)
            @(nColors) get_colormap_gradient(9/10*[1 0 0], [1 0 0], nColors)
            };
    else
        funcColorMaps = {
            @(nColors) get_colormap_gradient(9/10*[0 1 0.5], [0 1 0.5], nColors)
            @(nColors) get_colormap_gradient(9/10*[1 1 0], [1 1 0], nColors)
            };
    end
    overlayAlpha = 0.6;
    % dummy image with all zeros but same size as original anat image
    sharedParameters = {'plotLabels', false, 'plotTitle', false, ...
        'colorBar', 'off', 'overlayAlpha', overlayAlpha, ...
        'overlayMode', 'map', ...
        'overlayColorMaps', funcColorMaps
        };
else
    sharedParameters = {'plotLabels', false, 'plotTitle', false, ...
        'colorBar', 'off', 'overlayAlpha', overlayAlpha, ...
        'overlayThreshold', figOptions.thresholdRange, ...
        'overlayMode', 'map', ...
        'overlayColorMaps', funcColorMaps};
end

sharedParameters_tra = {sharedParameters{:}, ...
    'rotate90', 1};

sharedParameters_sag = {sharedParameters{:}, ...
    'rotate90', 2, 'sliceDimension', 1};

sharedParameters_cor = {sharedParameters{:}, ...
    'rotate90', 1, 'sliceDimension', 2};

refCropX = details.representation.fig_spm_subject.cropX;
refCropY = details.representation.fig_spm_subject.cropY;

%% load files
for f = 1:nFiles
    Y{f} = MrImage(files{f});
    % remove negative values (grey background) from anatomy
    if f==2, Y{f}.threshold(0.1, 'exclude'); end
    % adapt crop to number of samples in array
    if f == 1, Yref = Y{1}.copyobj(); end % to not crop twice
    [cropX, cropY] = adjust_crop(refCropX, refCropY, ...
        Yref, Y{f});
    Y{f} = Y{f}.select('x', cropX(1):cropX(2), 'y', cropY(1):cropY(2));
end

% Extend, if you use more underlays
displayRanges = [
    details.representation.fig_spm_subject.displayRangeFunctional
    figOptions.displayRangeStructural
    ];



% perform multiple comparison correction on activation maps
% remove small clusters (multiple comparison correction)
absY = abs(Y{nFiles});
idxContrast = 1;
maskActivation = absY.binarize(figOptions.thresholdRange(1));
switch figOptions.FWEcorrection
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


for f = 3:nFiles
    Y{f} = Y{f}.*maskActivation
end

views = figOptions.views;
nViews = numel(views);

sliceString = sprintf('_%s_doUseSmoothedMaps%d_thresholdMin%3.1f_FWEcorrection-%s', figOptions.contrastType, ...
    figOptions.doUseSmoothedMaps, figOptions.thresholdRange(1), ...
    figOptions.FWEcorrection);

if figOptions.doPlotMask
    sliceString = [sliceString '_mask'];
    
    % binarize into masks for edges
    for f = 3:nFiles
        Y{f} = Y{f}.binarize(figOptions.thresholdRange(1));
    end
    
end

switch figOptions.underlays
    case 'functional'
        idxUnderlayArray = 1;
    case 'structural'
        idxUnderlayArray = 2;
    case 'both'
        idxUnderlayArray = 1:nFiles-2; % last 2 images are contrasts
    case 'none'
        % e.g., for edges;
        idxUnderlayArray = 1;
end

nUnderlays = numel(idxUnderlayArray);
fh = zeros(nUnderlays*nViews,1); 
for f = 1:nUnderlays
    
    idxUnderlay = idxUnderlayArray(f);
    % resize underlay to actual resolution of activation map
    % otherwise, voxel shifts may occur
    resizedUnderlay = Y{idxUnderlay}.resize(Y{nFiles}.dimInfo);
    
    % black background
    if contains(figOptions.underlays, 'none')
        resizedUnderlay.data(:) = 0;
        sliceString = [sliceString '_noUnderlay'];
    end
    
    centerX = round(resizedUnderlay.dimInfo.x.nSamples/2);
    
    for iView = 1:nViews
        nRows = figOptions.nRows(iView);
        switch lower(views{iView})
            case 'tra'
                %% transverse plot
                fh((f-1)*nViews+iView) = resizedUnderlay.threshold(...
                    displayRanges(idxUnderlay,:)).plot_overlays(...
                    {Y{nFiles-1}, Y{nFiles}}, ...
                    'selectedSlices', figOptions.selectedSlice, ...
                    sharedParameters_tra{:}, ...
                    'nRows', nRows);
                
            case 'tra_zoom'
                %% zoomed transversal plot
                zI = resizedUnderlay.crop_all({'y', [10 70]}, Y(nFiles+[-1 0]));
                
                fh((f-1)*nViews+iView) = zI{1}.threshold(...
                    displayRanges(idxUnderlay,:)).plot_overlays(...
                    zI(end-1:end), ...
                    'selectedSlices', figOptions.selectedSlice, ...
                    sharedParameters_tra{:}, ...
                    'nRows', nRows);
            case 'sag'
                %%  sagittal plot
                zI = resizedUnderlay.crop_all({'x', [100 100]}, Y(nFiles+[-1 0]));
                fh((f-1)*nViews+iView) = zI{1}.threshold(...
                    displayRanges(idxUnderlay,:)).plot_overlays(...
                    zI(end-1:end), ...
                    sharedParameters_sag{:}, ...
                    'nRows', nRows);
            case 'sag_zoom'
                %% zoom sagittal
                zI = resizedUnderlay.crop_all({'y', [10 70]}, Y(nFiles + [-1 0]));
                fh((f-1)*nViews+iView) = zI{1}.threshold(...
                    displayRanges(idxUnderlay,:)).plot_overlays(...
                    zI(end-1:end), ...
                    'selectedSlices', 100, ...
                    sharedParameters_sag{:}, ...
                    'nRows', nRows);
            case 'cor_zoom'
                %% zoomed coronal view
                % zI = Y{f}.crop_all({'x', centerX+[-60 59], 'y', [12 43]}, Y(nFiles + [-1 0]));
                zI = resizedUnderlay.crop_all({'x', centerX+[-60 59], 'y', [37 37]}, Y(nFiles + [-1 0]));
                fh((f-1)*nViews+iView) = zI{1}.threshold(...
                    displayRanges(idxUnderlay,:)).plot_overlays(...
                    zI(end-1:end), ...
                    sharedParameters_cor{:}, ...
                    'nRows', nRows);
        end
        % final cosmetics and figure labelling
        axis off;
        set(gcf, 'Name', spm_file(...
            details.representation.fig_spm_subject.fileSaveArray{idxUnderlay}, ...
            'suffix', sprintf('_%s%s', views{iView},sliceString)));
    end
end