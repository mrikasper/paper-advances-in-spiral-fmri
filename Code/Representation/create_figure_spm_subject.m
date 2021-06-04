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

if nargin < 1
    options = spifi_get_analysis_options();
    figOptions = options.representation.fig_spm_subject;
end

details = spifi_get_subject_details(figOptions.iSubj, ...
    [], figOptions.iSess, figOptions.iRecon);

files = {
    details.preproc.func.biascorrected
    details.preproc.anat.biascorrected
    ''
    ''
    };

% funcColorMaps = {
%     @(x) split_brewermap(x, 'RdYlBu', 'pos')
%     @(x) split_brewermap(x, 'RdYlBu', 'neg')
%     };
% funcColorMaps = {
%     @(x) split_brewermap(x, 'parula', 'neg')
%     @(x) split_brewermap(x, 'parula', 'pos')
%     };
%

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
        %         funcColorMaps = {
        %             'winter'
        %             'autumn'
        %             };
        %         overlayAlpha = 1;
        funcColorMaps = {
            @(nColors) get_colormap_gradient(9/10*[0 1 1], [0 1 1], nColors)
            @(nColors) get_colormap_gradient(9/10*[1 1 0], [1 1 0], nColors)
            };
        overlayAlpha = 0.6;
end

nFiles = numel(files);

sharedParameters = {'plotLabels', false, 'plotTitle', false, ...
    'colorBar', 'off', 'overlayAlpha', overlayAlpha, ...
    'overlayThreshold', figOptions.thresholdRange, ...
    'overlayMode', 'map', ...
    'overlayColorMaps', funcColorMaps};

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
    if f==2, Y{f}.apply_threshold(0.1, 'exclude'); end
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



views = figOptions.views;
nViews = numel(views);

fh = zeros((nFiles-2)*nViews,1); % last 2 images are contrasts

sliceString = sprintf('_%s', figOptions.contrastType);

for f = 1:nFiles-2
    
    % resize underlay to actual resolution of activation map
    % otherwise, voxel shifts may occur
    resizedUnderlay = Y{f}.resize(Y{nFiles}.dimInfo);
    
    centerX = round(resizedUnderlay.dimInfo.x.nSamples/2);
    
    for iView = 1:nViews
        switch lower(views{iView})
            case 'tra'
                %% transverse plot
                fh((f-1)*nViews+iView) = resizedUnderlay.apply_threshold(displayRanges(f,:)).plot_overlays(...
                    {Y{nFiles-1}, Y{nFiles}}, ...
                    'selectedSlices', figOptions.selectedSlice, ...
                    sharedParameters_tra{:}, ...
                    'nRows', 1);
                
            case 'tra_zoom'
                %% zoomed transversal plot
                zI = resizedUnderlay.crop_all({'y', [10 70]}, Y(nFiles+[-1 0]));
                
                fh((f-1)*nViews+iView) = zI{1}.apply_threshold(displayRanges(f,:)).plot_overlays(...
                    zI(end-1:end), ...
                    'selectedSlices', figOptions.selectedSlice, ...
                    sharedParameters_tra{:}, ...
                    'nRows', 4);
            case 'sag'
                %%  sagittal plot
                zI = resizedUnderlay.crop_all({'x', [100 100]}, Y(nFiles+[-1 0]));
                fh((f-1)*nViews+iView) = zI{1}.apply_threshold(displayRanges(f,:)).plot_overlays(...
                    zI(end-1:end), ...
                    sharedParameters_sag{:}, ...
                    'nRows', NaN);
            case 'sag_zoom'
                %% zoom sagittal
                zI = resizedUnderlay.crop_all({'y', [10 70]}, Y(nFiles + [-1 0]));
                fh((f-1)*nViews+iView) = zI{1}.apply_threshold(displayRanges(f,:)).plot_overlays(...
                    zI(end-1:end), ...
                    'selectedSlices', 100, ...
                    sharedParameters_sag{:}, ...
                    'nRows',NaN);
            case 'cor_zoom'
                %% zoomed coronal view
                % zI = Y{f}.crop_all({'x', centerX+[-60 59], 'y', [12 43]}, Y(nFiles + [-1 0]));
                zI = resizedUnderlay.crop_all({'x', centerX+[-60 59], 'y', [37 37]}, Y(nFiles + [-1 0]));
                fh((f-1)*nViews+iView) = zI{1}.apply_threshold(displayRanges(f,:)).plot_overlays(...
                    zI(end-1:end), ...
                    sharedParameters_cor{:}, ...
                    'nRows', NaN);
        end
        % final cosmetics and figure labelling
        axis off;
        set(gcf, 'Name', spm_file(...
            details.representation.fig_spm_subject.fileSaveArray{f}, ...
            'suffix', sprintf('_%s%s', views{iView},sliceString)));
    end
end