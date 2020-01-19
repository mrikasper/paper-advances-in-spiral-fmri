function [fh, figOptions] = create_figure_spm_inout(figOptions)
% Creates activation plots of single subject for in and out part of trajectory
% overlaid on mean functional / structural image
%
%   fh = create_figure_spm_inout(figOptions)
%
% IN
%   figOptions  See also spifi_get_analysis_options:
%               options.representation.fig_...
% OUT
%   fh          vector of figure handles useful in the corresponding paper
%               figure
% EXAMPLE
%   create_figure_spm_inout
%
%   See also
%
% Author:   Lars Kasper
% Created:  2018-09-21
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
    figOptions = options.representation.fig_spm_inout;
end

fhAll = [];

views = figOptions.views;
nViews = numel(views);

for iRecon = 1:3
    details = spifi_get_subject_details(figOptions.iSubj, options, ...
        figOptions.iSess, iRecon);
    
    files = {
        details.preproc.func.biascorrected
        details.preproc.anat.biascorrected
        %details.summary.tcon % is not updated by smoothing options in, if
        % create_summary is not rerun
        details.glm.tcon
        };
    
    nFiles = numel(files);
    
    
    %% load files
    for f = 1:nFiles
        Y{f} = MrImage(files{f});
        % raised cosine k-space filter on mean spiral-in image
        if f == 1 && options.doKfilter && iRecon==1 
           Y{f} = Y{f}.kfilter('raised_cosine','2D', 'beta',0.7, ...
               'fractionFOV',0.7, 'doPlotFilter', false);
        end
    end
    
    sharedParameters = {'plotLabels', false, 'plotTitle', false, ...
        'overlayColorMaps', {'winter','autumn'}, ...
        'colorBar', 'off', 'overlayAlpha', 1, ...
        'overlayThreshold', figOptions.thresholdRange, ...
        'overlayMode', 'map', 'nRows',1};
    
    sharedParameters_tra = {sharedParameters{:}, ...
        'rotate90', 1};
    
    sharedParameters_sag = {sharedParameters{:}, ...
        'rotate90', 2, 'sliceDimension', 1};
    
    sharedParameters_cor = {sharedParameters{:}, ...
        'rotate90', 1, 'sliceDimension', 2};
    
    
    if isnumeric(figOptions.selectedSlice) % slice indices given explicitly
        selectedSlices = figOptions.selectedSlice;
        sliceString = sprintf('_%d', selectedSlices);
    else
        sliceString = ['_' figOptions.selectedSlice];
        switch lower(figOptions.selectedSlice)
            case 'max'
                [tMax, posMax] = Y{nFiles}.max;
                selectedSlices = posMax;
            case 'min'
                [tMin, posMin] = Y{nFiles}.min;
                selectedSlices = posMin;
            case 'maxextent'
                [selectedSlices, centre] = spifi_find_most_activated_slices(...
                    files{nFiles}, figOptions.thresholdRange(1));
            case 'maxextent_l'
                statMap = MrImage(files{nFiles});
                % use only one hemisphere to determine max
                iSliceHemisphereBorder = statMap.dimInfo.sample2index(0,1);
                statMapL = statMap.select('x', 1:iSliceHemisphereBorder);
                [selectedSlices, centre] = spifi_find_most_activated_slices(...
                    statMapL, figOptions.thresholdRange(1));
            case 'maxextent_r'
                statMap = MrImage(files{nFiles});
                
                % use only one hemisphere to determine max
                iSliceHemisphereBorder = statMap.dimInfo.sample2index(0,1);
                statMapR = statMap.select('x', iSliceHemisphereBorder+1:statMap.dimInfo.nSamples(1));
                [selectedSlices, centre] = spifi_find_most_activated_slices(...
                    statMapR, figOptions.thresholdRange(1));
                % add offset from selection of second half of voxels again
                selectedSlices(1) = selectedSlices(1) + iSliceHemisphereBorder;
        end
    end
    
    % Different display window for in and out part...
    % Extend, if you use more underlays
    switch iRecon
        case 1
            displayRanges = [
                figOptions.displayRangeIn
                figOptions.displayRangeStructural
                ];
        case 2
            displayRanges = [
                figOptions.displayRangeOut
                figOptions.displayRangeStructural
                ];
        case 3
            displayRanges = [
                figOptions.displayRangeCombined
                figOptions.displayRangeStructural
                ];
    end
    
    fh = zeros(nViews*(nFiles-1),1);
    crop_tra = details.representation.fig_spm_inout.crop_tra;
    %     crop_tra(1:2) = round(crop_tra(1:2)*190/240);
    %     crop_tra(3:4) = round(crop_tra(3:4)*230/292);
    crop_tra(1:2) = round(crop_tra(1:2)*160/240);
    crop_tra(3:4) = round(crop_tra(3:4)*192/292);
    for f = 1:nFiles-1
        
        resizedUnderlay = Y{f}.resize(Y{nFiles}.dimInfo);
        
        %% loop over all views, creating separate figures
        for iView = 1:nViews
            %% transversal view, whole slice
            switch lower(views{iView})
                case 'tra'
                    % zoom all images to transversal crop
                    zI = resizedUnderlay.crop_all({'x',crop_tra(1:2), 'y', crop_tra(3:4)}, ...
                        Y{nFiles});
                    
                    fh(iView+(nFiles-1)*(f-1)) = ...
                        zI{1}.apply_threshold(displayRanges(f,:)).plot_overlays(...
                        {zI{2}, zI{2}.*-1}, ...
                        sharedParameters_tra{:}, ...
                        'selectedSlices', selectedSlices);
                    
                    %% zoomed visual cortex transversal
                    % zoom all images to crop
                case 'tra_zoom'
                    zI = resizedUnderlay.crop_all({'x',crop_tra(1:2)+[5 -5], 'y', ...
                        crop_tra(3)+[0 70]}, ...
                        Y{nFiles});
                    
                    fh(iView+(nFiles-1)*(f-1)) = ...
                        zI{1}.apply_threshold(displayRanges(f,:)).plot_overlays(...
                        {zI{2}, zI{2}.*-1}, ...
                        sharedParameters_tra{:}, ...
                        'selectedSlices', selectedSlices(1)); % HACK: only 1 slice shown
                case 'sag'
                    zI = resizedUnderlay.crop_all({
                        'y', ...
                        crop_tra(3)+[0 70]}, ...
                        Y{nFiles});
                    
                    fh(iView+(nFiles-1)*(f-1)) = ...
                        zI{1}.apply_threshold(displayRanges(f,:)).plot_overlays(...
                        {zI{2}, zI{2}.*-1}, ...
                        sharedParameters_sag{:}, ...
                        'selectedSlices', selectedSlices);
                case 'cor'
                    zI = resizedUnderlay.crop_all({
                        'x', ...
                        crop_tra(1:2)+[25 -25]}, ...
                        Y{nFiles});
                    
                    fh(iView+(nFiles-1)*(f-1)) = ...
                        zI{1}.apply_threshold(displayRanges(f,:)).plot_overlays(...
                        {zI{2}, zI{2}.*-1}, ...
                        sharedParameters_cor{:}, ...
                        'selectedSlices', selectedSlices);
            end
            
            % final cosmetics and figure labelling
            axis off;
            set(gcf, 'Name', spm_file(...
                details.representation.fig_spm_group.fileSaveArray{f}, ...
                'suffix', sprintf('_%s%s', views{iView}, sliceString)));
        end
    end
    fhAll = [fhAll; fh]; % concat different in/out recon figs
end

fh = fhAll;