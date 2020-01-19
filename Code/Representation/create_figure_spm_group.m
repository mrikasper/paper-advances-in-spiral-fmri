function [fh, figOptions] = create_figure_spm_group(figOptions)
% Creates activation plots all subjects at max activation peak,
% overlaid on mean functional / structural
%
%   fh = create_figure_spm_group(figOptions)
%
% IN
%   figOptions  See also spifi_get_analysis_options:
%               options.representation.fig_...
% OUT
%   fh          vector of figure handles useful in the corresponding paper
%               figure
% EXAMPLE
%   create_figure_spm_group
%
%   See also
%
% Author:   Lars Kasper
% Created:  2019-09-01
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
    figOptions = options.representation.fig_spm_group;
end

idxSubjects = figOptions.iSubjects;
nSubjects = numel(idxSubjects);
for iSubj = 1:nSubjects
    idxSubject = idxSubjects(iSubj);
    
    details = spifi_get_subject_details(idxSubject, ...
        [], figOptions.iSess, figOptions.iRecon);
    
    files = {
        %details.summary.mean.realigned
        %details.summary.warpedanat
        details.preproc.func.biascorrected
        details.preproc.anat.biascorrected
        %details.summary.tcon % is not updated by smoothing options in, if
        % create_summary is not rerun
        details.glm.tcon
        };
    
    nFiles = numel(files);
    
    % prepare figure handle vector
    if iSubj == 1
        fh = zeros((nFiles-1)*numel(figOptions.iSubjects),1);
    end
    
    %% load files
    for f = 1:nFiles
        Y{f} = MrImage(files{f});
    end
    
    % Extend, if you use more underlays
    displayRanges = [
        figOptions.displayRangeFunctional
        figOptions.displayRangeStructural
        ];
    
    
    sharedParameters = {'plotLabels', true, 'plotTitle', false, ...
        'overlayColorMaps', {'winter','autumn'}, ...
        'colorBar', 'off', 'overlayAlpha', 1, ...
        'overlayThreshold', figOptions.thresholdRange, ...
        'overlayMode', 'map'};
    
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
    
    crop_tra = figOptions.crop_tra{idxSubject};
    views = figOptions.views;
    nViews = numel(views);
    for f = 1:nFiles-1
        
        % resize underlay to actual resolution of activation map
        resizedUnderlay = Y{f}.resize(Y{nFiles}.dimInfo);
        
        %% loop over all views, creating separate figures
        for iView = 1:nViews
            %% transversal view, whole slice
            switch lower(views{iView})
                case 'tra'
                    % zoom all images to transversal crop
                    zI = resizedUnderlay.crop_all({'x',crop_tra(1:2), 'y', crop_tra(3:4)}, ...
                        Y{nFiles});
                    
                    fh(f+(nFiles-1)*(iSubj-1)) = ...
                        zI{1}.apply_threshold(displayRanges(f,:)).plot_overlays(...
                        {zI{2}, zI{2}.*-1}, ...
                        sharedParameters_tra{:}, ...
                        'selectedSlices', selectedSlices(3));
                    
                    %% zoomed visual cortex transversal
                    % zoom all images to crop
                case 'tra_zoom'
                    zI = resizedUnderlay.crop_all({'x',crop_tra(1:2)+[5 -5], 'y', ...
                        crop_tra(3)+[0 70]}, ...
                        Y{nFiles});
                    
                    fh(f+(nFiles-1)*(iSubj-1)) = ...
                        zI{1}.apply_threshold(displayRanges(f,:)).plot_overlays(...
                        {zI{2}, zI{2}.*-1}, ...
                        sharedParameters_tra{:}, ...
                        'selectedSlices', selectedSlices(3));
                case 'sag'
                    zI = resizedUnderlay.crop_all({
                        'y', ...
                        crop_tra(3)+[0 70]}, ...
                        Y{nFiles});
                    
                    fh(f+(nFiles-1)*(iSubj-1)) = ...
                        zI{1}.apply_threshold(displayRanges(f,:)).plot_overlays(...
                        {zI{2}, zI{2}.*-1}, ...
                        sharedParameters_sag{:}, ...
                        'selectedSlices', selectedSlices(1));
                case 'cor'
                    zI = resizedUnderlay.crop_all({
                        'x', ...
                        crop_tra(1:2)+[25 -25]}, ...
                        Y{nFiles});
                    
                    fh(f+(nFiles-1)*(iSubj-1)) = ...
                        zI{1}.apply_threshold(displayRanges(f,:)).plot_overlays(...
                        {zI{2}, zI{2}.*-1}, ...
                        sharedParameters_cor{:}, ...
                        'selectedSlices', selectedSlices(2));
            end
            
            % final cosmetics and figure labeling
            axis off;
            set(gcf, 'Name', spm_file(...
                details.representation.fig_spm_group.fileSaveArray{f}, ...
                'suffix', sprintf('_%s%s', views{iView}, sliceString)));
        end
        
    end
end