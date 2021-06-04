function [fhArray, outputStats] = ...
    create_figure_activation_tissue_overlap_summary_all_subjects()
% Reports share of activated voxels in different ROIs/boundaries per
% subject
%
%   [fhArray , outputStats] = create_figure_activation_tissue_overlap_summary_all_subjects()
%
% IN
%
% OUT
%   fhArray         figure handles of created figures
%   outputStats     data used to generate histograms
% EXAMPLE
%   create_figure_activation_tissue_overlap_summary_all_subjects
%
%   See also

% Author:   Lars Kasper
% Created:  2021-04-22
% Copyright (C) 2015 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
%

options = spifi_get_analysis_options();
idxSubjectArray = 2:7;
for iSubj = 1:numel(idxSubjectArray)
    details = spifi_get_subject_details(idxSubjectArray(iSubj), options, 1, 1);
    X{iSubj}= MrImage(...
        details.representation.fig_specificity_activation_tissue_overlap.output_spm_with_tissue_rois);
    
    for iRoi = 1:numel(X{iSubj}.rois)
        nSignificantVoxelsInRoi(iRoi, iSubj) = X{iSubj}.rois{iRoi}.perVolume.nVoxels;
        %roiNameArray{iRoi} = X{iSubj}.rois{iRoi}.name;
    end
    percentageSignificantVoxelsInRoi(:,iSubj) =  ...
        nSignificantVoxelsInRoi(2:end,iSubj)/nSignificantVoxelsInRoi(1,iSubj)*100;
    
    idxUnidentified = 5;
    percentageSignificantVoxelsInIdentifiedRoi(:,iSubj) =  ...
        nSignificantVoxelsInRoi(2:end,iSubj)...
        /(nSignificantVoxelsInRoi(1,iSubj)-nSignificantVoxelsInRoi(idxUnidentified,iSubj))*100;
end

fhArray = [];

%% Stacked bar plot of share of voxels in different ROIs
fhArray(end+1,1) = tapas_physio_get_default_fig_params();
set(fhArray(end), 'Name', 'Histogram Stacked Barplot Activation Tissue Overlap per Subject');

roiNameArray = {'Activated Voxels', 'GM', 'WM', 'CSF', ...
    'Ambiguous', 'Pial Surface', ...
    'WM/GM Surface'};
idxPlottingOrder1 = [2 6 7 3 4 5];

% order - 1. because total voxel number was removed
hb = bar(idxSubjectArray, percentageSignificantVoxelsInRoi(idxPlottingOrder1-1, :)', 'stacked', 'FaceColor', 'flat');
% colormap trying to match line color in compute_spatial_specificity_activation_tissue_overlap
% not so simple, because of different color map binning there
barColorMap = jet(9);
barColorMap(1:2,:) = [];
for idxBar = 1:numel(hb)
    hb(idxBar).FaceColor = barColorMap(idxBar,:);
end


stringLegend1 = roiNameArray(idxPlottingOrder1);
% legend(stringLegend1, 'Location','NorthEastOutside');
% for paper plotting:
legend(stringLegend1, 'Location','SouthOutside', 'NumColumns', 6);
xlabel('Subjects');
ylabel('Portion of Activated Voxels');
ylim([0 100]);
ha = gca;
ha.YTick = 0:10:100;
set(ha, 'YTickLabel', ...
    cellfun(@(x) sprintf('%2.0f %%', x), num2cell(ha.YTick), 'UniformOutput', false));


%% Stacked bar plot of share of voxels in different ROIs, only considering
% clearly identified ROIs
fhArray(end+1,1) = tapas_physio_get_default_fig_params();
set(fhArray(end), 'Name', 'Histogram Stacked Barplot Activation Identified Tissues Overlap per Subject');

roiNameArray = {'Activated Voxels', 'GM', 'WM', 'CSF', ...
    'Ambiguous', 'Pial Surface', 'WM/GM Surface'};
idxPlottingOrder = [2 6 7 3 4];

% order - 1. because total voxel number was removed
hb = bar(idxSubjectArray, percentageSignificantVoxelsInIdentifiedRoi(idxPlottingOrder-1, :)', 'stacked', 'FaceColor', 'flat');
% colormap trying to match line color in compute_spatial_specificity_activation_tissue_overlap
% not so simple, because of different color map binning there
barColorMap = jet(9);
barColorMap(1:2,:) = [];
for idxBar = 1:numel(hb)
    hb(idxBar).FaceColor = barColorMap(idxBar,:);
end

stringLegend = roiNameArray(idxPlottingOrder);
legend(stringLegend, 'Location','NorthEastOutside');
xlabel('Subjects');
ylabel('Portion of Activated Voxels');
ylim([0 100]);
ha = gca;
ha.YTick = 0:10:100;
set(ha, 'YTickLabel', ...
    cellfun(@(x) sprintf('%2.0f %%', x), num2cell(ha.YTick), 'UniformOutput', false));

outputStats = struct('idxSubjectArray', idxSubjectArray, ...
    'nSignificantVoxelsInRoi', nSignificantVoxelsInRoi(idxPlottingOrder1-1,:), ...
    'stringLegend', {stringLegend1}, ...
    'percentageSignificantVoxelsInRoi', percentageSignificantVoxelsInRoi(idxPlottingOrder1-1,:), ...
    'percentageSignificantVoxelsInIdentifiedRoi', percentageSignificantVoxelsInIdentifiedRoi(idxPlottingOrder-1,:));