function [fh, roiArray, outputStats] = create_figure_contour_congruency_summary_all_subjects()
% Line plot summary (over Slices) of Contour Distances func/struct images
% over subject
%
%    [fh, roiArray, outputStats] = create_figure_contour_congruency_summary_all_subjects()
%
% IN
%
% OUT
%
% EXAMPLE
%   create_figure_contour_congruency_summary_all_subjects
%
%   See also

% Author:   Lars Kasper
% Created:  2021-04-21
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

verbosityLevel = 1; %1,2,3 for ever more plots

idxSlices = 2:34;

options = spifi_get_analysis_options();
idxRoi = 3;
idxSubjectArray = 2:7;
for iSubj = 1:numel(idxSubjectArray)
    idxSubj = idxSubjectArray(iSubj);
    details = spifi_get_subject_details(idxSubj, options, 1, 1);
    X{iSubj} = MrImage(details.representation.fig_specificity_contours.output_edge_distance_map);
    roiArray{iSubj} = X{iSubj}.rois{idxRoi};
    contourDistanceMean(:,iSubj) =  roiArray{iSubj}.perSlice.mean;
    contourDistanceMedian(:,iSubj) = roiArray{iSubj}.perSlice.median;
    contourDistanceSd(:,iSubj) =  roiArray{iSubj}.perSlice.sd;
    stringLegend{iSubj} = sprintf('Subject %d',idxSubj);
    
    
    for idxSlice = idxSlices
        contourDistanceZeroPortion(idxSlice,iSubj) = sum(roiArray{iSubj}.data{idxSlice}==0)/roiArray{iSubj}.perSlice.nVoxels(idxSlice);
        contourDistanceOnePortion(idxSlice,iSubj) = sum(roiArray{iSubj}.data{idxSlice}==1)/roiArray{iSubj}.perSlice.nVoxels(idxSlice);
    end
    
    % collect all distances
    allData = cell2mat(roiArray{iSubj}.data);
    % distanceBinValues = unique(allData);
    % distanceBinValues = [0:0.5:3];
    distanceBinEdges = [0 0.01 1.01 1.51 2.01 2.51 3.01 max(allData)];
    stringLegend2 = {'<= 0.0 voxels', '<= 1.0 voxels', '<= 1.5 voxels', ...
        '<= 2.0 voxels', '<= 2.5 voxels', '<= 3.0 voxels', '>   3.0 voxels'};
    
    figure;
    hh = histogram(allData, distanceBinEdges, 'Normalization', 'probability');
    histDistanceVolume(:,iSubj) = hh.Values*100;
    close(gcf);
end

fh = [];
fh(end+1,1) = tapas_physio_get_default_fig_params();
% set(fh(end), 'DefaultAxesColorOrder',    ...     
%     [0    0.4470    0.7410
%     0.8500    0.3250    0.0980
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840
%     ]);
set(fh(end), 'Name', 'Mean Contour Distance over Slices per subject');
plot(idxSlices, contourDistanceMean(idxSlices,:));
legend(stringLegend);
xlabel('Slices');
ylabel('Mean Contour Distance (voxels)');

if verbosityLevel >=3
    fh(end+1,1) = tapas_physio_get_default_fig_params();
    set(fh(end), 'Name', 'Median Contour Distance over Slices per subject');
    idxSlices = 2:34;
    plot(idxSlices, contourDistanceMedian(idxSlices,:));
    legend(stringLegend);
    xlabel('Slices');
    ylabel('Median Contour Distance (voxels)');
    
    fh(end+1,1) = tapas_physio_get_default_fig_params();
    set(fh(end), 'Name', 'SD of Contour Distance over Slices per subject');
    plot(idxSlices, contourDistanceSd(idxSlices,:));
    legend(stringLegend);
    xlabel('Slices');
    ylabel('SD of Contour Distance (voxels)');
end

%% plot of percentage of voxels with 1 or less distance
if verbosityLevel >=3
    fh(end+1,1) = tapas_physio_get_default_fig_params();
    set(fh(end), 'Name', 'Portion of Contour Distance 0 over Slices per subject');
    idxSlices = 2:34;
    plot(idxSlices, 100*contourDistanceZeroPortion(idxSlices,:));
    legend(stringLegend);
    xlabel('Slices');
    ylabel('Portion of Contour Distance equal to 0 (percent)');
end

if verbosityLevel >=2
    fh(end+1,1) = tapas_physio_get_default_fig_params();
    set(fh(end), 'Name', 'Portion of Contour Distance equal to 0 or 1 over Slices per subject');
    idxSlices = 2:34;
    plot(idxSlices, 100*(contourDistanceZeroPortion(idxSlices,:)+contourDistanceOnePortion(idxSlices,:)));
    legend(stringLegend);
    xlabel('Slices');
    ylabel('Portion of Contour Distance 0 or 1(percent)');
    
    fh(end+1,1) = tapas_physio_get_default_fig_params();
    set(fh(end), 'Name', 'Portion of Contour Distance 0 or 1 over Slices per subject');
    idxSlices = 2:34;
end

%% Stacked bar plot of volume distribution of distances per subject
fh(end+1,1) = tapas_physio_get_default_fig_params();
set(fh(end), 'Name', 'Histogram Stacked Barplot Contour Distances per Subject');
hb = bar(idxSubjectArray, histDistanceVolume, 'stacked', 'FaceColor', 'flat');
% colormap trying to match line color in compute_spatial_specificity_activation_tissue_overlap
% not so simple, because of different color map binning there
barColorMap = jet(9); 
barColorMap(1:2,:) = [];
for idxBar = 1:numel(hb)
    hb(idxBar).FaceColor = barColorMap(idxBar,:);
end
legend(stringLegend2, 'Location','NorthEastOutside');
xlabel('Subjects');
ylabel('Contour Distance Distribution');
ylim([0 100]);
ha = gca;
ha.YTick = 0:10:100;
set(ha, 'YTickLabel', ...
    cellfun(@(x) sprintf('%2.0f %%', x), num2cell(ha.YTick), 'UniformOutput', false));

%% report some outputs in command line
contourDistanceMean_all_subjects = cell2mat(cellfun(@(x) x.perVolume.mean, roiArray, 'UniformOutput', false))
contourDistanceMean_all_subjects_mean = mean(contourDistanceMean_all_subjects)
contourDistanceMean_all_subjects_mm = contourDistanceMean_all_subjects*0.8

histDistanceVolume

outputStats = struct('contourDistanceMean_all_subjects', contourDistanceMean_all_subjects, ...
    'contourDistanceMean', contourDistanceMean, ...
    'histDistanceVolume', histDistanceVolume);