function [fh, Y] = create_figure_summary_results_highres_all_subjects(...
    iSess,iRecon)
% For Figure S1 (supplementary), creates plots of all transversal slices
% for all subjects of high res spiral out sequence
% - mean functional
% - overlay mean functional with edges of anatomical reference
% - SNR, SD, CoV, images
% - overlay t-maps onto mean functional image
%
%
%
%   fh = create_figure_summary_results_all_subjects(iSess, iRecon)
%
% IN
%   iSess       1 = fmri1 (high res out, default)
%               2 = fmri2 (in/out)
%   iRecon      reconstruction used (e.g. PartIn or PartOut of In/Out
%               spiral)
%               fmriOut
%               1 = magnitude
%               2 = phase
%               fmriInOut
%               1 = In          4 = In (phase)
%               2 = Out         5 = Out (phase)
%               3 = Combined    6 = Combined (phase)
% OUT
%   fh          vector of figure handles useful in the corresponding paper
%               figure
% EXAMPLE
%   create_figure_summary_results_highres_all_subjects
%
%   See also
%
% Author:   Lars Kasper
% Created:  2018-09-19
% Copyright (C) 2018 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich

% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.

options = spifi_get_analysis_options();

idxSubjectArray = 2:7;
selectedSlices = 1:36;


if nargin < 1
    iSess = 1; % high res
end

if nargin < 2
    iRecon = 1; % 1=magn, 2 = phase
end

nSubjects = numel(idxSubjectArray);
for iSubj = 1:nSubjects
    idxSubject = idxSubjectArray(iSubj);
    details = spifi_get_subject_details(idxSubject, [], ...
        iSess, iRecon);
    
    files = {
        details.preproc.func.realigned
        details.preproc.biascorrect.biasfield
        details.preproc.anat.biascorrected
        details.glm.tcon_differential1
        };
    
    nFiles = numel(files);
    
    %% load files
    for f = 1:nFiles
        Y{f} = MrImage(files{f});
    end
    
    
    % perform phase masking and bias field correction
    %  if details.isPhase % mask phase images
    %         M = MrImage(regexprep(files{f},'_phase', ''), 'selectedVolumes',1);
    %         M = M.binarize(M.prctile(55)).imerode(strel('disk',2)).imfill('holes');
    %         Y{f} = Y{f}.*M;
    %     end
    
    nFiles = numel(Y);
    
    cropX = details.representation.fig_mean.cropX;
    cropY = details.representation.fig_mean.cropY;
    if details.isPhase
        displayRange = [-1 2]; % [prctile(Y{nFiles},3), prctile(Y{nFiles},97)];
    else
        displayRange = details.representation.fig_mean.displayRange;
    end
    
    
    %% shared parameters
    nRows = 6;
    thresholdRange = [3.2 8];
    
    sharedParameters = {'plotLabels', false, 'plotTitle', false, 'rotate90', 1, ...
        'nRows',6};
    sharedParameters_overlay = {sharedParameters{:}, ...
        'overlayColorMaps', {'winter','autumn'}, ...
        'colorBar', 'off', 'overlayAlpha', 1, ...
        'overlayThreshold', thresholdRange, ...
        'overlayMode', 'map'};
    %% compute derived statistical images
    X{1} = Y{1}.select('t', 1);         % first volume realigned
    X{2} = Y{1}.mean;                   % mean volume realigned
    X{3} = X{2}.*Y{2};                  % bias-field corrected mean realigned
    X{4} = Y{1}.snr;                    % SFNR image
    X{5} = Y{1}.std;                    % SD image
    X{6} = Y{1}.compute_coeff_var();    % CoV image
    
    % mask background for CoV
    M = X{2}.copyobj;
    M = M.binarize(M.prctile(50)).imfill('holes');
    nImages = numel(X);
    
    % other backgrounds do not have to be masked
    doMaskPlot = zeros(nImages,1);
    doMaskPlot(6) = 1;
    
    displayRanges = [
        0 0.8
        0 0.8
        0 0.8
        0 30
        0 0.05
        0 0.2
        ];
    
    %% Plot images derived from realigned functional timeseries
    if iSubj == 1
        fh = zeros(nSubjects,nImages+2);
    end
    
    nameArray = {
        'firstfunc_realigned'
        'meanfunc_realigned'
        'meanfunc_nobias'
        'SFNR'
        'SD'
        'CoV'
        'meanfunc_anat_edges'
        'tMap'
        };
    for iImage = 1:nImages
        plotY = X{iImage}.copyobj;
        if doMaskPlot(iImage)
            plotY = plotY.*M;
        end
        fh(iSubj,iImage) = plotY.plot('y', [cropY(1):cropY(2)], 'x', [cropX(1):cropX(2)], ...
            'z', selectedSlices, 'displayRange', ...
            displayRanges(iImage,:), sharedParameters{:});
        axis off
        set(gcf, 'Name', sprintf('%s %d %s', details.subjectId, iImage, nameArray{iImage}));
    end
    
    %% Plot overlay edges of structural onto functional
    figOptions = options.representation.fig_congruency;
    idxImageUnderlay = 3;
    idxFileOverlay = nFiles - 1;
    underlayImage = X{idxImageUnderlay}.select.apply_threshold(displayRanges(idxImageUnderlay,:));
    overlayImage = Y{idxFileOverlay}.apply_threshold(figOptions.thresholdEdgeOutliers);
    
    % crop different resolutions to same zoom
    zI = underlayImage.crop_all({'y', cropY, 'x', cropX}, overlayImage);
    
    fh(iSubj, nImages+1) = zI{1}.plot_overlays(zI{2}, ...
        'overlayMode',  'edge', ...
        'edgeThreshold', figOptions.thresholdEdge, ...
        'colorBar', 'off', sharedParameters{:});
    axis off;
    set(gcf, 'Name', sprintf('%s %d %s', details.subjectId, nImages+1, ...
        nameArray{nImages+1}));
    
    
    %% finally, plot overlay t-maps on functional as well
    idxImageUnderlay = 3;
    idxFileOverlay = nFiles;
    
    underlayImage = X{idxImageUnderlay}.apply_threshold(...
        displayRanges(idxImageUnderlay,:));
    overlayImage = Y{idxFileOverlay};
    zI = underlayImage.crop_all({'y', cropY, 'x', cropX}, overlayImage);
    
    fh(iSubj, nImages+2) = zI{1}.plot_overlays(...
        {zI{2}, zI{2}.*-1},  sharedParameters_overlay{:}, ...
        'selectedSlices', selectedSlices);
    axis off;
    set(gcf, 'Name', sprintf('%s %s: %d %s', details.subjectId, ...
        details_recon, nImages+2, ...
        nameArray{nImages+2}));
end

end