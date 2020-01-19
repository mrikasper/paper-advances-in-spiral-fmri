function [summaryTable, tabOptions] = create_table_tsnr_spm_group(tabOptions)
% Creates activation plots of single subject for in and out part of trajectory
% overlaid on mean functional / structural image
%
%   summaryTable = create_table_tsnr_spm_group(figOptions)
%
% IN
%   tabOptions  as in spifi_get_analysis_options:
%               options.representation.tab_...
% OUT
%   summaryTable          vector of table handles useful in the corresponding paper
%               table
% EXAMPLE
%   create_table_tsnr_spm_group
%
%   See also spifi_get_analysis_options
%
% Author:   Lars Kasper
% Created:  2018-10-21
% Copyright (C) 2018 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich

% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.

selectedSlices  = tabOptions.selectedSlices;
doPlotMasks     = tabOptions.doPlotMasks;
idxSubjects     = tabOptions.iSubjects;
idxPairs        = tabOptions.idxPairs;
threshold       = tabOptions.clusterFormingThreshold;

% pairs of recon
sessReconPairs = [
    1 1
    2 1
    2 2
    2 3
    ];

% names of sess/recon pairs
namePairs = {
    'high-res spiral-out (0.8 mm)'
    'in-part spiral in/out (1.5mm)'
    'out-part spiral in/out (1.5mm)'
    'signal-weighted combination spiral in/out (1.5mm)'
    };

nPairs = numel(idxPairs);
nSubjects = numel(idxSubjects);

for iPair = 1:nPairs
    idxPair = idxPairs(iPair);
    namePair = namePairs{idxPair};
    fprintf('\n\n\tEvaluating sequence %s, %d/%d\n\n', namePair, ...
        iPair, nPairs);
    for iSubj = 1:nSubjects
        iRow = (iPair-1)*nSubjects + iSubj;
        iSess = sessReconPairs(idxPair,1);
        iRecon = sessReconPairs(idxPair,2);
        idxSubject = idxSubjects(iSubj);
        
        details = spifi_get_subject_details(idxSubject, ...
            [], iSess, iRecon);
        
        fprintf('\n\t\tEvaluating subject %s, %d/%d\n\n', details.subjectId, ...
            iSubj, nSubjects);
        
        % tSNR, differential t-contrast, GM mask
        files = {
            details.summary.snr.realigned
            details.summary.tcon
            details.summary.mean.realigned
            details.preproc.anat.tpm_GM
            };
        
        nFiles = numel(files);
        
        %% load files
        for f = 1:nFiles
            Y{f} = MrImage(files{f});
        end
        
        %% Create GM mask and check displayed on the mean
        
        % gray matter mask
        %       M = Y{end}.Y{end}.copyobj.resize(Y{1}.dimInfo).binarize(0.8).imerode(strel('disk',1));
        % Note: First resize, then binarize, to avoid interpolation errors
        M = Y{end}.copyobj.resize(Y{1}.dimInfo).binarize(0.8);
        
        if doPlotMasks
            % plot with mean underlay
            Y{end-1}.apply_threshold([0 0.8]).plot_overlays(M, ...
                'overlayMode', 'edge', 'rotate90',1, 'selectedSlices', selectedSlices)
            set(gcf, 'Name', [details.subjectId ' - GM mask on Mean realigned fMRI'])
            
            % plot edges with tSNR underlay
            Y{1}.apply_threshold([0 30]).plot_overlays(M, ...
                'overlayMode', 'edge', 'rotate90',1, 'selectedSlices', selectedSlices)
            set(gcf, 'Name', [details.subjectId ' - GM mask on SFNR of realigned fMRI'])
            
            % plot mask-multiplied tSNR image
            plot(Y{1}.*M, 'z', selectedSlices, 'rotate90',1, ...
                'displayRange', [0 30], 'colorBar', 'on')
            set(gcf, 'Name', [details.subjectId ' - GM mask times SFNR of realigned fMRI'])
            
            % plot t-contrast with GM mask
            Y{2}.plot_overlays(M, ...
                'overlayMode', 'mask', ...
                'overlayAlpha', 0.5, ...
                'rotate90',1, ...
                'selectedSlices', selectedSlices)
            set(gcf, 'Name', [details.subjectId ' - GM mask on differential t-contrast'])
        end
        
        %% Extract SFNR data from GM ROI
        % remove incomplete cut-off slices (realignment!) from stats,
        % because different cut in functional and structural possible
        SNR = Y{1}.select('z', selectedSlices);
        M2 = M.select('z', selectedSlices);
        SNR.extract_rois(M2);
        SNR.compute_roi_stats();
        
        if doPlotMasks
            SNR.rois{1}.plot;
        end
        
        % gather for later table
        SFNR_mean(iRow,1) = SNR.rois{1}.perVolume.mean;
        SFNR_SD(iRow,1) = SNR.rois{1}.perVolume.sd;
        RowNames{iRow,1} = sprintf('%s - %s', namePair, details.subjectId);
        
        
        %% Extract data from T-map
        idxContrast = 1;
        nMinVoxels = get_fwe_clustersize(details.glm.spm_mat, idxContrast);
        % maximum peak over both contrasts
        TMAP = Y{2};
        SPM_T_max(iRow,1) = TMAP.abs.max;
        
        
        % note: for t-maps, this takes positive and negative peaks into account
        Msignificant = TMAP.copyobj.abs.binarize(threshold);
        switch tabOptions.FWEcorrection
            case 'none'
            case 'cluster'
                % remove small clusters that don't survive multiple
                % comparison correction
                Msignificant = Msignificant.remove_clusters([1 nMinVoxels-1]);
            case 'peak'
                % TODO
        end
        % remove small clusters (multiple comparison correction)
        
        TMAP.extract_rois(Msignificant);
        SPM_T_nVoxels(iRow,1) =  TMAP.rois{1}.perVolume.nVoxels;
        SPM_T_volume(iRow,1) =  TMAP.rois{1}.perVolume.nVoxels.*...
            prod(TMAP.dimInfo.resolutions); % normalized for differing resolutions
    end
end

summaryTable = table(SFNR_mean, SFNR_SD, SPM_T_max, SPM_T_nVoxels, ...
    SPM_T_volume, ...
    'RowNames', RowNames);