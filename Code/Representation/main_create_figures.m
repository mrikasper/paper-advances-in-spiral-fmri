% Script main_create_figures
% Generates all figures of Spiral Functional Imaging Paper (SPIFI)
%
%  main_create_figures
%
%
%   See also
%
% Author:   Lars Kasper
% Created:  2018-09-19
% Copyright (C) 2018 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spifi_setup_paths();
options = spifi_get_analysis_options();
paths = options.paths();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop to create all paper figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0 = Table 1;
% Supplementary Figures: 
%   11 = Figure S1 (all reported means and maps for high-res spiral Out)
%   12 = Figure S2; (all reported means and maps for spiral in/out/combined)
%   14 =  Figure S4 (supplementary video smoothed/unsmoothed/cropped recon)
idxFiguresToPlot    = 17; % [0:12];
% false or true; true saves them to Results\Summary\PipelineXX\Figures\FigureYY
doSaveFigures       = false;

% configurations of subjects, session (high-res out or in/out) and recon
% (magnitude/phase and in/out/combined (for in/out)) taken for display
% idxSubj, idxSess, idxRecon
idxSubjSessReconArray = [
    2, 1, 1
    2, 2, 1
    2, 2, 2
    ];

nConfigs = size(idxSubjSessReconArray,1);
% change to subset of [2:7] to create same plots for other subjects
% include 12 for cropped-to-1mm reconstruction of SPIFI_0002
% NOTE: for Fig. 21 and 22, all subjects have to be included in order to
% update the output figure with newly computed values
idxSubjectArray = 2%2:7%[2:7];
for idxSubj = idxSubjectArray
    isFirstRun = idxSubj == idxSubjectArray(1);
    isLastRun = idxSubj == idxSubjectArray(end);
    
    % idxSubj, idxSess, idxRecon
    idxSubjSessReconArray = [
        idxSubj, 1, 1
        idxSubj, 2, 1
        idxSubj, 2, 2
        idxSubj, 1, 2
        ];
    
    for iConfig = [1]%nConfigs
        nFigures = numel(idxFiguresToPlot);
        fh = {};
        ft = {}; % table handle
        for iFigure = 1:nFigures
            idxFigure = idxFiguresToPlot(iFigure);
            iSubj = idxSubjSessReconArray(iConfig,1);
            iSess = idxSubjSessReconArray(iConfig,2);
            iRecon = idxSubjSessReconArray(iConfig,3);
            details = ...
                spifi_get_subject_details(iSubj, options, iSess, iRecon);
            
            switch idxFigure
                case 0
                    ft{1} = create_table_tsnr_spm_group(options.representation.tab_tsnr_spm_group);
                case 1
                    fh{1} = create_figure_trajectory_gradients(options.representation.fig_traj);
                case 2
                    figOpts = options.representation.fig_mean;
                    figOpts.iSubj = iSubj;
                    figOpts.iSess = iSess;
                    figOpts.iRecon = iRecon;
                    fh{2} = create_figure_mean(figOpts);
                case 3
                    figOpts = options.representation.fig_snr;
                    figOpts.iSubj = idxSubjSessReconArray(iConfig,1);
                    figOpts.iSess = iSess;
                    figOpts.iRecon = iRecon;
                    fh{3} = create_figure_snr(figOpts);
                case 4
                    figOpts = options.representation.fig_congruency;
                    figOpts.iSubj = idxSubjSessReconArray(iConfig,1);
                    figOpts.iSess = iSess;
                    figOpts.iRecon = iRecon;
                    fh{4} = create_figure_congruency(figOpts);
                case 5 % New Fig. 5, contour congruency, for reviewer 1
                    figOpts = details.representation.fig_specificity_contours;
                    figOpts.idxTissue = 1;
                    figOpts.iSubj = idxSubj;
                    [fh{5}, distMapX, figOpts] = compute_spatial_specificity_contour_congruency(figOpts);
                    
                    if isLastRun % need stats of all subjects from plot above to compute this
                        [fhTmp, roiArray, outputStats5] = create_figure_contour_congruency_summary_all_subjects();
                        fh{5} = [fh{5};fhTmp];
                    end
                case 6
                    figOpts = options.representation.fig_spm_subject;
                    figOpts.iSubj = idxSubjSessReconArray(iConfig,1);
                    figOpts.iSess = iSess;
                    figOpts.iRecon = iRecon;
                    if ~isfile(details.summary.mean.raw) % rerun creation of summary
                        spifi_run_analysis(figOpts.iSubj, figOpts.iSess, figOpts.iRecon, 12, options)
                    end
                    fh{6} = create_figure_spm_subject(figOpts);
                case 7
                    figOpts = options.representation.fig_spm_group;
                    figOpts.iSubjects = idxSubjectArray;
                    figOpts.iSess = iSess;
                    figOpts.iRecon = iRecon;
                    % since those are anyway plotted for for all subjects, do for first run of
                    % this loop only
                    if isFirstRun
                        fh{7} = create_figure_spm_group(figOpts);
                    end
                case 8 % New Fig. 8, activation overlap, for reviewer 1
                    figOpts = options.representation.fig_specificity_activation_tissue_overlap;
                    figOpts.iSubj = idxSubj;
                    figOpts.doUseSmoothedMaps = true;
                    figOpts.doUseVisualCortexMask = false;

                    [fh{8}, Y, figOpts] = compute_spatial_specificity_activation_tissue_overlap(figOpts);
                    
                    if isLastRun % need stats of all subjects from plot above to compute this
                        options.representation.fig_specificity_activation_tissue_overlap = figOpts;
                        [fhTmp, outputStats8] = create_figure_activation_tissue_overlap_summary_all_subjects(options);
                        fh{8} = [fh{8};fhTmp];
                    end
                case 9
                    figOpts = options.representation.fig_spm_inout;
                    figOpts.iSubj = iSubj;
                    %figOpts.iSess = iSess;
                    %figOpts.iRecon = iRecon;
                    fh{9} = create_figure_spm_inout(figOpts);
                case 11 % Figure S1 (supplementary)
                    % since those are anyway plotted for for all subjects, do for first run of
                    % this loop only
                    if isFirstRun
                        iSess = 1; % high res spiral out
                        iRecon = 1 % magnitude
                        fh{11} = create_figure_summary_results_all_subjects(iSess, iRecon);
                    end
                case 12 % Figure S2 (supplementary)
                    % since those are anyway plotted for for all subjects, do for first run of
                    % this loop only
                    if isFirstRun
                        iSess = 2; % spiral in/out
                        fh{12} = [];
                        for iRecon = 1:3 % magnitude
                            fhTmp = create_figure_summary_results_all_subjects(iSess, iRecon);
                            fh{12} = [fh{12};fhTmp];
                        end
                    end
                case 14 % Figure S4 (supplementary video smoothed/unsmoothed/cropped recon)
                    if isFirstRun
                        fh{14} = create_figure_activation_smoothed_unsmoothed_cropped();
                    end
                case 17 % Figure S7 (supplementary figure zoomed activation overlay smoothed/unsmoothed)
                        % also for figure 8A (smoothed/unsmoothed overlay
                        % mask)
                    if isFirstRun
                        fh{17} = create_figure_activation_overlay_smoothed_unsmoothed();
                    end
                    
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save Figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if doSaveFigures
        iFormatArray = 4; % 3 = highres, 4 = ultra-high res (slow)
        
        nFigures = numel(fh);
        for iFigure = 1:nFigures
            % save Figure-related plots in specific figure folders
            pathFigs = fullfile(paths.figures, sprintf('Figure%d', iFigure));
            % also works for vector of figure handles in fh{iFigure}!
            save_plots_abstract(fh{iFigure}, pathFigs, iFormatArray)
        end
    end
end

%% Also save raw subject tables, by sequence/recon combination
hasTable = ismember(0, idxFiguresToPlot);

if hasTable
    for iSpiralReconPair = 1:4
        idxRows = (1:6) + 6*(iSpiralReconPair-1);
        writetable(rows2vars(ft(idxRows,:)), fullfile(paths.figures, ...
            sprintf('Table1_raw_%d.xlsx', iSpiralReconPair)));
    end
end