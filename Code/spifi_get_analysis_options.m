function options = spifi_get_analysis_options()

idxPipeline = 2;
paths = spifi_get_paths(idxPipeline);

options.paths = paths;


% idxProcessingStepsArray
% doDeletePreviousResults = ismember(-1,idxProcessingStepsArray);
% doSetupPaths    = ismember(0, idxProcessingStepsArray);
% doCopyRawData   = ismember(1, idxProcessingStepsArray);
% doComputeBehav  = ismember(2, idxProcessingStepsArray);
% doPreprocFunctional = ismember(3, idxProcessingStepsArray);
% doCombineEchoes = idxSess==2 && idxRecon == 3 && doPreprocFunctional;
% doCoregAnatomical = ismember(4, idxProcessingStepsArray);
% doSegmentAnatomical = ismember(5, idxProcessingStepsArray);
% doBiasCorrectFunctional = ismember(6, idxProcessingStepsArray);
% doWarpImagesToMNI = ismember(7, idxProcessingStepsArray)
% doComputePhysio = ismember(8, idxProcessingStepsArray); % needs realignment params
% doComputeGlm    = ismember(9, idxProcessingStepsArray);
% doReportTaskContrasts = ismember(10, idxProcessingStepsArray);
% doReportPhysioContrasts = ismember(11, idxProcessingStepsArray);
% doComputeSummary = ismember(12, idxProcessingStepsArray);
% doCreateFigures = ismember(13, idxProcessingStepsArray);
options.idxProcessingStepsArray = [0:13];

% [1 2]; % 1:2 fMRI Out fmriInOut w/ Philips Achieva Gradients
% [3:6]; % 3:4 etaBeta 200/600; 5:6 etaBeta 100/1200
options.session.iTrajArray = [1, 2, 1, 2];
options.session.recons = {
     {'','phase'}
     {'In', 'Out', 'Combined', 'In_phase', 'Out_phase', 'Combined_phase'}
     {'', 'phase'}
     {'In', 'Out', 'Combined', 'In_phase', 'Out_phase', 'Combined_phase'}
     };
 
options.session.names = {...
    'fmri1', 'fmri2', ...
    'fmri3', 'fmri4'};

%% options for preproc and glm
options.pipelineId =  sprintf('Pipeline%02d', idxPipeline);

% input image for segmentation:
% T1: 't1';
% ME (multi-echo spin-warp): 'raw', 'mean', 'coreg'
options.preproc.segment.source = 'coreg';
options.preproc.realign.quality = 0.7;
options.preproc.realign.interpolationOrder = 3;
options.preproc.coreg.stateFuncForRef = 'realigned'; % raw, realigned, smoothed
options.preproc.biascorrect.source = 'mean_realigned'; % mean_realigned, coreg, raw, realigned, smoothed

options.glm.maskThreshold = 0.05; % very liberal do have no false negatives through masking
options.glm.doUsePhysio = true;
options.glm.doUseSmoothed = true;

% reported contrasts:
% 1 = ULLR, URLL, Vis FX, 
% 4 = Button Left, Right, All
options.glm.idxContrastArray = 1:8;
options.preproc.smooth.fwhm = [.8 .8 .8];

%% Options for representation (paper figures)

options.representation.fig_traj.titles = options.session.names;
options.representation.fig_traj.iTrajArray = options.session.iTrajArray;
options.representation.fig_traj.iRunArray = 1:2; % what will be plotted
options.representation.fig_traj.doRemoveAxes = true; % for copying into PPT

%% mean
% cropY{iSubj}(iSess,:) 
% crops per subject and session
cropY = {
    [15 255; 10 168]% 1
    [15 255; 10 168]% 2
    [25 265; 16 174]% 3
    [25 265; 10 168]% 4
    [15 255; 10 168]% 5
    [30 270; 20 178]% 6
    [15 255; 10 168]% 7
    };

% cropX{iSubj}(iSess,:) 
% crops per subject and session
cropX = {
    [25 220; 17 147]% 1
    [25 220; 17 147]% 2
    [20 215; 13 143]% 3
    [25 220; 17 147]% 4
    [25 220; 17 147]% 5
    [25 220; 17 147]% 6
    [25 220; 17 147]% 7
    };

% displayRange{iSess}{iRecon) 
% probably not subject specific, but recon and session?
displayRange = {
    [0.0 0.8; 0 pi*3/4] % fmri1; fmri1_phase
    [0.0 1.6; 0.0 1.6; 0.0 1.6; -pi pi; -pi pi; -pi pi] % fmri2_In;Out;Mean; then phase
};

fig_mean.iSubj = 2;
fig_mean.iSess = 1;
fig_mean.iRecon = 1;
fig_mean.source = 'raw'; % 'raw', 'realigned', 'smoothed'
fig_mean.cropX = cropX;
fig_mean.cropY = cropY;
fig_mean.selectedSlices = { 
    [10 17 24 31]
    [7  14 21 28]
    [1:3:36]
    [1:36]
        }; % 1:3:36
fig_mean.displayRange = displayRange;

options.representation.fig_mean = fig_mean;

%% snr
fig_snr.iSubj = 7;
fig_snr.iSess = 1;
fig_snr.iRecon = 1;
fig_snr.source = 'realigned'; % 'raw', 'realigned', 'smoothed'
fig_snr.selectedSlices = [10 17 24 31];
fig_snr.excludedVolumes = [1:5];
fig_snr.displayRangeSnr = [0 30; 0 10];
fig_snr.displayRangeSd = [0 0.05; 0 0.5];
fig_snr.displayRangeCoeffVar = [0 0.2; 0 0.8];
fig_snr.cropY = cropY;
fig_snr.cropX = cropX;

options.representation.fig_snr = fig_snr;

%% congruency
fig_congruency.iSubj = 2;
fig_congruency.iSess = 1;
fig_congruency.iRecon = 1;
% coreg uses rigid-body coregisteread mean ME anatomical as underlay
fig_congruency.source = 'coreg'; %'raw'; % 'raw', 'realigned', 'coreg'
fig_congruency.selectedSlice = 14; %18;
fig_congruency.displayRangeFunctional = displayRange;
fig_congruency.displayRangeStructural = [0.1 1.3];
% before edge creation, values below this threshold will be equalized to
% not show up in the edges (e.g., background voxels)
fig_congruency.thresholdEdgeOutliers = 0.15; % 0.005;
fig_congruency.thresholdEdge = 0.0015; % 0.005;
% if true, edge is not overlayed, but put in natural resolution on black background
fig_congruency.doCreateEdgeSeparately = true;
% if true, mean is plotted as underlay, if false, first volume is plotted
fig_congruency.doPlotFunctionalMean = true;
fig_congruency.cropY = cropY;
fig_congruency.cropX = cropX;

options.representation.fig_congruency = fig_congruency;

%% single subject SPM results
cropYVisualCortex = {
    [15 105; 10 168]% 1
    [15 105; 10 168]% 2
    [25 115; 16 174]% 3
    [15 105; 10 168]% 4
    [15 105; 10 168]% 5
    [20 110; 10 168]% 6
    [15 105; 10 168]% 7
    };
fig_spm_subject.iSubj = 2;
fig_spm_subject.iSess = 1;
fig_spm_subject.iRecon = 1;
fig_spm_subject.selectedSlice = [10 17 24 31];%[10 17 24 31 6 15 18 22]; %same slices as 4 and some more

% best slices
% [6 9 10 14 15 18 22 24];

%[
%     [10 17 24 31]-2
%     [10 17 24 31]
%     [10 17 24 31]+2
%     ]
    %1:36; %7:27%[16 22 28 34];%1:3:36;%17;
fig_spm_subject.displayRangeFunctional = displayRange;
fig_spm_subject.displayRangeStructural = [0.1 1.3];
% thresholds: p = [0.05 0.01 0.001 0.05(FWE)] 
%             t = [1.66 2.38 3.21  5.96]
fig_spm_subject.thresholdRange = [3.2 8];%[3 10];
fig_spm_subject.cropY = cropY; %cropYVisualCortex;
fig_spm_subject.cropX = cropX;
% choose either block wise or differential contrast ULRR vs URLL
% 'differential'    ULRR - URLL and URLL - ULLR contrasts ([1 -1], [-1 1])
% 'condition'       ULLR and URLL individual contrasts ([1 0], [0 1])
fig_spm_subject.contrastType = 'condition'; %'differential'; 'condition'
fig_spm_subject.views = {'tra_zoom'};%{'tra', 'tra_zoom', 'sag', 'sag_zoom', 'cor_zoom'};

options.representation.fig_spm_subject = fig_spm_subject;

%% SPM group results
options.representation.fig_spm_group.iSubjects = 2:7;
options.representation.fig_spm_group.selectionBlock = {'x', 0:100, 'y', 25:120};

% slice indices, e.g., [3:34]; 
%'max' or 'min' t value; 
% 'maxextent' - three cross-sections with most activated voxels
% 'maxextent_L' - maximum activation left hemisphere
% 'maxextent_R' - maximum activation right hemisphere
options.representation.fig_spm_group.selectedSlice = 'maxextent_R';
options.representation.fig_spm_group.views = {'sag'};%{'tra', 'tra_zoom', 'cor', 'sag'};
options.representation.fig_spm_group.crop_tra = {
    [1 240 1 292] 
    [30 220 20 255] 
    [32 210 28 266]
    [28 215 28 266] 
    [25 212 27 245] 
    [25 225 30 270] 
    [25 225 10 265] 
    };

options.representation.fig_spm_group.displayRangeFunctional = [0 0.8];
options.representation.fig_spm_group.displayRangeStructural = [0.1 1.3];
options.representation.fig_spm_group.thresholdRange = [3.2 8]; %[2 6];%[3 10];

%% comparison in/out/combined activation?
options.representation.fig_spm_inout.iSubj = 2;
options.representation.fig_spm_inout.iSess = 2; % 2 = InOut_1, 4 = InOut_2
options.representation.fig_spm_inout.selectedSlice =  [10 17 24 31]; %[10 17 24 31 6 15 18 22];
options.representation.fig_spm_inout.displayRangeIn = [0 1.4];
options.representation.fig_spm_inout.displayRangeOut = [0 1.4];
options.representation.fig_spm_inout.displayRangeCombined = [0 1.4];
options.representation.fig_spm_inout.displayRangeStructural = [0.1 1.3];
options.representation.fig_spm_inout.thresholdRange = [3.2 8];
options.representation.fig_spm_inout.views = {'tra', 'tra_zoom'};%{'tra', 'tra_zoom', 'cor', 'sag'};
% k-filter (raised cosine) for spiral-in underlay against ringing?
options.representation.fig_spm_inout.doKfilter = false;

%% Activation spatial specificity
options.representation.fig_specificity.selectionBlock = {'x', 0:100, 'y', 25:120, 'z', 18};
options.representation.fig_specificity.displayRange = [0 .8];
options.representation.fig_specificity.thresholdRange = [3.2 8];

%% Table with tSNR and dufferential t-statistics (SPM) summary for whole group
% 1 = high res spiral out;          2 = in-part spiral in/out
% 3 = out-part spiral in/out;       4 = combined images spiral in/out
tab_tsnr_spm_group.idxPairs         = 1:4; 
tab_tsnr_spm_group.iSubjects        = [2:7]; % subjects included in table
tab_tsnr_spm_group.doPlotMasks      = false; % Plot masks for ROI analysis as overlays
% slices used for ROI analysis to avoid margin effects cut-off slices of realignment
tab_tsnr_spm_group.selectedSlices   = 5:34; 
tab_tsnr_spm_group.clusterFormingThreshold = 3.2; % t value
tab_tsnr_spm_group.FWEcorrection = 'cluster'; % 'none'; 'cluster', 'peak'; %TODO: peak!
options.representation.tab_tsnr_spm_group = tab_tsnr_spm_group;

%% 
options.batches.preproc_glm = 'mb_preproc_glm.m';
options.batches.preproc = 'mb_preproc_quick.m';
options.batches.glm = 'mb_glm.m';
options.batches.segment = 'mb_segment.m';
options.batches.physio = 'mb_physio.m';
options.batches.warp2mni = 'mb_normalize.m';
