function spifi_run_analysis(idxSubj, idxSess, idxRecon, idxProcessingStepsArray, options)
% Function spifi_run_analysis
% fMRI preproc and statistical analysis for a single subject, run and recon
% choice. Different preprocessing and statistical analysis steps may be
% defined
%
% USE
%  spifi_run_analysis(idxSubj, idxSess, idxRecon, idxProcessingStepsArray, options)
%
% IN
%   idxSubj     subject ID (1-7)
%   idxSess     session to be extracted (e.g. fMRIOut of fMRIInOut)
%               1   fMRIOut_1   (in paper high-res spiral out 0.8 mm)
%               2   fMRIInOut_1 (in papaer spiral in/out 1.5 mm)
%               3   fMRIOut_2
%               4   fMRIInOut_2
%   idxRecon    reconstruction used (e.g. PartIn or PartOut of In/Out
%               spiral)
%               fmriOut
%               1 = magnitude
%               2 = phase
%               fmriInOut
%               1 = In          4 = In (phase)
%               2 = Out         5 = Out (phase)
%               3 = Combined    6 = Combined (phase)
%   idxProcessingStepsArray
%               index vector, e.g., [1:12] of processing steps to be
%               executed
%                doDeletePreviousResults    -1 
%                doSetupPaths               0
%                doCopyRawData              1
%                doComputeBehav             2
%                doPreprocFunctional        3
%                doCombineEchoes = idxSess==2 && idxRecon == 3 && doPreprocFunctional;
%                doCoregAnatomical          4
%                doSegmentAnatomical        5
%                doBiasCorrectFunctional	6
%                doWarpImagesToMNI          7
%                doComputePhysio            8 % needs realignment params
%                doComputeGlm               9
%                doReportTaskContrasts      10
%                doReportPhysioContrasts 	11
%                doComputeSummary           12
%                doCreateFigures            13
%   options     See also spifi_get_analysis_options, subject-independent
%               settings

% Author:   Lars Kasper
% Created:  2019-06-07
% Copyright (C) 2019 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.

if nargin < 5 || isempty(options)
    options = spifi_get_analysis_options();
end
details = spifi_get_subject_details(idxSubj, options, idxSess, idxRecon);
paths = details.paths;
files = details.files;

if nargin < 4
    idxProcessingStepsArray = options.idxProcessingStepsArray;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% idxProcessingStepsArray
doDeletePreviousResults = ismember(-1,idxProcessingStepsArray);
doSetupPaths    = ismember(0, idxProcessingStepsArray);
doCopyRawData   = ismember(1, idxProcessingStepsArray);
doComputeBehav  = ismember(2, idxProcessingStepsArray);
doPreprocFunctional = ismember(3, idxProcessingStepsArray);
doCombineEchoes = idxSess==2 && idxRecon == 3 && doPreprocFunctional;
doCoregAnatomical = ismember(4, idxProcessingStepsArray);
doSegmentAnatomical = ismember(5, idxProcessingStepsArray);
doBiasCorrectFunctional = ismember(6, idxProcessingStepsArray);
doWarpImagesToMNI = ismember(7, idxProcessingStepsArray);
doComputePhysio = ismember(8, idxProcessingStepsArray); % needs realignment params
doComputeGlm    = ismember(9, idxProcessingStepsArray);
doReportTaskContrasts = ismember(10, idxProcessingStepsArray);
doReportPhysioContrasts = ismember(11, idxProcessingStepsArray);
doComputeSummary = ismember(12, idxProcessingStepsArray);
doCreateFigures = ismember(13, idxProcessingStepsArray);


scanInfo = spifi_get_scaninfo(idxSess);

%%
nameSess = details.session;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup paths and spm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doSetupPaths
    spifi_setup_paths();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Delete previous results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doDeletePreviousResults
    rmdir(paths.subject, 's');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copy Raw Data for results subject folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doCopyRawData
    copyfile(paths.raw, paths.subject);
    pause(10);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute behavioral regressors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileMultipleConditions = details.glm.multcond;
fileBehav = fullfile(paths.behav, files.behav);
if doComputeBehav
    spifi_compute_behav(fileBehav, fileMultipleConditions, scanInfo);    
end


%% Echo combination for in/out sessions
if doCombineEchoes
   spifi_combine_echoes(details.preproc.combine_echoes) 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathPreproc = details.preproc.root;
preprocOptions = options.preproc;
templateBatch = options.batches.preproc;
savedBatch = details.batches.preproc;
pfxFmri = ['^' details.recon '\.'];

if doPreprocFunctional
    matlabbatch = spifi_run_preproc(pathPreproc, pfxFmri, scanInfo, ...
        preprocOptions, templateBatch, savedBatch);
end


% co-register mean of ME images to mean fMRI and reslice
if doCoregAnatomical
    spifi_coreg_anat_to_func(details.preproc.coreg);
end


% co-register mean of ME images to mean fMRI and reslice
if doSegmentAnatomical
    fileAnat = details.preproc.segment.source;
    templateBatch = options.batches.segment;
    savedBatch = details.batches.segment;
    matlabbatch = spifi_run_segment_anat(fileAnat, templateBatch, savedBatch);
end

% needs bias field from segmentation
if doBiasCorrectFunctional
    spifi_biascorrect_func(details.preproc.biascorrect);
end

% needs forward deformation (warp) field from segmentation
if doWarpImagesToMNI
   spifi_warp2mni_func_anat(details.preproc.warp2mni);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determining multiple regressors w/ or w/o Physiological noise correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filePhysio = fullfile(paths.physlog, [nameSess '_phys.log']);
fileRp = details.preproc.rpfile;
fileMultipleRegressors = details.glm.multreg;
filePhysioOut = details.glm.physio_mat;
templateBatch = options.batches.physio;
savedBatch = details.batches.physio;

if doComputePhysio
    matlabbatch = spifi_run_physio(filePhysio, fileRp, fileMultipleRegressors, ...
        filePhysioOut, templateBatch, savedBatch);    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SPM glm analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% will be done by tapas_physio_overlay_contrasts below
doCreateContrastReportWithJob = false;

if ~options.glm.doUsePhysio
    fileMultipleRegressors = details.preproc.rpfile;
end

pathGlm = details.glm.session;

if options.glm.doUseSmoothed
    [~,pfxNiftiPreproc,~] = fileparts(details.preproc.func.smoothed);
else
    [~,pfxNiftiPreproc,~] = fileparts(details.preproc.func.realigned);
end
pfxNiftiPreproc = ['^' pfxNiftiPreproc '\.'];
templateBatch = options.batches.glm;
savedBatch = details.batches.glm;
maskThreshold = options.glm.maskThreshold;

if doComputeGlm
    matlabbatch = spifi_run_glm(fileMultipleConditions, fileMultipleRegressors, ...
       pathPreproc, pfxNiftiPreproc, scanInfo, pathGlm, templateBatch, savedBatch, ...
       doCreateContrastReportWithJob, maskThreshold);
end

%% Report nuisance and task contrasts using PhysIO Toolbox
if doReportTaskContrasts && ~doCreateContrastReportWithJob
     args = tapas_physio_overlay_contrasts(...
        'fileReport', details.glm.contrasts_ps, ...
        'fileSpm', details.glm.spm_mat, ...
        'idxContrasts', options.glm.idxContrastArray, ...
        'fileStructural', details.preproc.anat.biascorrected, ...
        'threshold', 0.001, ...
        'correction', 'none')   
end

if doReportPhysioContrasts
    args = tapas_physio_report_contrasts(...
        'fileReport', details.glm.contrasts_ps, ...
        'fileSpm', details.glm.spm_mat, ...
        'filePhysIO', details.glm.physio_mat, ...
        'fileStructural', details.preproc.anat.biascorrected, ...
        'reportContrastThreshold', 0.001, ...
        'reportContrastCorrection', 'none')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute summary images for faster plotting and data transfer later on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if doComputeSummary
    spifi_compute_summary(details);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create figure (parts) for publication
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if doCreateFigures
    main_create_figures();
end