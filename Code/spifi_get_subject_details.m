function details = spifi_get_subject_details(iSubj, options, iSess, iRecon)
% subject-specific parameter options
% USE
%   details = spifi_get_subject_details(iSubj, options, iSess, iRecon)
% IN
%   iSubj       subject ID (1-7)
%   options     See also spifi_get_analysis_options, subject-independent
%               settings
%   iSess       session to be extracted (e.g. fMRIOut of fMRIInOut)
%               1   fMRIOut_1
%               2   fMRIInOut_1
%               3   fMRIOut_2
%               4   fMRIInOut_2
%   iRecon      reconstruction used (e.g. PartIn or PartOut of In/Out
%               spiral)
%               fmriOut
%               1 = magnitude
%               2 = phase
%               fmriInOut
%               1 = In          4 = In (phase)
%               2 = Out         5 = Out (phase)
%               3 = Combined    6 = Combined (phase)
if nargin < 4
    iRecon = 1;
end

if nargin < 3
    iSess = 2;
end

if nargin < 2 || isempty(options)
    options = spifi_get_analysis_options();
end

if nargin < 1
    iSubj = 2;
end

subjectId = sprintf('SPIFI_%04d', iSubj);

details.subjectId = subjectId;
details.reconId = options.session.recons{iSess}{iRecon};
details.isPhase = ~isempty(strfind(details.reconId ,'phase'));
details.session = options.session.names{iSess};
if isempty(details.reconId)
    details.recon = details.session;
else
    details.recon = sprintf('%s_%s', details.session, details.reconId);
end
details.files.func = [details.recon '.nii'];
details.files.anat = 'me1.nii';
details.files.physlog = [details.session '_phys.log'];
details.files.behav = [details.session '_behav.log'];

%% raw data (RawFmri folder)
details.paths.raw = fullfile(options.paths.data.root, subjectId);
details.paths.subject = fullfile(options.paths.results.root, subjectId);
details.paths.summary = fullfile(options.paths.results.summary , subjectId);
details.paths.scandata = fullfile(details.paths.subject, 'scandata');
details.paths.physlog = fullfile(details.paths.subject, 'physlog');
details.paths.behav = fullfile(details.paths.subject, 'behavior');

details.batches.root = fullfile(details.paths.subject, 'batches', details.recon);
details.batches.preproc_glm = fullfile(details.batches.root, options.batches.preproc_glm);
details.batches.preproc = fullfile(details.batches.root, options.batches.preproc);
details.batches.glm = fullfile(details.batches.root, options.batches.glm);
details.batches.physio = fullfile(details.batches.root, options.batches.physio);
details.batches.segment = fullfile(details.batches.root, options.batches.segment);
details.batches.warp2mni = fullfile(details.batches.root, options.batches.warp2mni);

details.behav.root = fullfile(details.paths.subject, 'behavior');
details.glm.root = fullfile(details.paths.subject, 'glm');

if options.glm.doUseSmoothed
    details.glm.session = fullfile(details.glm.root, ...
        sprintf('s%02d_%s',options.preproc.smooth.fwhm(1)*10, details.recon));
else
    details.glm.session = fullfile(details.glm.root, ...
        sprintf('unsmoothed_%s',details.recon));
end

details.glm.multcond = fullfile(details.behav.root, ...
    [details.session '_multiple_conditions.mat']);
details.glm.multreg = fullfile(details.paths.physlog, ...
    [details.session '_multiple_regressors.mat']);
details.glm.physio_mat = fullfile(details.paths.physlog, ...
    [details.session '_physio.mat']);
details.glm.contrasts_ps = fullfile(details.glm.session, ...
    [details.recon '_contrasts.ps']);
details.glm.spm_mat = fullfile(details.glm.session, 'SPM.mat');

details.glm.tcon_condition1 = fullfile(details.glm.session, 'spmT_0007.nii');
details.glm.tcon_condition2 = fullfile(details.glm.session, 'spmT_0008.nii');
details.glm.tcon_differential1 = fullfile(details.glm.session, 'spmT_0001.nii'); % differential (symmetric) contrast of conditions
details.glm.tcon_differential2 = fullfile(details.glm.session, 'spmT_0002.nii'); % differential (symmetric) contrast of conditions
details.glm.tcon = fullfile(details.glm.session, 'spmT_0001.nii'); % differential (symmetric) contrast of conditions
details.glm.fcon = fullfile(details.glm.session, 'spmF_0003.nii'); % Vis Fx F-contrast

details.preproc.root = details.paths.scandata;

details.preproc.func.raw = fullfile(details.preproc.root, ...
    details.files.func);

details.preproc.func.realigned = fullfile(details.preproc.root, ...
    ['ra' details.files.func]);
details.preproc.func.mean_realigned = fullfile(details.preproc.root, ...
    ['meana' details.files.func]);


details.preproc.rpfile = fullfile(details.preproc.root, ...
    sprintf('rp_a%s.txt', details.recon));

isStandardSmoothing = isequal(options.preproc.smooth.fwhm(1),0.8);

if isStandardSmoothing
    details.preproc.func.smoothed = fullfile(details.preproc.root, ...
        ['sra' details.files.func]);
else
    details.preproc.func.smoothed = fullfile(details.preproc.root, ...
        sprintf('s%02dra%s', options.preproc.smooth.fwhm(1)*10, details.files.func));
end


%% anatomical files before/after preproc
details.preproc.anat.raw = fullfile(details.preproc.root, ...
    details.files.anat);
details.preproc.anat.mean = fullfile(details.preproc.root, ...
    ['mean' details.files.anat]);

%% coregistration of anatomical
source = options.preproc.coreg.stateFuncForRef;
details.preproc.coreg.ref = details.preproc.func.(source);
details.preproc.coreg.source = details.preproc.anat.raw;
details.preproc.coreg.mean = details.preproc.anat.mean;

[~, refName] = fileparts(details.preproc.coreg.ref);
details.preproc.coreg.out = spm_file(details.preproc.anat.mean, ...
    'prefix', sprintf('cor2%s_', refName));

details.preproc.anat.coreg = details.preproc.coreg.out;

%% echo combination
% run details for sessions to be combine!
switch iRecon
    case 3 % echo combination
        detailsEcho1 = spifi_get_subject_details(iSubj, options, iSess,1);
        detailsEcho2 = spifi_get_subject_details(iSubj, options, iSess,2);
        details.preproc.combine_echoes.source = {
            detailsEcho1.preproc.func.raw
            detailsEcho2.preproc.func.raw
            };
    otherwise % no combination needed
        details.preproc.combine_echoes.source = {
            details.preproc.func.raw
            details.preproc.func.raw
            };
end
details.preproc.combine_echoes.out = details.preproc.func.raw;

%% segmentation/bias-correction of anatomical
source = options.preproc.segment.source;
switch lower(source)
    case {'raw', 'coreg', 'mean'}
        details.preproc.segment.source = details.preproc.anat.(source);
    case 't1'
        error('Segmenting %s not implemented yet', source);
end

details.preproc.anat.biascorrected = spm_file(details.preproc.segment.source, ...
    'prefix', 'm');
details.preproc.anat.biasfield = spm_file(details.preproc.segment.source, ...
    'prefix', 'BiasField_');
details.preproc.anat.tpm_GM = spm_file(details.preproc.segment.source, ...
    'prefix', 'c1');
details.preproc.anat.tpm_WM = spm_file(details.preproc.segment.source, ...
    'prefix', 'c2');
details.preproc.anat.tpm_CSF = spm_file(details.preproc.segment.source, ...
    'prefix', 'c3');
details.preproc.anat.forwardDeformationField = spm_file(details.preproc.segment.source, ...
    'prefix', 'y_');
details.preproc.anat.inverseDeformationField = spm_file(details.preproc.segment.source, ...
    'prefix', 'iy_');

%% bias correct functional
source = options.preproc.biascorrect.source;
details.preproc.biascorrect.source = details.preproc.func.(source);
details.preproc.biascorrect.biasfield = details.preproc.anat.biasfield;
details.preproc.biascorrect.out = spm_file(details.preproc.biascorrect.source, ...
    'prefix', 'm');

details.preproc.func.biascorrected = details.preproc.biascorrect.out;

%% warp to standard space (MNI)
details.preproc.warp2mni.templateBatch = options.batches.warp2mni;
details.preproc.warp2mni.saveBatch = details.batches.warp2mni;

details.preproc.warp2mni.source = {
    details.preproc.biascorrect.out
    details.preproc.anat.biascorrected
    details.glm.tcon_differential1
    };
details.preproc.warp2mni.forwardDeformationField = ...
    details.preproc.anat.forwardDeformationField;
details.preproc.warp2mni.out = cellfun(@(x) spm_file(x, 'prefix', 'w'), ...
    details.preproc.warp2mni.source, 'UniformOutput', false);

details.preproc.anat.warp2mni = details.preproc.warp2mni.out{1};
details.preproc.func.warp2mni = details.preproc.warp2mni.out{2};
details.glm.warp2mni.tcon_differential1 =  details.preproc.warp2mni.out{3};

%% transferred summary results
details.summary.root = details.paths.summary;

details.summary.firstfunc = fullfile(details.summary.root, ...
    regexprep(details.files.func, '\.nii', '_dyn0001\.nii'));

details.summary.tcon = fullfile(details.summary.root, ...
    regexprep(details.files.func, '\.nii', '_spmT_0001\.nii'));

details.summary.glm.warp2mni.tcon_differential1 = spm_file(...
    details.summary.tcon, 'prefix', 'warped_');
details.summary.warped_anat = fullfile(details.summary.root, ...
    'warped_struct.nii');

 details.summary.anat.warp2mni = fullfile(details.summary.root, ...
    ['warped_mean_' details.recon '.nii']);

details.summary.func.biascorrected = fullfile(details.summary.root, ...
    ['biascorrected_mean_' details.recon '.nii']);

 details.summary.func.warp2mni = fullfile(details.summary.root, ...
    ['warped_struct_coreg2_' details.recon '.nii']);

%% create file names for summary measures (mean/std/snr) of different processing states
preprocStateArray = {'raw', 'realigned', 'smoothed'};
for iState = 1:numel(preprocStateArray)
    preprocState = preprocStateArray{iState};
    details.summary.mean.(preprocState) = fullfile(details.summary.root, ...
        ['mean_' details.recon '_' preprocState '.nii']);
    details.summary.std.(preprocState) = fullfile(details.summary.root, ...
        ['std_' details.recon '_' preprocState '.nii']);
    details.summary.snr.(preprocState) = fullfile(details.summary.root, ...
        ['snr_' details.recon '_' preprocState '.nii']);
end

%% representation/subject-specific

%% fig_mean
fig_mean = options.representation.fig_mean;
details.representation.fig_mean.displayRange = ...
    fig_mean.displayRange{iSess}(iRecon,:);
details.representation.fig_mean.cropY = ...
    fig_mean.cropY{iSubj}(iSess,:);
details.representation.fig_mean.cropX = ...
    fig_mean.cropX{iSubj}(iSess,:);
details.representation.fig_mean.fileSaveArray = {
    [details.subjectId '_' details.recon '_' fig_mean.source '_vol0001']
    [details.subjectId '_' details.recon '_' fig_mean.source '_mean']
    };
details.representation.fig_mean.pathSave = fullfile(options.paths.figures, ...
    'FigureMean');

%% fig_snr
fig_snr = options.representation.fig_snr;
details.representation.fig_snr.cropY = ...
    fig_snr.cropY{iSubj}(iSess,:);
details.representation.fig_snr.cropX = ...
    fig_snr.cropX{iSubj}(iSess,:);
details.representation.fig_snr.fileSaveArray = {
    [details.subjectId '_' details.recon '_' fig_snr.source '_snr']
    [details.subjectId '_' details.recon '_' fig_snr.source '_sd']
    [details.subjectId '_' details.recon '_' fig_snr.source '_coeffVar']
    };
details.representation.fig_mean.pathSave = fullfile(options.paths.figures, ...
    'FigureMean');

%% fig_congruency
details.representation.fig_congruency.displayRangeFunctional = ...
    options.representation.fig_congruency.displayRangeFunctional{iSess}(iRecon,:);
details.representation.fig_congruency.cropY = ...
    options.representation.fig_congruency.cropY{iSubj}(iSess,:);
details.representation.fig_congruency.cropX = ...
    options.representation.fig_congruency.cropX{iSubj}(iSess,:);
details.representation.fig_congruency.fileSaveArray = {
    [details.subjectId '_' details.recon '_func']
    [details.subjectId '_' details.recon '_anat']
    [details.subjectId '_' details.recon '_overlay']
    };


%% fig_spm_subject
fig_spm_subject = options.representation.fig_spm_subject;

details.representation.fig_spm_subject.displayRangeFunctional = ...
    options.representation.fig_spm_subject.displayRangeFunctional{iSess}(iRecon,:);
details.representation.fig_spm_subject.cropY = ...
    fig_spm_subject.cropY{iSubj}(iSess,:);
details.representation.fig_spm_subject.cropX = ...
    fig_spm_subject.cropX{iSubj}(iSess,:);

details.representation.fig_spm_subject.pathSave = fullfile(options.paths.figures, ...
    'Figure5');
details.representation.fig_spm_subject.fileSaveArray = {
    [details.subjectId '_' details.recon '_meanfunc']
    [details.subjectId '_' details.recon '_anat']
    };

%% fig_spm_group
details.representation.fig_spm_group.fileSaveArray = {
    [details.subjectId '_' details.recon '_meanfunc']
    [details.subjectId '_' details.recon '_anat']
    };

%% fig_spm_inout
details.representation.fig_spm_inout.crop_tra = ...
     options.representation.fig_spm_group.crop_tra{iSubj};
details.representation.fig_spm_inout.fileSaveArray = {
    [details.subjectId '_' details.recon '_meanfunc']
    [details.subjectId '_' details.recon '_anat']
    };

[~,~] = mkdir(details.behav.root);
[~,~] = mkdir(details.glm.session);
[~,~] = mkdir(details.preproc.root);
[~,~] = mkdir(details.paths.summary);
[~,~] = mkdir(details.batches.root);