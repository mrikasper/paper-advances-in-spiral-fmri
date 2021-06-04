%% model specification
TR = 3.128
physio_path = ...
fileExplicitMask = '';% fullfile(data_path, 'masks', 'brain_mask.nii');

% timing in seconds
data.glm.timingUnits = 'secs';
% repetition time
data.glm.repetitionTime = TR;
% model derivatives
data.glm.hrfDerivatives = [0 0]; %[1 1];
% noise model FAST
data.glm.serialCorrelations = 'FAST';
% add conditions
% specify block length first
block_length = 18;
% specify first condition
first_condition = [1 5 7 11 13 17];
first_condition_onsets = first_condition*block_length;
% specify second condition
second_condition = [2 4 8 10 14 16];
second_condition_onsets = second_condition*block_length;
% add to glm
data.glm.conditions.names = {'simple', 'complex'};
data.glm.conditions.onsets = {second_condition_onsets, first_condition_onsets};
% add durations
data.glm.conditions.durations = {block_length, block_length};
% add confound regressors
confound_regressors = load(fullfile(physio_path,'multiple_regressors.txt'));      
data.glm.regressors.other = confound_regressors;
% add an explicit mask
data.glm.explicitMasking = fileExplicitMask;
% turn of inplicit masking threshold;
data.glm.maskingThreshold = -Inf;

% compute stat images
data.compute_stat_images;
% estimate
data.specify_and_estimate_1st_level;