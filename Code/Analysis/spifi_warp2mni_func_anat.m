function [meanX, biasCorrectedX] = spifi_warp2mni_func_anat(opts)
% Biascorrects mean functional image by applying previously calculated
% biasfield of segmentation (of anatomical image)
%
%   [meanX, coregX] = spifi_warp2mni_func_anat(opts)
%
% IN
%   opts   details.preproc.warp2mni
% OUT
%
% EXAMPLE
%   spifi_warp2mni_func_anat
%
%   See also

% Author:   Lars Kasper
% Created:  2019-11-15
% Copyright (C) 2019 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich

run(opts.templateBatch);

matlabbatch{1}.spm.spatial.normalise.write.subj.def{1} = ...
    opts.forwardDeformationField;
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = ...
    opts.source(isfile(opts.source)); % only warp existing files
tnufmri_convert_batch_mat_to_m(matlabbatch, opts.saveBatch)
spm_jobman('interactive', matlabbatch);
spm_jobman('run', matlabbatch);