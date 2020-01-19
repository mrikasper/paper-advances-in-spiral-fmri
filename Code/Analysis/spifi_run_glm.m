function matlabbatch = spifi_run_glm(fileMultipleConditions, fileMultipleRegressors, ...
    pathPreproc, pfxNiftiPreproc, scanInfo, pathGlm, templateBatch, savedBatch, ...
    doCreateContrastReportWithJob, maskThreshold)
% Sets up and estimates GLM and creates task contrasts
%
%  matlabbatch = spifi_run_glm(fileMultipleConditions, fileMultipleRegressors, ...
%        pathPreproc, pfxNiftiPreproc, scanInfo, pathGlm, templateBatch, savedBatch, ...
%        doCreateContrastReportWithJob);
%
% IN
%
% OUT
%
% EXAMPLE
%   spifi_run_glm
%
%   See also

% Author:   Lars Kasper
% Created:  2019-06-16
% Copyright (C) 2019 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
%

% delete existing directory
if isfolder(pathGlm)
    rmdir(pathGlm, 's')
end

run(templateBatch)
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = scanInfo.TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = scanInfo.nSlices;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = round(scanInfo.nSlices/2);
matlabbatch{1}.spm.stats.fmri_spec.dir = {pathGlm};
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(...
    spm_select('ExtFPList', pathPreproc, pfxNiftiPreproc, Inf));
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {fileMultipleConditions};

matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {fileMultipleRegressors};
matlabbatch{1}.spm.stats.fmri_spec.mthresh = maskThreshold;

if ~doCreateContrastReportWithJob
    matlabbatch(4) = [];
end
% TODO with timestamp!
tnufmri_convert_batch_mat_to_m(matlabbatch,savedBatch)
spm_jobman('interactive', matlabbatch);
spm_jobman('run', matlabbatch);