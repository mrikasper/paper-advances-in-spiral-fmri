function matlabbatch = spifi_run_physio(filePhysio, fileRp, fileMultipleRegressors, ...
    filePhysioOut, templateBatch, savedBatch)
% Computes multile regressors file (nuisance) from physiological recordings
% and realignment parameters using PhysIO toolbox
%
%   matlabbatch = spifi_run_physio(filePhysio, fileRp, fileMultipleRegressors, ...
%         filePhysioOut, templateBatch, savedBatch)
%
% IN
%
% OUT
%
% EXAMPLE
%   spifi_run_physio
%
%   See also

% Author:  Lars Kasper
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
[fp,fn,ext] = fileparts(filePhysioOut);
filePhysioOut = [fn ext];
run(templateBatch);
matlabbatch{1}.spm.tools.physio.log_files.cardiac = {filePhysio};
matlabbatch{1}.spm.tools.physio.log_files.respiration = {filePhysio};
matlabbatch{1}.spm.tools.physio.save_dir = {fp};
matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = ...
    fileMultipleRegressors;
matlabbatch{1}.spm.tools.physio.model.output_physio = ...
    filePhysioOut;
matlabbatch{1}.spm.tools.physio.model.movement.yes.file_realignment_parameters = ...
    {fileRp};
tnufmri_convert_batch_mat_to_m(matlabbatch, savedBatch)
spm_jobman('interactive', matlabbatch);
spm_jobman('run', matlabbatch);
