function matlabbatch = spifi_run_segment_anat(fileAnat, templateBatch, savedBatch)
% Segments anatomical image and saves bias-corrected version, as well as
% bias field and tissue probability maps
%
%   matlabbatch = spifi_run_segment_anat(fileAnat, templateBatch, savedBatch)
%
% IN
%
% OUT
%
% EXAMPLE
%   spifi_run_segment_anat
%
%   See also

% Author:  Lars Kasper
% Created:  2019-08-25
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
run(templateBatch);
matlabbatch{1}.spm.spatial.preproc.channel.vols{1} = fileAnat;
tnufmri_convert_batch_mat_to_m(matlabbatch, savedBatch);
spm_jobman('interactive', matlabbatch);
spm_jobman('run', matlabbatch);
