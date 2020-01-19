function matlabbatch = spifi_run_preproc(pathPreproc, pfxFmri, scanInfo, ...
    preprocOptions, templateBatch, fileSavedBatch)
% runs preprocessing of functional data from existing template batch (slice
% timing, realignment and smoothing
%
%  matlabbatch = spifi_run_preproc(pathPreproc, pfxFmri, scanInfo, ...
%         preprocOptions, templateBatch, fileSavedBatch);
%
% IN
%
% OUT
%
% EXAMPLE
%   spifi_run_preproc
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

isPhase = find_string(pfxFmri, '_phase');

run(templateBatch);

matlabbatch{1}.spm.temporal.st.scans = {cellstr(...
    spm_select('ExtFPList', pathPreproc, pfxFmri, Inf))};
matlabbatch{1}.spm.temporal.st.nslices = scanInfo.nSlices;
matlabbatch{1}.spm.temporal.st.tr = scanInfo.TR;
matlabbatch{1}.spm.temporal.st.ta = scanInfo.TR*(1-1/scanInfo.nSlices);
matlabbatch{1}.spm.temporal.st.so = [scanInfo.nSlices:-1:1];
matlabbatch{1}.spm.temporal.st.refslice = round(scanInfo.nSlices/2);

matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = ...
    preprocOptions.realign.interpolationOrder;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = ...
    preprocOptions.realign.quality;

% realign magnitude and apply to phase!
if isPhase
    %% realign magnitude instead and apply to magnitude, but use existing ...
    % rp_files, if they exist
    pfxFmriMagn = regexprep(pfxFmri, '_phase', '');
    
    rpFileMagn = fullfile(pathPreproc, ['rp_a' pfxFmriMagn(2:end-2)' '.txt']);
    if ~exist(rpFileMagn)
        matlabbatch{2}.spm.spatial.realign.estwrite.data = ...
            {cellstr(spm_select('ExtFPList', pathPreproc, pfxFmriMagn, Inf))};
        matlabbatch(3) = [];
    else
        matlabbatch(2:3) = [];
    end
    
    tnufmri_convert_batch_mat_to_m(matlabbatch, fileSavedBatch)
    spm_jobman('interactive', matlabbatch);
    spm_jobman('run', matlabbatch);
    
    
    %% do smoothing separately afterwards
    run(templateBatch);
    
    if preprocOptions.smooth.fwhm > 0
        matlabbatch{3}.spm.spatial.smooth.fwhm = preprocOptions.smooth.fwhm;
    else
        % remove smoothing completely
        matlabbatch(3) = [];
    end
    
    matlabbatch{3}.spm.spatial.smooth.data = ...
        {cellstr(spm_select('ExtFPList', pathPreproc, ['^ra' pfxFmri(2:end)], Inf))};
    
    matlabbatch(1:2) = [];
    
    tnufmri_convert_batch_mat_to_m(matlabbatch, fileSavedBatch)
    spm_jobman('interactive', matlabbatch);
    spm_jobman('run', matlabbatch);
    
else
    
    if preprocOptions.smooth.fwhm > 0
        matlabbatch{3}.spm.spatial.smooth.fwhm = preprocOptions.smooth.fwhm;
    else
        % remove smoothing completely
        matlabbatch(3) = [];
    end
    
    tnufmri_convert_batch_mat_to_m(matlabbatch, fileSavedBatch)
    spm_jobman('interactive', matlabbatch);
    spm_jobman('run', matlabbatch);
end