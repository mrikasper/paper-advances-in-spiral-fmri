function spifi_compute_summary(details)
%Computes summary images for later comparisons/figure creation (mean, snr)
%
%   spifi_compute_summary(details)
%
% IN
%
% OUT
%
% EXAMPLE
%   spifi_compute_summary
%
%   See also

% Author:   Lars Kasper
% Created:  2019-06-19
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

% copy the trivial files
[~,~] = mkdir(details.summary.root);

preprocStateArray = {'raw', 'realigned', 'smoothed'};
for iState = 1:numel(preprocStateArray)
    preprocState = preprocStateArray{iState};
    fileFunctional = details.preproc.func.(preprocState);
    Y = MrImage(fileFunctional);
    
    % mask phase images
    if details.isPhase
        fileMask = regexprep(fileFunctional, '_phase', '');
        M = MrImage(fileMask);
        M = mean(M, 't');
        M = M.binarize(prctile(M,60)).imfill('holes');
        Y = Y.*M;
    end
    
    Y.mean.save('fileName', details.summary.mean.(preprocState))
    Y.std.save('fileName', details.summary.std.(preprocState))
    Y.snr.save('fileName', details.summary.snr.(preprocState))
end

copyfile(details.glm.tcon, details.summary.tcon);


% copy warped and unbiased stuff, don't mind, if it does not exist
try
    copyfile(details.preproc.anat.warp2mni, details.summary.anat.warp2mni);
    copyfile(details.preproc.func.biascorrected, details.summary.func.biascorrected);
    copyfile(details.preproc.func.warp2mni, details.summary.func.warp2mni);
    copyfile(details.glm.warp2mni.tcon_differential1, ...
        details.summary.glm.warp2mni.tcon_differential1);
end
