function fh = create_figure_activation_smoothed_unsmoothed_cropped()
% Creates zoomed sagittal/coronal/transverse views of activation related to
% different reconstruction (cropping), processing (smoothing) and thresholding
% choices, overlayed on anatomical image (mean ME)
%
%   output = create_figure_activation_smoothed_unsmoothed_cropped(input)
%
% IN
%
% OUT
%
% EXAMPLE
%   create_figure_activation_smoothed_unsmoothed_cropped
%
%   See also

% Author:   Lars Kasper
% Created:  2021-10-10
% Copyright (C) 2021 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
%

% 2 = full resolution, 12 = cropped spiral to 1mm, zero-filled
idxSubjectArray = [2 12];
doSmoothingArray = [false, true];
activationThresholdArray = [2.38 3.2];

iSess = 1;
iRecon = 1;

fh = [];
for idxSubj = idxSubjectArray
    for doSmoothing = doSmoothingArray
        for activationThreshold = activationThresholdArray
            options = spifi_get_analysis_options();
           
            fig_spm_subject = options.representation.fig_spm_subject;
            fig_spm_subject.iSubj = idxSubj;
            fig_spm_subject.iSess = iSess;
            fig_spm_subject.iRecon = iRecon;    
            
            % thresholds: p = [0.05 0.01 0.001 0.05(FWE)]
            %             t = [1.66 2.38 3.21  5.96]
            fig_spm_subject.thresholdRange = [activationThreshold 8];%[2 8]%[3 10];
            % choose either block wise or differential contrast ULRR vs URLL
            % 'differential'    ULRR - URLL and URLL - ULLR contrasts ([1 -1], [-1 1])
            % 'condition'       ULLR and URLL individual contrasts ([1 0], [0 1])
            fig_spm_subject.contrastType = 'differential'; %'differential'; 'condition'
            fig_spm_subject.underlays = 'structural';
            fig_spm_subject.views = {'tra_zoom', 'sag', 'cor_zoom'};%{'tra', 'tra_zoom', 'sag', 'sag_zoom', 'cor_zoom'};
            fig_spm_subject.doUseSmoothedMaps = doSmoothing;
            fhTmp = create_figure_spm_subject(fig_spm_subject);
            fh = [fh; fhTmp];
        end
    end
end