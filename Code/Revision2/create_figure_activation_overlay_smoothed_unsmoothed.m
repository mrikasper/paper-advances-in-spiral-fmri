function fh = create_figure_activation_overlay_smoothed_unsmoothed()
% Creates zoomed sagittal/coronal/transverse views of activation related to
% different reconstruction (cropping), processing (smoothing) and thresholding
% choices, overlayed on anatomical image (mean ME)
%
%   output = create_figure_activation_overlay_smoothed_unsmoothed()
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
idxSubjectArray = [2];%[2 12];

% set true for overlay cropped/uncropped and select uncropped subject (2)
% above
doCompareToCropped1mm = true;%true;

iSess = 1;
iRecon = 1;

% doSmoothing, activationThreshold, FWECorrection, doCreateEdgesSeparately, underlays
parameterArray = {
    false, 2.38,    'none',     false,  'none'
    false, 2.38,    'none',     true,  'none'
    false, 3.2,     'cluster',  true,  'none'
    false, 3.2,     'cluster',  false,  'none'
    true,  3.2,     'cluster',  true,   'structural'
    false, 3.2,     'none',  false,  'structural'
    false, [3.2 3.2 2.38 2.38],     {'cluster', 'cluster', 'none', 'none'}, ...
    false,  'structural'
    };


nParameterConfigs = size(parameterArray, 1);
fh = [];
for idxSubj = idxSubjectArray
    for iConfig = 5%nParameterConfigs
        doSmoothing             = parameterArray{iConfig, 1};
        activationThreshold     = parameterArray{iConfig, 2};
        FWEcorrection           = parameterArray{iConfig, 3};
        doPlotMask              = parameterArray{iConfig, 4};
        underlays               = parameterArray{iConfig, 5};
        
        options = spifi_get_analysis_options();
        
        fig_spm_subject = options.representation.fig_spm_subject;
        fig_spm_subject.iSubj = idxSubj;
        fig_spm_subject.iSess = iSess;
        fig_spm_subject.iRecon = iRecon;
        
        fig_spm_subject.selectedSlice = 3:34;
        
        % thresholds: p = [0.05 0.01 0.001 0.05(FWE)]
        %             t = [1.66 2.38 3.21  5.96]
        if numel(activationThreshold) > 1
            fig_spm_subject.thresholdRange = activationThreshold(:);
            fig_spm_subject.thresholdRange(:,2) = 8;
        else
            fig_spm_subject.thresholdRange = [activationThreshold 8];%[2 8]%[3 10];
        end
        % choose either block wise or differential contrast ULRR vs URLL
        % 'differential'    ULRR - URLL and URLL - ULLR contrasts ([1 -1], [-1 1])
        % 'condition'       ULLR and URLL individual contrasts ([1 0], [0 1])
        fig_spm_subject.contrastType = 'differential'; %'differential'; 'condition'
        fig_spm_subject.underlays = underlays;
        fig_spm_subject.views = {'tra_zoom'};%{'tra', 'tra_zoom', 'sag', 'sag_zoom', 'cor_zoom'};
        fig_spm_subject.nRows = 8;
        fig_spm_subject.FWEcorrection = FWEcorrection; % 'none'; 'cluster', 'peak'; %TODO: peak!
        fig_spm_subject.doUseSmoothedMaps = doSmoothing;
        fig_spm_subject.doPlotMask = doPlotMask;
        fig_spm_subject.displayRangeStructural = [0.2 1.3];
        %fhTmp = create_figure_spm_subject(fig_spm_subject);
        fig_spm_subject.underlays = 'none';
        fig_spm_subject.doCompareToCropped1mm = doCompareToCropped1mm;
        fhTmp = create_figure_spm_subject_smoothed_unsmoothed(fig_spm_subject);
        fh = [fh; fhTmp];
    end
end