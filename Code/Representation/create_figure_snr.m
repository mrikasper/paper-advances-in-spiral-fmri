function fh = create_figure_snr(figOptions)
% Creates plots of time series snr and sd (after realignment)
%
%   fh = create_figure_snr(figOptions)
%
% IN
%   figOptions  See also spifi_get_analysis_options:
%               options.representation.fig_...
% OUT
%   fh          vector of figure handles useful in the corresponding paper
%               figure
% EXAMPLE
%   create_figure_snr
%
%   See also
%
% Author:   Lars Kasper
% Created:  2018-09-19
% Copyright (C) 2018 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich

% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.

if nargin < 1
    options = spifi_get_analysis_options();
    figOptions = options.representation.fig_snr;
end

details = spifi_get_subject_details(figOptions.iSubj, [], ...
    figOptions.iSess, figOptions.iRecon);

files = {
    details.preproc.func.(figOptions.source)
    };
nFiles = numel(files);

%% parameter collection
cropX = details.representation.fig_snr.cropX;
cropY = details.representation.fig_snr.cropY;

nRows = min(6, round(numel(figOptions.selectedSlices)/4));
    
sharedParameters = {'plotLabels', false, 'plotTitle', false, 'rotate90', 1, ...
    'x', cropX(1):cropX(2), ...
    'y', cropY(1):cropY(2), ...
    'nRows', nRows, ...
    'colorBar', 'on'};


%% load files
for f = 1:nFiles
    Y{f} = MrImage(files{f});
    % masking for SNR plots, mask to avoid background effect in plots
    M = MrImage(regexprep(files{f},'_phase', ''), 'selectedVolumes',1);
    M = M.binarize(M.prctile(55)).imerode(strel('disk',2)).imfill('holes');
    if details.isPhase
        Y{f} = Y{f}.*M;
    end
end


fh = zeros(0,1);
for f = 1:nFiles
    % different ranges for magnitude and phase
    displayRangeSnr = figOptions.displayRangeSnr(1+details.isPhase,:);
    displayRangeSd = figOptions.displayRangeSd(1+details.isPhase,:);
    displayRangeCoeffVar = figOptions.displayRangeCoeffVar(1+details.isPhase,:);
    
    reducedY = Y{f}.select('invert', true, 't', figOptions.excludedVolumes);
    fh(end+1) = reducedY.snr.plot('z', figOptions.selectedSlices, 'displayRange', ...
        displayRangeSnr, sharedParameters{:});
    set(gcf, 'Name',  details.representation.fig_snr.fileSaveArray{1});
    axis off;
    
    fh(end+1) = reducedY.std.plot('z', figOptions.selectedSlices, 'displayRange', ...
       displayRangeSd, sharedParameters{:});
    set(gcf, 'Name',  details.representation.fig_snr.fileSaveArray{2});
    axis off;
    
    maskedY = reducedY.*M;
    fh(end+1) = maskedY.compute_coeff_var.plot('z', figOptions.selectedSlices, 'displayRange', ...
       displayRangeCoeffVar, sharedParameters{:});
    set(gcf, 'Name',  details.representation.fig_snr.fileSaveArray{3});
    axis off;
    
    
end