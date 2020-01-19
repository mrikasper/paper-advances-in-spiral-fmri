function [fh, Y] = create_figure_mean(figOptions)
% Creates plots of single volume and mean for spiral out and in/out time
% series
%
%   fh = create_figure_mean(figOptions)
%
% IN
%   figOptions  See also spifi_get_analysis_options:
%               options.representation.fig_...
% OUT
%   fh          vector of figure handles useful in the corresponding paper
%               figure
% EXAMPLE
%   create_figure_mean
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
    figOptions = options.representation.fig_mean;
end
details = spifi_get_subject_details(figOptions.iSubj, [], ...
    figOptions.iSess, figOptions.iRecon);

files = {
    details.preproc.func.(figOptions.source)
    };

nFiles = numel(files);

%% load files
for f = 1:nFiles
    Y{f} = MrImage(files{f});
    if details.isPhase % mask phase images
        M = MrImage(regexprep(files{f},'_phase', ''), 'selectedVolumes',1);
        M = M.binarize(M.prctile(55)).imerode(strel('disk',2)).imfill('holes');
        Y{f} = Y{f}.*M;
    end
    Y{nFiles+f} = Y{f}.mean;
end

nFiles = numel(Y);

cropX = details.representation.fig_mean.cropX;
cropY = details.representation.fig_mean.cropY;
if details.isPhase
    displayRange = [-1 2]; % [prctile(Y{nFiles},3), prctile(Y{nFiles},97)]; 
else
    displayRange = details.representation.fig_mean.displayRange;
end
sharedParameters = {'plotLabels', false, 'plotTitle', false, 'rotate90', 1};

%% different slice sets can be plotted (e.g. upper slices/lower slices
selectedSlices = figOptions.selectedSlices;
if ~iscell(selectedSlices)
    selectedSlices = {selectedSlices};
end

nSliceSets = numel(selectedSlices);

fh = zeros(nFiles,nSliceSets);
for f = 1:nFiles
    for iSliceSet = 1:nSliceSets
        nRows = min(6, round(numel(selectedSlices{iSliceSet})/4));
        fh(f,iSliceSet) = Y{f}.plot('y', [cropY(1):cropY(2)], 'x', [cropX(1):cropX(2)], ...
            't', 1, 'z', selectedSlices{iSliceSet}, 'displayRange', ...
            displayRange, 'nRows', nRows, sharedParameters{:});
        axis off
        set(gcf, 'Name',  sprintf('%s_sliceSet%d', ...
            details.representation.fig_mean.fileSaveArray{f}, iSliceSet));
    end
end