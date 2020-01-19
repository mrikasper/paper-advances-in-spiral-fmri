function fh = create_figure_congruency(figOptions)
% Creates plots of spiral functional mean, structural image and overlay edge of structural
% NOTE: For multi-echo data, mean of all echo images is presented
%
%   fh = create_figure_congruency(figOptions)
%
% IN
%   figOptions  See also spifi_get_analysis_options:
%               options.representation.fig_...
% OUT
%   fh          vector of figure handles useful in the corresponding paper
%               figure
% EXAMPLE
%   create_figure_congruency
%
%   See also
%
% Author:   Lars Kasper
% Created:  2018-09-20
% Copyright (C) 2018 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich

% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.

iDimEcho = 4;

if nargin < 1
    options = spifi_get_analysis_options();
    figOptions = options.representation.fig_congruency;
end
details = spifi_get_subject_details(figOptions.iSubj, [], ...
    figOptions.iSess, figOptions.iRecon);

switch figOptions.source
    case {'raw', 'realigned'}
        files = {
            details.preproc.func.(figOptions.source)
            details.preproc.anat.coreg
            };
    case 'coreg'
        files = {
            details.preproc.coreg.ref
            details.preproc.anat.coreg
            };
end

nFiles = numel(files);

nameArray = {'func', 'anat', 'overlay', 'edge'};

cropXRef = figOptions.cropX{figOptions.iSubj}(figOptions.iSess,:);
cropYRef = figOptions.cropY{figOptions.iSubj}(figOptions.iSess,:);
%% load files
for f = 1:nFiles
    Y{f} = MrImage(files{f});
    if f == 1, Yref = Y{1}.select('t',1).copyobj(); end % to not crop twice
    [cropX, cropY] = adjust_crop(cropXRef, cropYRef, Yref, Y{f});
    Y{f} = Y{f}.select('x', cropX(1):cropX(2), 'y', cropY(1):cropY(2));
end

doPlotSeparate = true;
if doPlotSeparate
    fh = zeros(3*(nFiles-1),1);
else
    fh = zeros(nFiles-1,1);
end

displayRangeFunctional = details.representation.fig_congruency.displayRangeFunctional;
sharedParameters = {'plotLabels', false, 'plotTitle', false, 'rotate90', 1};
for f = 1:nFiles-1
    if figOptions.doPlotFunctionalMean
        Y{f} = Y{f}.mean('t');
        nameArray{1} = 'meanfunc';
    end
    
    if doPlotSeparate
        fh(f+(nFiles-1)) = Y{nFiles}.mean(iDimEcho).plot('z', ...
            figOptions.selectedSlice, 'displayRange', ...
            figOptions.displayRangeStructural, ...
            sharedParameters{:});
        axis off
        set(gcf, 'name', sprintf('%s %s %s', details.subjectId, details.recon, nameArray{2})),
        fh(f+2*(nFiles-1)) = Y{f}.select('t',1).plot('z', ...
            figOptions.selectedSlice, 'displayRange', ...
            displayRangeFunctional, sharedParameters{:});
        axis off;
        set(gcf, 'name', sprintf('%s %s %s', details.subjectId, details.recon, nameArray{1})),
    end
    if figOptions.doCreateEdgeSeparately
        % dummy image with all zeros but same size as original anat image
        Z = Y{nFiles}.select('t',1).copyobj;
        Z.data(:) = 0;
        fh(f) = Z.apply_threshold(displayRangeFunctional).plot_overlays(...
            mean(Y{nFiles}.apply_threshold(figOptions.thresholdEdgeOutliers), iDimEcho), 'selectedSlices', figOptions.selectedSlice, 'overlayMode',  'edge', ...
            'edgeThreshold', figOptions.thresholdEdge, ...
            'colorBar', 'off', sharedParameters{:});
        set(gcf, 'name', sprintf('%s %s %s', details.subjectId, details.recon, nameArray{4})),
    else
        fh(f) = Y{f}.select('t',1).apply_threshold(displayRangeFunctional).plot_overlays(...
            mean(Y{nFiles}.apply_threshold(figOptions.thresholdEdgeOutliers), iDimEcho), 'selectedSlices', figOptions.selectedSlice, 'overlayMode',  'edge', ...
            'edgeThreshold', figOptions.thresholdEdge, ...
            'colorBar', 'off', sharedParameters{:});
        set(gcf, 'name', sprintf('%s %s %s', details.subjectId, details.recon, nameArray{3})),
    end
    axis off
end