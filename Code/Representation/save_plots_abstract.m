function save_plots_abstract(fh, pathFigs, iFormatArray)
% saves plots as wanted .fig and .png (using export_fig)
%
%   save_plots_abstract(fh, pathFigs)
%
% IN
%   iFormatArray, vector,   1 = .png
%                           2 = .fig
%                           3 = high-res (600 DPI) png
% OUT
%
% EXAMPLE
%   save_plots_abstract
%
%   See also
%
% Author:   Lars Kasper
% Created:  2016-11-10
% Copyright (C) 2016 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich

% This file is part of the Zurich fMRI Methods Evaluation Repository, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.

if nargin < 2
    pathFigs = pwd;
end

if nargin < 1
    fh = gcf;
end

if nargin < 3
    iFormatArray = [3 2];
end

saveAsFig = ismember(2, iFormatArray);
saveAsPng = ismember(1, iFormatArray);
saveAsHighResPng = ismember(3, iFormatArray);
saveAsHighestResPng = ismember(4, iFormatArray);

nFigs = numel(fh);
if nFigs>0 && ~exist(pathFigs, 'dir')
    mkdir(pathFigs);
end

for iFig = 1:nFigs
    currentFig = fh(iFig);
    stringTitle = str2fn(get(currentFig, 'Name'));
    
    if saveAsPng
        windowStyle = get(currentFig, 'WindowStyle');
        set(currentFig, 'WindowStyle', 'normal');
        export_fig(fullfile(pathFigs, stringTitle), ...
            '-png',  '-transparent', '-a1', '-q101', ...
            currentFig);
        set(currentFig, 'WindowStyle', windowStyle);
        
    end
    
    if saveAsFig
        save_fig('fh', currentFig, 'imageType', 'fig', 'pathSave', pathFigs, ...
            'doPrefixFigNumber', false);
    end
    
    if saveAsHighResPng
        tic
        windowStyle = get(currentFig, 'WindowStyle');
        set(currentFig, 'WindowStyle', 'normal');
        export_fig(fullfile(pathFigs, stringTitle), ...
            '-png',  '-transparent', '-a2', '-r300', ...
            currentFig);
        set(currentFig, 'WindowStyle', windowStyle);
        toc
    end
    
    if saveAsHighestResPng
        tic
        windowStyle = get(currentFig, 'WindowStyle');
        set(currentFig, 'WindowStyle', 'normal');
        export_fig(fullfile(pathFigs, stringTitle), ...
            '-png',  '-transparent', '-a4', '-q101', '-r600', ...
            currentFig);
        set(currentFig, 'WindowStyle', windowStyle);
        toc
    end
    
    
end