function [qualityMeasuresArray, fh] = spifi_assess_motion(idxSubjArray)
% Assesses motion quality parameters for different subjects using PhysIO
% and the rp*.txt file from SPM
%
%   [qualityMeasuresArray, fh] = spifi_assess_motion(idxSubjArray)
%
% IN
%
% OUT
%
% EXAMPLE
%   spifi_assess_motion()
%
%   See also

% Author:   Lars Kasper
% Created:  2021-09-17
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
if nargin < 1
    idxSubjArray = 2:7;
end

options = spifi_get_analysis_options();

doPlot = true;
doSavePlot = true;


nSubjects = numel(idxSubjArray);

for iSubj = 1:nSubjects
    idxSubj = idxSubjArray(iSubj);
    details = spifi_get_subject_details(idxSubj, options);
    
    physio = tapas_physio_new;
    movement = physio.model.movement;
    verbose = physio.verbose;
    
    movement.file_realignment_parameters = details.preproc.rpfile;
    
    [R, movement, verbose] = tapas_physio_create_movement_regressors(...
        movement, verbose);
    
    if doPlot
        fh = tapas_physio_plot_movement_outliers_fd(R, movement.quality_measures, ...
            movement.censoring, movement.censoring_threshold);
        
        
        if iSubj == 1
            fh2 = tapas_physio_get_default_fig_params();
            fh3 = tapas_physio_get_default_fig_params();
        end
        
        %% Figure with realignment-parameters
        figure(fh2);
        hs2(iSubj) = subplot(3,2,iSubj);
        copyobj(fh.Children(end).Children, hs2(iSubj))
        ylim([-1 1]*1.0);
        
        if iSubj == 1
            legend();
        end
        title(['Subject ' int2str(idxSubj)]);
        
        if iSubj == nSubjects
            try % only for R2019b and newer, but not critical to plot
                sgtitle(fh.Children(end).Title.String);
            end
            set(fh2, 'Name',fh.Children(end).Title.String);
        end
        
        if iSubj >= nSubjects -1
            xlabel('Scan Volumes');
        end
        
        %% Figure with framewise displacement
        figure(fh3);
        hs3(iSubj) = subplot(3,2,iSubj);
        copyobj(fh.Children(end-2).Children, hs3(iSubj))
        ylim([0 1]*2.0);
        
        if iSubj == 1
            legend();
        end
        title(sprintf('Subject %d: %s', idxSubj, fh.Children(end-2).Title.String{2}));
        
        if iSubj == nSubjects
            try % only for R2019b and newer, but not critical to plot
                sgtitle(fh.Children(end-2).Title.String{1});
            end
            set(fh3, 'Name',fh.Children(end-2).Title.String{1})
        end
        
        if iSubj >= nSubjects -1
            xlabel('Scan Volumes');
        end
        
        close(fh);
    end
    qualityMeasuresArray(iSubj) = movement.quality_measures;
    
end
%% Output FD statistics
meanFD = [qualityMeasuresArray.meanFD];
fprintf('\n\tMean Framewise Displacement (FD) per subject:\n');
fprintf('\t  Subject %d: %4.3f mm\n ', [idxSubjArray; meanFD]);
fprintf('\n\tMean and standard deviation of mean Framewise Displacement (FD) over subjects:\n\t  %4.2f +/- %4.2f mm\n\n', ...
    mean(meanFD), std(meanFD));

%% Save plots
if doPlot
    fh = [fh2;fh3];
    if doSavePlot
        pathFigs = options.paths.figures;
        save_plots_abstract(fh, pathFigs, 3)
    end
else
    fh = [];
end
