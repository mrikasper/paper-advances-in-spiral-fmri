function fh = create_figure_trajectory_gradients(options)
%Plots nominal gradient waveform and trajectory for spiral out and in/out,
%as used in the ReadInPatch for gradient execution
%
%   fh = create_figure_trajectory_gradients(options)
%
% IN
%
% OUT
%
% EXAMPLE
%   create_figure_trajectory_gradients
%
%   See also spifi_create_trajectories

% Author:   Lars Kasper
% Created:  2019-07-07
% Copyright (C) 2019 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich

if nargin < 1
    options = spifi_get_analysis_options();
    options = options.representation.fig_traj;
    options.doRemoveAxes = false;
end

%% derived parameters from options, e.g, file names of saved trajectories
paths = spifi_get_paths();

idSubject = 'SPIFI';
pathData = paths.data.trajectories;

pathCodeTrajCreation = fullfile(paths.code.analysis.root, 'Spirals');

fileNameIndex = fullfile(pathCodeTrajCreation, ...
    sprintf('index_gradient_files_%s.m', idSubject));

idArray = read_index_gradient_files(fileNameIndex);
iTrajToPlotArray = options.iTrajArray(options.iRunArray);

idsTrajToPlot = idArray(iTrajToPlotArray);
nTrajs = numel(idsTrajToPlot);

% load gradients
for iTraj = 1:nTrajs
    idTraj = idsTrajToPlot(iTraj);
    gwi = GradientWriterInput();
    gradobj{iTraj} = GradientWriter(GradientSystem(),gwi);
    gradobj{iTraj}.loadFromTxt(fullfile(pathData, idSubject, num2str(idTraj), ...
        'gradients.txt'));
end

%% plot gradients
for iTraj = 1:nTrajs
    stringTitle = sprintf('%s: Gradients', options.titles{options.iRunArray(iTraj)});
    fh(iTraj) = figure('Name', stringTitle);
    gradobj{iTraj}.plot('g')
    title(stringTitle);
    xLimits(iTraj,:) = xlim;
end


%% adjust plot for figure export
colors = [
    0   0   0
    91 155 213
    ];
colors = (colors+1)/256;

for iTraj = 1:nTrajs
    figure(fh(iTraj));
    xlim([min(xLimits(:,1)), max(xLimits(:,2))]);
    hLine = findobj(gca, 'Type','Line')
    set(hLine, 'LineWidth', 2)
    set(hLine, 'Color', colors(iTraj,:))
    set(hLine(2), 'LineStyle', '-.')
    if options.doRemoveAxes
        axis off;
        title ''
    end
end