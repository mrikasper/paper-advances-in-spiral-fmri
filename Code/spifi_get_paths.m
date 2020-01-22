function paths = spifi_get_paths(idxPipeline)
% returns environment-dependent (Mac/Win/Unix/Cluster) paths of SPIFI-project
% NOTE: expects ans environment variable 'BIOTOPE' for Linux to identify on
%       which cluster (EULER, IBT) you currently are.

if nargin < 1
    idxPipeline = 2;
end


%% #MOD# Change the following two paths to match your configuration.
% Note: All other paths are set automatically relative to those ones, or
% are irrelevant for the analysis/figure creation part

% raw fmri data, also behav/physlog, as transferred from FAIR repository
paths.data.root = 'Z:/kasperla_ibtnas03/SpiralFmri';      

% root folder where results, intermediate data and figures will be written to
paths.results.root =  'Z:/kasperla_ibtnas01/SPIFI/Results'; 


%% Start automatic configurations

if ispc % LarsPC
        projectRoot = 'C:\Users\kasperla\Documents\Projects\SPIFI\';
    paths.data.root = fullfile(projectRoot, 'Data'); % raw coil/monitoring data, also behav/physlog, as transfered from scanner
    paths.results.root = fullfile(projectRoot, 'Results');
    paths.code.recon.root = 'C:\Users\kasperla\Documents\Code\Recon';
    
elseif ismac
    
    paths.data.root = ''; % raw coil/monitoring data, also behav/physlog, as transfered from scanner
    paths.results.root =  '';
    paths.code.recon.root       = '/home/kasperla/Documents/code/matlab/recon_svn/trunk';
    
else % UNIX
    
    switch upper(getenv('BIOTOPE'))
        case 'EULER'
            paths.home = '/cluster/home/kasperla';
        case 'IBT'
            paths.home = '/home/kasperla';
    end
    
    paths.data.root = fullfile(paths.home, 'TransferSPIFI'); % raw coil/monitoring data, also behav/physlog, as transfered from scanner
    paths.results.root = fullfile(paths.home, 'ResultsSPIFI');
    paths.code.recon.root = fullfile(paths.home, 'code', 'matlab', 'recon_svn_trunk');
    
end % system-dependent paths


%% Generic paths, relative to this project

paths.code.analysis.root    = fileparts(mfilename('fullpath'));
paths.code.tnufmri          = fullfile(paths.code.analysis.root, 'Toolboxes', 'tnufmri');
paths.code.spm              = fullfile(paths.code.analysis.root, 'Toolboxes', 'spm12');
paths.code.physio           = fullfile(paths.code.analysis.root, 'Toolboxes', 'PhysIO', 'code');
paths.code.uniqc            = fullfile(paths.code.analysis.root, 'Toolboxes', 'UniQC');
paths.code.uniqc_tasks      = fullfile(paths.code.analysis.root, 'Toolboxes', 'UniQC-Tasks');
paths.code.recon.trajectories = fullfile(paths.code.recon.root, 'utils', 'nominalTrajectory');

paths.data.trajectories = fullfile(paths.data.root, 'Trajectories');

%% Create dependent paths
% paths.analysis.summary = fullfile( paths.analysis.results, 'Summary-v0.2');
pipelineId =  sprintf('Pipeline%02d', idxPipeline);
paths.results.summary = fullfile( paths.results.root, 'Summary', pipelineId);
paths.results.root = fullfile(paths.results.root, pipelineId); % current pipeline

paths.code.recon.project =  fullfile(paths.code.analysis.root, ...
    '..', '..', 'recon-scripts-spifi');


%% paths for output figures
paths.figures = fullfile(paths.results.summary, 'Figures');