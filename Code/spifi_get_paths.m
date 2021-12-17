function paths = spifi_get_paths(idxPipeline)
% returns environment-dependent (Mac/Win/Unix/Cluster) paths of SPIFI-project
% NOTE: expects ans environment variable 'BIOTOPE' for Linux to identify on
%       which cluster (EULER, IBT) you currently are.

if nargin < 1
    idxPipeline = 2;
end

fs = filesep;

%% #MOD# Change the following two paths to match your configuration.
% Note: All other paths are set automatically relative to those ones, or
% are irrelevant for the analysis/figure creation part
% Note2: For more flexibility, use the case/switch conditional below to set
% up username-specific environments

% raw fmri data, also behav/physlog, as transferred from FAIR repository
paths.data.root = 'Z:/kasperla_ibtnas03/SpiralFmri';

% root folder where results, intermediate data and figures will be written to
paths.results.root =  'Z:/kasperla_ibtnas01/SPIFI/Results';

% not needed for typical analysis, leave empty
paths.code.recon.root = '';

%% Start automatic configurations
username = get_username();

switch username
    case 'kasperla'
        
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
       
        %% #MOD UNCOMMENT THE FOLLOWING LINES AND ADD YOUR USERNAME AND FOLDERS
        
        %     case 'yourname'
        %         % raw fmri data, also behav/physlog, as transferred from FAIR repository
        %         paths.data.root = 'Z:/kasperla_ibtnas03/SpiralFmri';
        %
        %         % root folder where results, intermediate data and figures will be written to
        %         paths.results.root =  'Z:/kasperla_ibtnas01/SPIFI/Results';
        %
        %         % not needed for typical analysis, leave empty
        %         paths.code.recon.root = '';
    otherwise
        error('unknown user, set files in spifi_get_paths for your environment')
        
end

%% Generic paths, relative to this project

paths.code.analysis.root    = fileparts(mfilename('fullpath'));
paths.code.spm              = fullfile(paths.code.analysis.root, 'Toolboxes', 'spm12');
paths.code.spm_anatomy_toolbox = fullfile(paths.code.spm, 'toolbox', 'Anatomy');
paths.code.physio           = fullfile(paths.code.analysis.root, 'Toolboxes', 'TAPAS', 'PhysIO', 'code');
paths.code.uniqc            = fullfile(paths.code.analysis.root, 'Toolboxes', 'TAPAS', 'UniQC');
paths.code.recon.trajectories = fullfile(paths.code.recon.root, 'utils', 'nominalTrajectory');

paths.data.trajectories = fullfile(paths.data.root, 'Trajectories');
paths.data.spm_anatomy_toolbox.maps_early_visual = strcat(...
    paths.code.spm_anatomy_toolbox, ...
    fs, 'PMaps', fs, ...
    {
    'Visual_hOc1.nii'
    'Visual_hOc2.nii'
    'Visual_hOc3d.nii'
    'Visual_hOc3v.nii'
    });

paths.data.spm_anatomy_toolbox.mask_visual_cortex = fullfile...
    (paths.code.analysis.root, 'Revision2', 'ROI_VisualCortex_MNI.img');

%% Create dependent paths
% paths.analysis.summary = fullfile( paths.analysis.results, 'Summary-v0.2');
pipelineId =  sprintf('Pipeline%02d', idxPipeline);
paths.results.summary = fullfile( paths.results.root, 'Summary', pipelineId);
paths.results.root = fullfile(paths.results.root, pipelineId); % current pipeline

paths.code.recon.project =  fullfile(paths.code.analysis.root, ...
    '..', '..', 'recon-scripts-spifi');


%% paths for output figures
paths.figures = fullfile(paths.results.summary, 'Figures');