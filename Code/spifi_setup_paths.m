function paths = spifi_setup_paths(mode)
% setup paths for SPIFI fMRI project (Analysis and recon), set pwd and some matlab environment

if nargin < 1 
	% determines certain defaults for matlab plotting
	mode = 'spm'; % 'recon', 'qc','analysis' or 'spm'; 
end

fprintf('\n\tSetting up paths for SPIFI, mode: %s...\n', mode);
paths = spifi_get_paths();

%% Add code paths
fprintf('\tRestoring default paths, adding new paths for SPIFI...');

restoredefaultpath;

addpath(paths.code.spm);
addpath(genpath(paths.code.analysis.root));
addpath(genpath(paths.code.physio));
addpath(genpath(paths.code.uniqc));
addpath(genpath(paths.code.recon.root));
addpath(genpath(paths.code.recon.project));

fprintf('done.\n');

%% Set pwd and plotting options
switch mode
	case 'spm'
		set(0, 'DefaultFigureWindowStyle', 'normal');
        fprintf('\tSetting up SPM paths...')
        spifi_setup_spm();
        fprintf('done.\n');
		cd(paths.code.analysis.root);
	case {'analysis', 'qc'}
		set(0, 'DefaultFigureWindowStyle', 'docked');
		cd(paths.code.analysis.root);
case 'recon'
		cd(paths.code.recon.project);
		set(0, 'DefaultFigureWindowStyle', 'docked');
end

fprintf('\n\tFinished setup for SPIFI. Have fun!\n\n\n');