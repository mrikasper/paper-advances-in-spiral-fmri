% Script coregister_boundaries_fsl_bbr
% Coregisters functional mean spiral (after bias-field correction) to mean 
% multi-echo (6-deg rigid body spm_coreg to rafmri1) using
% FSL FLIRT boundary-based registration
%
%  coregister_boundaries_fsl_bbr
%
% GOAL: To show that resulting deformation fields are small and boundaries
% are accurately represented by spirals
%
%   See also
 
% Author:   Saskia Bollmann & Lars Kasper
% Created:  2020-07-13
% Copyright (C) 2020 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3. 
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For SPIFI_0007:

pathData = '/Users/kasperla/Documents/Projects/SPIFI/Results/Pipeline02/SPIFI_0007/scandata';

cmdString = {...
  'mkdir fsl;' ...
  'cp cor2rafmri1_meanme1.nii fsl/t1.nii'...
  'cp mmeanafmri1.nii fsl/fmri.nii;'...
  'mv fsl/t1.nii struct.niil;'...
  'mv fsl/fmri.nii fsl/func.nii;'...
  'mv struct.nii fsl/'...
  };

% TODO: call to windows subsystem linux?
% system(cmdString);

%% FSL

% first (recommended) brain extraction on structural
% bet struct.nii struct_brain -m -Z

% FLIRT command line
% flirt -in func.nii -ref struct.nii -out bbr_func -cost bbr -omat affine_trafo.txt

% GUI command
% /usr/local/fsl/bin/flirt -in /mnt/c/Users/kasperla/Documents/Projects/SPIFI/Results/Pipeline02/SPIFI_0007/scandata/fsl/func.nii -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain -out /mnt/c/Users/kasperla/Documents/Projects/SPIFI/Results/Pipeline02/SPIFI_0007/scandata/fsl/bbr_func2.nii.gz -omat /mnt/c/Users/kasperla/Documents/Projects/SPIFI/Results/Pipeline02/SPIFI_0007/scandata/fsl/bbr_func2.mat -bins 256 -cost normmi -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp sinc -sincwidth 7 -sincwindow hanning



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UniQC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = MrImageSpm4D('struct.nii')
Y = MrImageSpm4D('func.nii')
Y.apply_threshold([0 0.8]).plot_overlays(X, 'overlayMode', 'edge')
Y.apply_threshold([0 0.8]).plot_overlays(X, 'overlayMode', 'edge')
X.edge('sobel',0.05).plot
Y.edge('sobel',0.06).plot

[xTPMs, ~, xBiasField] = X.segment
[yTPMs, ~, yBiasField] = Y.segment