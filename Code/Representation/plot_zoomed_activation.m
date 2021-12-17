% Script plot_zoomed_activation
% Trying out position (not index-)based zooming, on differential vs block
% regressors
%
%  plot_zoomed_activation
%
%
%   See also
 
% Author:   Lars Kasper
% Created:  2019-09-05
% Copyright (C) 2019 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3. 
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd C:\Users\kasperla\Documents\Projects\SPIFI\Results\Pipeline02\SPIFI_0007\scandata
X =MrImage('mcor2rafmri2_Combined_meanme1.nii')
cd C:\Users\kasperla\Documents\Projects\SPIFI\Results\Pipeline02\SPIFI_0007\glm\s08_fmri2_Combined
Y  = MrImage('spmT_0001.nii')
Z1 = MrImage('spmT_0013.nii')
Z2 = MrImage('spmT_0014.nii')

%% Plot overview
X.threshold([0.1 1.4]).plot_overlays({Y,Y.*-1}, 'selectedSlices',18:21, 'overlayThreshold', [3 8], 'overlayMode', 'map', 'colorBar', 'off', 'overlayColorMaps', {'autumn','winter'}, 'rotate90',1,'nRows',2)
X.threshold([0.1 1.4]).plot_overlays({Z1,Z2}, 'selectedSlices',18:21, 'overlayThreshold', [3 8], 'overlayMode', 'map', 'colorBar', 'off', 'overlayColorMaps', {'autumn','winter'}, 'rotate90',1,'nRows',2)

%% Create and apply zoom
cropSamplesY = X.dimInfo.index2sample([10; 80],2);
zX = X.select('type', 'samples', 'y', cropSamplesY(1):cropSamplesY(2))
zY = Y.select('type', 'samples', 'y', cropSamplesY(1):cropSamplesY(2))
zZ1 = Z1.select('type', 'samples', 'y', cropSamplesY(1):cropSamplesY(2))
zZ2 = Z2.select('type', 'samples', 'y', cropSamplesY(1):cropSamplesY(2))
zX.threshold([0.1 1.4]).plot_overlays({zY,zZ.*-1}, 'selectedSlices',1:36, 'overlayThreshold', [3 8], 'overlayMode', 'map', 'colorBar', 'off', 'overlayColorMaps', {'autumn','winter'}, 'rotate90',1, 'nRows',9)
zX.threshold([0.1 1.4]).plot_overlays({zY,zY.*-1}, 'selectedSlices',1:36, 'overlayThreshold', [3 8], 'overlayMode', 'map', 'colorBar', 'off', 'overlayColorMaps', {'autumn','winter'}, 'rotate90',1, 'nRows',9)
