% function output = create_subfigure_congruency_in_out(input)
% compute difference between spiral in and out and create respective subfigure
%
%   output = create_subfigure_congruency_in_out(input)
%
% IN
%
% OUT
%
% EXAMPLE
%   create_subfigure_congruency_in_out
%
%   See also
 
% Author:   Lars Kasper
% Created:  2019-07-31
% Copyright (C) 2019 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3. 
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
%

%% set paths
cd C:\Users\kasperla\Documents\Projects\SPIFI\Code\spiral_fmri_7t_paper\Code
spifi_setup_paths();

%% Load data
cd C:\Users\kasperla\Documents\Projects\SPIFI\Results\Pipeline02\SPIFI_0002\scandata
X = MrImage('meanafmri2_In.nii')
Y = MrImage('meanafmri2_Out.nii')

%% try out plotting
plot(Y-X,'z',14,'colorBar', 'on')
plot((abs(Y-X).^.5),'z',14,'colorBar', 'on', 'plotLabels', false, 'rotate90', 1); axis off; title('')