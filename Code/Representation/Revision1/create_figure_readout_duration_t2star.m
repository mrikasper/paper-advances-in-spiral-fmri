function fhArray = create_figure_readout_duration_t2star()
%Compares full spiral readout to cropped one (1mm, 37ms) to show effect of
%longer readout under T2*
%
%   fhArray = create_figure_readout_duration_t2star()
%
% Note: This needs to run SPIFI_0007_MainReconstructImages of
% recon-scripts-spifi (git@tnurepository.ethz.ch:laymm/recon-scripts-spifi.git)
% at commit 8a79c10886f27fc83419d720ddf936f9cf026a69 before
%
% IN
%
% OUT
%
% EXAMPLE
%   create_figure_readout_duration_t2star
%
%   See also
 
% Author:   Lars Kasper
% Created:  2020-06-18
% Copyright (C) 2020 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3. 
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
%

X = MrImage
X.read_data_from_graphics_handle(21)
Y = MrImage
Y.read_data_from_graphics_handle(14)
X.plot
Y.plot
plot(abs(X./max(X)-Y./max(Y)))