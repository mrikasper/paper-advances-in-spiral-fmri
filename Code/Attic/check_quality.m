function qcImages = check_quality(fileNameNifti)
% Computes Quality Control measures for nifti time series and plots the
% relevant ones
%
%  qcImages  = check_quality(fileNameNifti)
%
% IN
%
% OUT
%   qcImages    cell(nImages,1) of MrImage, relevant quality control
%               images, e.g., mean/tsnr/maxdiffabs etc.
%
% EXAMPLE
%   check_quality
%
%   See also
%
% Author:   Lars Kasper
% Created:  2018-05-23
% Copyright (C) 2018 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3. 
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
%
% $Id: new_function2.m 354 2013-12-02 22:21:41Z kasperla $

Y = MrImage(fileNameNifti);

Y.mean.plot()
Y.std.plot()
Y.snr.plot()

plot(abs(Y-mean(Y)), 'useSlider', 'true');

plot(image2k(Y-mean(Y)), 'useSlider', 'true');



qcImages{1} = Y;