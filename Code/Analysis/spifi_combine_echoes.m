function combinedX = spifi_combine_echoes(combineOpts)
% Computes mean of all multi-echo echoes and coregisters it to the mean
% realigned functional image
%
%   combinedX = spifi_combine_echoes(combineOpts)
%
% IN
%   combineOpts   details.preproc.combine_echoes
% OUT
%
% EXAMPLE
%   spifi_coreg_anat_to_func
%
%   See also
 
% Author:   Lars Kasper
% Created:  2019-08-04
% Copyright (C) 2019 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich

X = MrImageSpm4D(combineOpts.source{1});
Y = MrImageSpm4D(combineOpts.source{2});
meanX = mean(X);
meanY = mean(Y);

combinedX = (X.*meanX + Y.*meanY)./(meanX+meanY);
combinedX.save('fileName', combineOpts.out);