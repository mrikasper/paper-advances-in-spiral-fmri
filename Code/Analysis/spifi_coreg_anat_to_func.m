function [meanX, coregX] = spifi_coreg_anat_to_func(coregOpts)
% Computes mean of all multi-echo echoes and coregisters it to the mean
% realigned functional image
%
%   [meanX, coregX] = spifi_coreg_anat_to_func(coregOpts)
%
% IN
%   coregOpts   details.preproc.coreg
% OUT
%
% EXAMPLE
%   spifi_coreg_anat_to_func
%
%   See also
 
% Author:   Lars Kasper
% Created:  2019-07-07
% Copyright (C) 2019 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich

X = MrImageSpm4D(coregOpts.source);
meanX = mean(X);
meanX.save('fileName', coregOpts.mean);

Y = MrImageSpm4D(coregOpts.ref);

% reslice to resolution of anatomical, because this will be used as a
% reference
Z = Y.mean('t').resize(meanX.dimInfo);

coregX = meanX.copyobj;

coregX.coregister_to(Z);
coregX.save('fileName', coregOpts.out);