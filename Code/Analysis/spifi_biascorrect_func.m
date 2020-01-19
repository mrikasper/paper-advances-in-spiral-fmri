function [meanX, biasCorrectedX] = spifi_biascorrect_func(biascorrectOpts)
% Biascorrects mean functional image by applying previously calculated
% biasfield of segmentation (of anatomical image)
%
%   [meanX, coregX] = spifi_biascorrect_func(biascorrectOpts)
%
% IN
%   biascorrectOpts   details.preproc.coreg
% OUT
%
% EXAMPLE
%   spifi_biascorrect_func
%
%   See also
 
% Author:   Lars Kasper
% Created:  2019-07-07
% Copyright (C) 2019 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich

X = MrImage(biascorrectOpts.source);
Y = MrImage(biascorrectOpts.biasfield);

% it's a multiplication, not division somehow, so rather the inverse bias
% field...resampling is taken care of by UniQC
biasCorrectedX = X.*Y;

biasCorrectedX.save('fileName', biascorrectOpts.out);