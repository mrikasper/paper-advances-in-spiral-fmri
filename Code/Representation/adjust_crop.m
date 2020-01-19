function [cropX2, cropY2] = adjust_crop(cropX1, cropY1, Y1, Y2)
%rescales plot according to dimInfo scaling
%
%  [cropX2, cropY2] = adjust_crop(cropX1, cropY1, Y1, Y2);
%
% IN
%   Y1,2    MrImage
%   cropX1,Y1 [1,2] vector of crop pixel index start/end
%
% OUT
%
% EXAMPLE
%   adjust_crop
%
%   See also

% Author:   Lars Kasper
% Created:  2019-06-30
% Copyright (C) 2019 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich

cropX2 = round(cropX1/Y1.dimInfo.x.nSamples...
    .*Y2.dimInfo.x.nSamples);
cropY2 = round(cropY1/Y1.dimInfo.y.nSamples...
    .*Y2.dimInfo.y.nSamples);