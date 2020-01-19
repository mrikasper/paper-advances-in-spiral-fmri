function split_map = split_brewermap(N, scheme, posneg, funcColorMap)
%splits brewer map to be only used for positive or negative parts of the map
% takes either upper or lower half of brewer map, but turns around lower
% part
%
%  split_map = split_brewermap(N, scheme, posneg)
%
% IN
%   scheme  either a scheme from brewermap or a color map handle
%   posneg  'pos' or 'neg' to get upper half or flipped lower part of
%           brewermap
%   funcColorMap
%           default: brewermap, could be other colormaps to split as well,
%           e.g., jet, hsv
% OUT
%
% EXAMPLE
%   split_brewermap
%
%   See also brewermap
 
% Author:   Lars Kasper
% Created:  2019-08-31
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

switch scheme
    case {'jet', 'parula'}
        funcColorMap = str2func(scheme);
        split_map = funcColorMap(2*N);
    otherwise
        split_map = brewermap(2*N, scheme);
end

switch posneg
    case 'pos'
        split_map = split_map((N+1):end,:);
    case 'neg'
        split_map = flip(split_map(1:N,:));
end