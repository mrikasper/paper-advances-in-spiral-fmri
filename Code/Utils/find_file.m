function foundFile = find_file(searchPattern, isFolder)
% finds first existing file matching a search pattern (windows-style, with
% *)
%
%   foundFile = find_file(searchPattern)
%
% IN
%
% OUT
%
% EXAMPLE
%   find_file
%
%   See also
%
% Author:   Lars Kasper
% Created:  2018-09-20
% Copyright (C) 2018 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich

% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
if nargin < 2
    isFolder = false;
end

x = dir(regexprep(searchPattern,'**','*'));

% only take folder/file, if specified
x = x([x(:).isdir] == isFolder);

if numel(x) > 0
    foundFile = fullfile(x(1).folder, x(1).name);
else
    error('File not found %s', searchPattern);
end