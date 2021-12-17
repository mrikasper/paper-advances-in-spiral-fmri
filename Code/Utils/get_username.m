function username = get_username()
% returns username on Linux, Mac or Windows
%
%   username = get_username();
%
% IN
%
% OUT
%
% EXAMPLE
%   get_username
%
%   See also
 
% Author:   Saskia Bollmann & Lars Kasper
% Created:  2021-12-16
% Copyright (C) 2021 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3. 
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
%

if ispc
   username = getenv('username'); 
else
   [~, username] = system('whoami');
   username(end) = [];
end