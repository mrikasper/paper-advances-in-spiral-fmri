function fh = plot_rp(rp)
%Plots SPM realignment parameters
%
%    fh = plot_rp(rp)
%
% IN
%
% OUT
%
% EXAMPLE
%   plot_rp
%
%   See also
%
% Author:   Lars Kasper
% Created:  2016-12-19
% Copyright (C) 2016 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the Zurich fMRI Methods Evaluation Repository, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
%
% $Id: plot_rp.m 921 2016-12-19 14:30:03Z kasperla $
lgstr={'x (mm)', 'y (mm)', 'z (mm)','pitch (deg)', 'roll (deg)', 'yaw (deg)'};
titstr = sprintf('Realignment Parameters');
fh = figure('Name', titstr);
lineWidth = 3;
set(gcf,'DefaultAxesFontSize', 20);
set(gcf, 'DefaultTextFontSize', 20);
try
    load colors
catch
    colors = get(0, 'DefaultAxesColorOrder');
end

hs(1) = subplot(2,1,1);
plot(rp(:,1:3), 'LineWidth', lineWidth);
legend(lgstr{1:3});
title('Translation (mm)');

hs(2) = subplot(2,1,2);
plot(rp(:,4:6)*180/pi, 'LineWidth', lineWidth);
legend(lgstr{4:6});
title('Rotation (deg)')

linkaxes(hs, 'x');
suptitle(titstr);