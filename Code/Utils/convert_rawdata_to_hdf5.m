%function output = convert_rawdata_to_hdf5(input)
% Converts MR Raw data (coil, traj, sense/b0 maps) to hdf5, to be used by
% MRIReco.jl and RRSG Challenge reference implementations
%
%   output = convert_rawdata_to_hdf5(input)
%
% IN
%
% OUT
%
% EXAMPLE
%   convert_rawdata_to_hdf5
%
%   See also

% Author:   Lars Kasper, based on code by Franz Patzig, IBT
% Created:  2020-06-09
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

saveCase = 'singleSlice_spiral0p8mm'; %'testSpiral'; 'singleSlice_spiral0p8mm'

pathSave = 'C:/Users/kasperla/Documents/Code/BRAIN-To/mrirecon-julia-spirals';

switch saveCase
    case 'testSpiral'
        fileMat = fullfile('C:\Users\kasperla\Documents\Projects\SPIFI\Results', ...
            'SPIFI_0007_ExportDataForPublication_testSpiral_testSpiral4il1mmSENSE_m400b424m400s1.mat');
        nInterleaves = 4;
    otherwise
        fileMat = fullfile('C:\Users\kasperla\Documents\Projects\SPIFI\Results', ...
            'SPIFI_0007_ExportDataForPublication_singleSlice_spiral0p8mm_m400b424m400s1.mat');
        nInterleaves = 1;
end
load(fileMat, 'trajectory', 'rawdata', 'senseMaps', 'b0Map_rad_per_s', 'b0Shift_rad_per_s');


nChannels = 32;
nDims = 3;

signal = rawdata;
traj = trajectory.k_xyz(:,1:nDims);
sense.data = senseMaps;

k_tmp = reshape(traj, length(traj)/nInterleaves, nInterleaves, nDims);
data = reshape(signal, length(traj)/nInterleaves, nInterleaves, nChannels);

% dimension permutation to initial RRSG challenge format
% rawdata: ch, il, read was order
% traj: il, read, kdim(x,y,z)
% Julia order: read, il, ch and kdim, read, il
doPermuteToInitialRRSG = true;
if doPermuteToInitialRRSG
    data = permute(data, [3, 2, 1]);
    k_tmp = permute(k_tmp, [2 1 3]);
end

sizeRawData = size(data);
sizeRawData(end+1:4) = 1; % make 4-dim

delete('data.h5');
h5create('data.h5','/trajectory', size(k_tmp))
h5write('data.h5', '/trajectory', k_tmp)
h5create('data.h5','/rawdata_r', sizeRawData)
h5write('data.h5', '/rawdata_r', real(data))
h5create('data.h5','/rawdata_i', sizeRawData)
h5write('data.h5', '/rawdata_i', imag(data))


%% Create SENSE Map h5
delete('sense.h5');
h5create('sense.h5','/senseMaps_r', size(sense.data))
h5write('sense.h5', '/senseMaps_r', real(sense.data))
h5create('sense.h5','/senseMaps_i', size(sense.data))
h5write('sense.h5', '/senseMaps_i', imag(sense.data))


%% Create B0 map h5
delete('b0.h5')
h5create('b0.h5','/b0Map', size(b0Map_rad_per_s))
h5write('b0.h5', '/b0Map', b0Map_rad_per_s+b0Shift_rad_per_s)

mkdir(fullfile(pathSave, saveCase))
movefile('*.h5', fullfile(pathSave, saveCase), 'f')