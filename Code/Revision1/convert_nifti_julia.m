function fileOut = convert_nifti_julia(fileIn, fileOut)
%Uses Matlab's niftiwrite to convert UniQC niftis to read with Nifti.jl
%package (to avoid "mappedarray not defined" error
%
%   fileOut = convert_nifti_julia(fileIn, fileOut)
%
% IN
%
% OUT
%
% EXAMPLE
%   convert_nifti_julia
%
%   See also
 
% Author:   Lars Kasper
% Created:  2021-06-29
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

if nargin < 1
    % default set for exported nifti maps
    pathIn = 'D:\SPIFI\ExportETHResearchCollection\SPIFI_0007_MapsForReconstruction';
    pathOut = fullfile(pathIn, 'ConversionNiftiWriteForJulia');
    fileIn = {
        'b0MapRaw1_Hz.nii'
        'b0MapSmoothedResizedSpiralGeometry1_Hz.nii'
        'coilSensitivityMaps_phase.nii'
        'b0MapSmoothed1_Hz.nii'
        'coilSensitivityMaps_magnitude.nii'
        };
    fileOut = strcat(pathOut, filesep, fileIn);
    fileIn = strcat(pathIn, filesep, fileIn);
end


if ~iscell(fileIn)
    fileIn = {fileIn};
end

if ~iscell(fileOut)
    fileOut = {fileOut};
end

for iFile = 1:numel(fileIn)
    % Load a NIfTI image using its .nii file name.
    info = niftiinfo(fileIn{iFile});
    V = niftiread(info);
  
    % Edit the description of the file.
    info.Description = 'Modified using MATLAB R2019b niftiwrite';
  
    % Write the image to a .nii file.
    niftiwrite(V, fileOut{iFile}, info);
end
