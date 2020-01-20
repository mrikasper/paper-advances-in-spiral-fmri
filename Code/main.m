% Script main
% Runs preprocessing and statistical analysis for all subjects and fMRI
% runs, for different image reconstructions extractions (e.g.,
% magnitude/phase, spiral in/out)
%
%  main
%
%
%   See also
 
% Author:   Lars Kasper
% Created:  2020-01-20
% Copyright (C) 2020 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3. 
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters (Modify to your needs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% idxSubj     subject ID (1-7)
% numeric index array referring to subjects SPIFI_0001 to SPIFI_0007
idxSubjectArray = 2:7

%   idxSess     session to be extracted (e.g. fMRIOut of fMRIInOut)
%               1   fMRIOut_1   (in paper high-res spiral out 0.8 mm)
%               2   fMRIInOut_1 (in papaer spiral in/out 1.5 mm)
%               3   fMRIOut_2
%               4   fMRIInOut_2
idxSessArray = 1:2;

%   idxRecon    reconstruction used (e.g. PartIn or PartOut of In/Out
%               spiral)
%               fmriOut
%               1 = magnitude
%               2 = phase
%               fmriInOut
%               1 = In          4 = In (phase)
%               2 = Out         5 = Out (phase)
%               3 = Combined    6 = Combined (phase)
%
% cell array referring to which recons are analyzed for each session (e.g.
%   magnitude for high-res spiral out,
%   magnitude of spiral-in, -out part and combined for spiral in/out
idxReconPerSessArray = {
    1
    1:3
    };

%   idxProcessingStepsArray
%               index vector, e.g., [1:12] of processing steps to be
%               executed
%                doDeletePreviousResults    -1 
%                doSetupPaths               0
%                doCopyRawData              1
%                doComputeBehav             2
%                doPreprocFunctional        3
%                doCombineEchoes = idxSess==2 && idxRecon == 3 && doPreprocFunctional;
%                doCoregAnatomical          4
%                doSegmentAnatomical        5
%                doBiasCorrectFunctional	6
%                doWarpImagesToMNI          7
%                doComputePhysio            8 % needs realignment params
%                doComputeGlm               9
%                doReportTaskContrasts      10
%                doReportPhysioContrasts 	11
%                doComputeSummary           12
%                doCreateFigures            13
idxProcessingStepsArray = 0:13;

%   options     See also spifi_get_analysis_options, subject-independent
%               settings, e.g. smoothing kernel, selected slices for
%               plotting
%
options = spifi_get_analysis_options();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run analysis in loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idxSubj = idxSubjectArray
    for idxSess = idxSessArray
        for idxRecon = idxReconPerSessArray{idxSess}
        spifi_run_analysis(idxSubj, idxSess, idxRecon, ...
            idxProcessingStepsArray, options);
        end
    end
end