function [onsets, durations, names, idxEvents] = ...
    get_multiple_conditions_visual(fileLog, varargin)
% Extracts paradigm logfile data and saves mat-file in SPM-compatible format
%
%   [onsets, durations, names] = ...
%       get_multiple_conditions_visual(fileLog, fileMultipleConditions, ...
%       timeUnit)
%
% IN
%   fileLog     cogent logfile created by run of main_paradigm_visual
%               e.g. 'Quad_localizer_FQSM_001_2015_01_14_12_54_48.log'
%
% Optional parameters as 'ParameterName', ParameterValue-pairs;
%
%   fileMultipleConditions
%               output mat-filename that variables shall be saved to
%               (if not given, data will not be saved)
%   timingSource  default: cycleStartLogfile
%               'cycleStartLogfile'         reads visual block timing from
%                                           cogent-logfile lines
%               'cycleStartMatVariable'     loads mat-variable with
%                                           cycle-start values directly
%                                           (created as cogent paradigm
%                                           output as well)
%               'scanTriggerCount'          creates block timing from
%                                           counted scan triggers and scan
%                                           timing information in mat-file
%                                           (variable 'inputArguments')
%
%   timeUnit    'scans' or 'seconds' (default)
%               determines whether onsets/durations-values are given in
%               seconds, or in fractions of scans
%               (where start of 1st scan = 0, start of 2nd scan = 1 etc.)
%
% OUT
%   onsets      cell(1, nConditions) of stimulus onset vectors
%   durations   cell(1, nConditions) of stimulus duration vectors
%   names       cell(1, nConditions) of condition name strings
%   idxEvents   row indices (ignoring header lines) of occuring events, struct containing:
%               idxButtonLeft
%               idxButtonRight
%               timeEventMilliSeconds
%               idxScanTrigger
%               idxCycleStart
%
% EXAMPLE
%   [onsets, durations, names] = get_multiple_conditions_visual(...
%       'logs/visual_FQSM_001_2015_01_14_12_54_48.log', ...
%       'logs/'multiple_conditions.mat', ...
%       'seconds');
%
%   See also main_paradigm_visual
%
% Author:   Lars Kasper
% Created:  2015-01-21
% Copyright (C) 2015 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the Zurich fMRI Methods Evaluation Repository, which is released
% under the terms of the GNU General Public Licence (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.
%
% $Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set and update input parameters
if nargin < 1
    fileLog = fullfile('logs','visual_FQSM_007_2015_12_18_13_54_54.log');
end

defaults.isVerbose              = true;
defaults.fileMultipleConditions = 'multiple_conditions.mat';
defaults.timingSource           = 'cycleStartLogfile';
defaults.timeUnit               = 'seconds';

inputArguments = propval(varargin, defaults);
strip_fields(inputArguments);

fileMat = regexprep(fileLog, '\.log', '\.mat');
hasMatInfoFile = exist(fileMat, 'file');

% use timing information from mat-file to compute stimulus onsets from scan
% onsets
if hasMatInfoFile
    load(fileMat)
end

if ~exist('inputArguments', 'var')
    % create default set of input parameters, if none were saved
    % NOT 5, since we want asynchrony to TR to get different onset delays and sub-sample the HRF
    inputArguments.nTrPerBlock                    = 4.975;
    inputArguments.trSeconds                      = 3;
    
    % number of cycles all block are repeated
    inputArguments.nDummies                       = 0;
    inputArguments.nCycles                        = 8;
    inputArguments.codeSerialButtonLeft           = 52;
    inputArguments.codeSerialButtonRight          = 51;%49;
    inputArguments.codeSerialTrigger              = 53;
    inputArguments.doIncludeFixationBlocks        = true;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read log file, generated by cogent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% example lines from log-file:
%
% 38402	[38401]	:	Byte	53	COM1	at	38401       # scan trigger
% 107504	[69102]	:	Byte	49	COM1	at	39008   # button 1 pressed
% 107504	[0]	:	Cycle	1	Start	at	41401       # stim cycle start
% 107505	[1]	:	Byte	53	COM1	at	41402       # button 2 pressed

fid = fopen(fileLog);

%Remove lines until COGENT START found
tline = [];
while isempty(strfind(tline, 'COGENT START'))
    tline = fgetl(fid);
end

% ignore error lines, starting with ERR
columnsLog = ...
    textscan(fid,'%d [%d] : %s %d %s at %d', 'CommentStyle', 'ERR');

fclose(fid);

timeLogMilliSeconds                 = columnsLog{1};

durationSinceLastLogMilliSeconds    = columnsLog{2};
stringEvent1                        = columnsLog{3};
codeEvent                           = columnsLog{4};
stringEvent2                        = columnsLog{5};
timeEventMilliSeconds               = double(columnsLog{6});

% remove last row with COGENT STOP...
idxRowRemove = find_string(stringEvent1, 'COGENT');
timeLogMilliSeconds(idxRowRemove)                 = [];
durationSinceLastLogMilliSeconds(idxRowRemove)    = [];
stringEvent1(idxRowRemove)                        = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set names of conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = { ...
    'UpperLeftLowerRight (ULLR)', ...
    'UpperRightLowerLeft (URLL)', ...
    'Button Presses Left', ...
    'Button Presses Right'
    };

nConditions = numel(names);
onsets      = cell(1, nConditions);
durations   = cell(1, nConditions);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set durations of conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch timingSource
    % fixed block lengths for stimulations
    case {'scanTriggerCount', 'cycleStartLogfile'}
        durationBlock =  ...
            inputArguments.nTrPerBlock * inputArguments.trSeconds;
        durations(1:2) = {durationBlock};
end

% events have zero duration
durations(3:4) ={0};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract onsets of conditions from event-codes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idxButtonLeft = find(codeEvent == inputArguments.codeSerialButtonLeft);
idxButtonRight = find(codeEvent == inputArguments.codeSerialButtonRight);
idxScanTrigger = find(codeEvent == inputArguments.codeSerialTrigger);
idxCycleStart = find_string(stringEvent1, 'Cycle');

idxEvents = struct('idxButtonLeft', idxButtonLeft, ...
    'idxButtonRight',idxButtonRight , 'idxScanTrigger', idxScanTrigger, ...
    'idxCycleStart', idxCycleStart, 'timeEventMilliSeconds', timeEventMilliSeconds);

% change to trigger count, if no cycle start logged
if isempty(idxCycleStart)
    timingSource = 'scanTriggerCount';
end

hasScanTriggers = ~isempty(idxScanTrigger);

if hasScanTriggers
    time1stScanMilliSeconds = timeEventMilliSeconds(idxScanTrigger(1));
else
    time1stScanMilliSeconds = timeStartParadigmSeconds*1000;
end


% block count different, if fixation included
nBlocksPerCycle = 2;
if isfield(inputArguments, 'doIncludeFixationBlocks') && ...
        inputArguments.doIncludeFixationBlocks
    nBlocksPerCycle = nBlocksPerCycle*2;
end
nTrPerCycle = nBlocksPerCycle*inputArguments.nTrPerBlock;

switch timingSource
    case 'cycleStartLogfile'
        onsets{1} = timeEventMilliSeconds(idxCycleStart);
        onsets{2} = onsets{1} + 1000*inputArguments.trSeconds*nTrPerCycle/2;
        
    case 'scanTriggerCount'
        
        if ~isfield(inputArguments, 'nTrFixationStart')
            % deprecated, old version
            idxScansOnsetCondition{1} = inputArguments.nDummies + ...
                (1:nTrPerCycle:(inputArguments.nCycles*nTrPerCycle))';
        else
            idxScansOnsetCondition{1} = inputArguments.nTrFixationStart + ...
                (1:nTrPerCycle:(inputArguments.nCycles*nTrPerCycle))';
        end
        
        
        idxScansOnsetCondition{2} = idxScansOnsetCondition{1} + ...
            nTrPerCycle/2;
        
        % Computation of onset times, also for non-integer block multiple
        % via last trigger plus fraction of duration to next trigger
        deltaTimeScanTrigger = diff(timeEventMilliSeconds(idxScanTrigger));
        for c = 1:2
            lastScanOnsets = floor(idxScansOnsetCondition{c});
            onsets{c} = timeEventMilliSeconds(...
                idxScanTrigger(lastScanOnsets)) + ...
                mod(idxScansOnsetCondition{c},1).* ...
                deltaTimeScanTrigger(lastScanOnsets);
        end
        
end

onsets{3} = timeEventMilliSeconds(idxButtonLeft);
onsets{4} = timeEventMilliSeconds(idxButtonRight);

% time relative to 1st scan, convert to seconds
onsets = cellfun(@(x) (x - time1stScanMilliSeconds)/1000, ...
    onsets, 'UniformOutput', false);



if isVerbose
    plot_multiple_conditions(onsets, durations, names);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove empty conditions to avoid error in SPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idxEmptyCondition = find(cellfun(@isempty, onsets));

names(idxEmptyCondition) = [];
onsets(idxEmptyCondition) = [];
durations(idxEmptyCondition) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save file multiple_conditions-file, if specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(fileMultipleConditions)
    save(fileMultipleConditions, 'names', 'onsets', 'durations');
end