%% Generating a spiral ISMRMRD data set
% Based on test_create_undersampled_dataset from ISMRMRD/examples/matlab
%
% Data from 32 coils from 36 slice object with 100 repetitions.
%

clear all;

%% input parameters
pathRecon = 'D:\SPIFI\Results';%C:\Users\kasperla\Documents\Projects\SPIFI\Results';
pathOutput = 'C:\Users\kasperla\Documents\Projects\SPIFI\Results\ISMRMRD';
saveCase = 'startMiddleEndVolume';
doDeleteExistingISMRMRDFile = true;

%% derived parameters
% for geometry of dataset, use test data:
switch(saveCase)
    case 'testSpiral'
        rwIn.dataId = 'testSpiral4il1mmSENSE_m400b424m400s1';
    otherwise
        rwIn.dataId = 'spiral0p8mm_m400b424m400s1';
end

dataType = '';
interleaves = 1;
R = 4; % accelerationFactor

switch saveCase
    case 'testSpiral'
        sli = 18;
        dyn = 2;
        interleaves = 1:4;
        R = 1;
    case 'singleSlice'
        sli = 18;
        dyn = 1;
    case 'singleVolume'
        sli = 1:36;
        dyn = 1;
    case 'startMiddleEndVolume'
        sli = 1:36;
        dyn = [1 50 100];
        dataType = '_data';
    case 'completeRun'
        sli = 1:36;
        dyn = 1:100;
end

% data from recon6, right before iterativeReconstruction
fileRecon = sprintf('SPIFI_0007_ExportDataForPublication_%s_%s%s.mat', ...
    saveCase, rwIn.dataId, dataType);

fileGeom = fullfile(pathRecon, ...
    'SPIFI_0007_ExportDataForPublication_testSpiral_testSpiral4il1mmSENSE_m400b424m400s1.mat');


% Output file Name
filename = fullfile(pathOutput, ...
    regexprep(fileRecon, '\.mat', '\.h5'));
if doDeleteExistingISMRMRDFile
    delete(filename);
end

fileRecon = fullfile(pathRecon, fileRecon);

%%
load(fileGeom, 'geometry');
load(fileRecon, 'rawdata', 'trajectory');

%%
dset = ismrmrd.Dataset(filename);

nInterleaves = numel(interleaves);
nSlices = numel(sli);
nReps = numel(dyn);
nCoils = 32;
nSamples = size(rawdata,1);

%%
% It is very slow to append one acquisition at a time, so we're going
% to append a block of acquisitions at a time.
% In this case, we'll do it one repetition at a time to show off this
% feature.  Each block has nSlices*nInterleaves aquisitions

acqblock = ismrmrd.Acquisition(nSlices*nInterleaves);

% Set the header elements that don't change
acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = nSamples;
acqblock.head.center_sample(:) = 0; % center of k-space
acqblock.head.active_channels(:) = nCoils;
acqblock.head.trajectory_dimensions(:) = 3; % rotated mps
acqblock.head.sample_time_us(:) = 1.8;


% measurement, phase and slice dir in x,y,z coordinates
rotationPhilipsXYZToPatient = geometry.rot_traj_xyz_to_mps';

acqblock.head.read_dir  = repmat(rotationPhilipsXYZToPatient(:,1),[1 nSlices*nInterleaves]);
acqblock.head.phase_dir = repmat(rotationPhilipsXYZToPatient(:,2),[1 nSlices*nInterleaves]);
acqblock.head.slice_dir = repmat(rotationPhilipsXYZToPatient(:,3),[1 nSlices*nInterleaves]);


% Loop over the acquisitions, set the header, set the data and append
for rep = 1:nReps
    fprintf('Appending repetition %d/%d\n', rep, nReps);
    for iInterleaf = 1:nInterleaves
        
        
        for iSlice = 1:nSlices
            
            
            % acquisition number w/i acquisition block made out of
            % interleaves and slices
            acqno = (iInterleaf-1)*nSlices + iSlice; % matlab counting!
            
            % Set the header elements that change from acquisition to the next
            % c-style counting
            scan_counter = (rep-1)*nSlices*nInterleaves + acqno - 1;
            
            acqblock.head.scan_counter(acqno) = scan_counter;
            % Note next entry is k-space encoded line number (not acqno which
            % is just the sequential acquisition number)
            acqblock.head.idx.slice(acqno) = sli(iSlice) - 1;
            acqblock.head.idx.kspace_encode_step_1(acqno) = interleaves(iInterleaf) - 1;
            acqblock.head.idx.repetition(acqno) = dyn(rep) - 1;
            
            % Set the flags
            acqblock.head.flagClearAll(iSlice);
            if iSlice == 1
                %            acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', acqno);
                %            acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', acqno);
                acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', acqno);
            elseif iSlice==nSlices
                %            acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);
                %            acqblock.head.flagSet('ACQ_LAST_IN_SLICE', acqno);
                acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', acqno);
            end
            
            
            % fill the data, respect convention:
            %          size(acqRead.traj{1})
            %               ans =
            %                   3        1284
            %               size(acqRead.data{1})
            %               ans =
            %                   1284          36
            acqblock.data{iSlice} = rawdata(:,:, scan_counter+1);
            acqblock.traj{iSlice} = (trajectory.k_xyz(:,:,scan_counter+1)).';
        end % interleaf loop
        
    end % slice loop
    
    if rep == nReps % last acquisition in measurement
        acqblock.head.flagSet('ACQ_LAST_IN_MEASUREMENT', iSlice);
    end
    
    % Append the acquisition block
    dset.appendAcquisition(acqblock);
    
end % rep loop

% TODO: Last in Measurement


%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill the xml header %
%%%%%%%%%%%%%%%%%%%%%%%%
% We create a matlab struct and then serialize it to xml.
% Look at the xml schema to see what the field names should be

header = [];

% subject Information
header.subjectInformation.patientName = 'SPIFI_0007';
header.subjectInformation.patientWeight_kg = '70';
header.subjectInformation.patientID = 'SPIFI_0007';
header.subjectInformation.patientBirthdate = '19830101';
header.subjectInformation.patientGender = 'M';

% study info
header.studyInformation.studyDate='20171124';
header.studyInformation.studyID='SPIFI';
header.studyInformation.studyDescription='SPIral Functional Imaging';

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 298023839; % 7T, subject-specific!

% Measurement Information
header.measurementInformation.measurementID = 'SPIFI_0007';
header.measurementInformation.patientPosition = 'HFS';
header.measurementInformation.seriesDate = '20171124'; % pseudonymized, only year correct
header.measurementInformation.protocolName = 'SPIFI_0007_2017_11_24_Subject7';
header.measurementInformation.seriesDescription = 'Spiral_0p8mm_R4_TE20_Run1';

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.systemVendor = 'Philips';
header.acquisitionSystemInformation.systemModel = 'Achieva';
header.acquisitionSystemInformation.institutionName = 'Insitute for Biomedical Engineering, ETH Zurich and University of Zurich, University Hospital Zurich';
header.acquisitionSystemInformation.receiverChannels = nCoils;
header.acquisitionSystemInformation.systemFieldStrength_T = 7;
header.acquisitionSystemInformation.stationName = 'Radiology';

% Sequence parameters
header.sequenceParameters.TR = 3.306;
header.sequenceParameters.TE = 20;
header.sequenceParameters.flipAngle_deg = 90;
header.sequenceParameters.sequence_type = 'GRE';

% The Encoding (Required)
header.encoding.trajectory = 'spiral';
% use intended FOV here, see https://gadgetron.discourse.group/t/ismrmrd-file-generation/141/5
header.encoding.encodedSpace.fieldOfView_mm.x = 190;
header.encoding.encodedSpace.fieldOfView_mm.y = 230;
header.encoding.encodedSpace.fieldOfView_mm.z = 36;
% taken from Recon7 converter:
header.encoding.encodedSpace.matrixSize.x = nSamples;% 288/R;
header.encoding.encodedSpace.matrixSize.y = 1; %288/R;
header.encoding.encodedSpace.matrixSize.z = 36; %36/R;

% Recon Space
% (in this case same as encoding space)
header.encoding.reconSpace.fieldOfView_mm.x = 190;
header.encoding.reconSpace.fieldOfView_mm.y = 230;
header.encoding.reconSpace.fieldOfView_mm.z = 36;
header.encoding.reconSpace.matrixSize.x = 240;
header.encoding.reconSpace.matrixSize.y = 292;
header.encoding.reconSpace.matrixSize.z = 36;

% Encoding Limits
% phase encoding direction, 1 il, see https://gadgetron.discourse.group/t/ismrmrd-file-generation/141/5
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = 0;
header.encoding.encodingLimits.kspace_encoding_step_1.center = 0;
header.encoding.encodingLimits.slice.minimum = 0;
header.encoding.encodingLimits.slice.maximum = nSlices-1;
header.encoding.encodingLimits.slice.center = floor((nSlices-1)/2);
header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = nReps-1;
header.encoding.encodingLimits.repetition.center = floor((nReps-1)/2);

% Parallel imaging info
% See https://gadgetron.discourse.group/t/ismrmrd-file-generation/141/5
header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 = R ;
header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2 = 1 ;
header.encoding.parallelImaging.calibrationMode = 'other' ;

% more info on our spiral
header.encoding.trajectoryDescription.identifier = 'ArchimedeanSpiral_0.8mm_R4_Lustig2008Design';
header.encoding.trajectoryDescription.userParameterLong.value = 0;
header.encoding.trajectoryDescription.userParameterLong.name = 'unused';
header.encoding.trajectoryDescription.userParameterDouble.value = 0;
header.encoding.trajectoryDescription.userParameterDouble.name = 'unused';

% from recon7 headers...to be adapted, maybe MPS now...
header.encoding.trajectoryDescription.userParameterString(1).name = 'Trajectory monitoring system vendor';
header.encoding.trajectoryDescription.userParameterString(1).value = 'IBT Prototype NI 7T';
header.encoding.trajectoryDescription.userParameterString(2).name = 'Monitoring software version';
header.encoding.trajectoryDescription.userParameterString(2).value = '2017.1.0000';

header.encoding.trajectoryDescription.userParameterString(3).name = 'kBasisType';
header.encoding.trajectoryDescription.userParameterString(3).value = 'RotatedToSubjectSliceGeometry';
header.encoding.trajectoryDescription.userParameterString(4).name = 'kBasisLabels';
header.encoding.trajectoryDescription.userParameterString(4).value = 'Measurement, Phase, Slice';
% header.encoding.trajectoryDescription.userParameterString(3).name = 'kBasisType';
% header.encoding.trajectoryDescription.userParameterString(3).value = 'sphericalHarmonicsAndSecondOrderSymmetricGradientConcomitant';
% header.encoding.trajectoryDescription.userParameterString(4).name = 'kBasisLabels';
% header.encoding.trajectoryDescription.userParameterString(4).value = 'B0, X, Y, Z, XY, ZY, 3Z^2-R^2, XZ, X^2-Y^2, 3YX^2-Y^3, XYZ, 5Z^2-R^2, 5Z^3-3ZR^2, 5Z^2-XR^2, X^2Z-YR^2, X^3-3XY^2, Z^2, X^2+Y^2, XZ, YZ';
% header.encoding.trajectoryDescription.userParameterString(5).name = 'kBasisMask';
% header.encoding.trajectoryDescription.userParameterString(5).value = '1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1';
% header.encoding.trajectoryDescription.userParameterString(6).name = 'kBasisUnits';
% header.encoding.trajectoryDescription.userParameterString(6).value = 'rad, rad/m, rad/m, rad/m, rad/m^2, rad/m^2, rad/m^2, rad/m^2, rad/m^2, rad/m^3, rad/m^3, rad/m^3, rad/m^3, rad/m^3, rad/m^3, rad/m^3, rad/m^2, rad/m^2, rad/m^2, rad/m^2';

% IBT-code specific flags
% header.userParameters.userParameterLong(1).syncPreScansAvailable = 0;
% header.userParameters.userParameterLong(2).isPartialFourier = 0;
% header.userParameters.userParameterDouble(1).triggerToAcquisitionDelay = 0; % separately merged...


%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

%% Write the dataset
dset.close();
