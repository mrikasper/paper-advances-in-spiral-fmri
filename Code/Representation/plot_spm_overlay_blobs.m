% Script plot_spm_overlay_blobs
% mimick behavior of interactively opening CheckReg and adding coloured
% blobs of activation
%
%  plot_spm_overlay_blobs
%
%
%   See also

% Author:   Lars Kasper
% Created:  2019-09-06
% Copyright (C) 2019 Institute for Biomedical Engineering
%                    University of Zurich and ETH Zurich
%
% This file is part of the TAPAS UniQC Toolbox, which is released
% under the terms of the GNU General Public License (GPL), version 3.
% You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version).
% For further details, see the file COPYING or
%  <http://www.gnu.org/licenses/>.



for idxSubj = 2%:7
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Define and load data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    idxSess = 1;
    idxRecon = 1;
    
    contrastCombination = 'differentialT';%'differentialT';'separateT';'F','InOutComb', 'InOut'
    rotationParameters = [0 0 0]; %[-pi/3 0 0]; % pitch roll yaw in rad
    
    % position for good view
    centre = [0 -65 0];
    
    % interesting position native space
    %centre = [5.3183 -65.1380 1.2470];
    %interesting position after -pi/3 rotation
    %centre = [0.2240 -49.3819 -66.9862];
    
    % zoomed FOV in mm
    zoom_mm = 40;
    
    details = spifi_get_subject_details(idxSubj, [], idxSess, idxRecon);
    
    %% Setup basic plot
    spm_check_registration(details.preproc.anat.biascorrected);
    
    %% Mimick blob creation from spm_orthviews('context_menu','add_c_blobs',1)
    
    colours = [
        0 1 1
        1 1 0
        1 0 0
        ];
    
    
    switch contrastCombination
        case 'differentialT'
            Ics = [1 2]; % contrasts to query
        case 'separateT'
            Ics = [7 8]; % contrasts to query
        case 'F'
            Ics = 3;
        case 'InOutComb' % compares
            Ics = [1 1 1];
            idxSessions = [2 2 2];
            idxRecons = [1 2 3];
        case 'InOut' % compares
            Ics = [1 1 2 2];
            idxSessions = [2 2 2 2];
            idxRecons = [1 2 1 2];
            colours = [
                0 1 1
                1 1 0
                0 1 1
                1 1 0
                ];
    end
    
    
    
    details = spifi_get_subject_details(idxSubj, [], idxSess, idxRecon);
    nContrasts = numel(Ics);
    for c = 1:nContrasts
        
        switch contrastCombination
            case {'InOutComb', 'InOut'} % different GLM for comparison
                details = spifi_get_subject_details(idxSubj, [], idxSessions(c), idxRecons(c));
        end
        
        %% query SPM
        xSPM.swd = details.glm.session;
        xSPM.Ic = Ics(c); % contrast to evaluate
        xSPM.Im = []; % explicit mask
        xSPM.thresDesc = 'none'; % multiple comparison description
        xSPM.u = 0.001; % threshold for display, computation by u = spm_uc(u,df,STAT,R,n,S);
        xSPM.k = 0; % cluster extent threshold
        [SPM,xSPM] = spm_getSPM(xSPM);
        
        %% add blobs
        % FORMAT spm_orthviews('AddColouredBlobs',handle,XYZ,Z,mat,colour,name)
        % Add blobs from a pointlist to the image specified by the handle(s)
        % handle   - image number to add blobs to
        % XYZ      - blob voxel locations
        % Z        - blob voxel intensities
        % mat      - matrix from voxels to millimeters of blob.
        % colour   - the 3 vector containing the colour that the blobs should be
        % name     - a name for this blob
        % Several sets of blobs can be added in this way, and it uses full colour.
        % Although it may not be particularly attractive on the screen, the colour
        % blobs print well.
        
        
        handle = 1;
        
        spm_orthviews('AddColouredBlobs',handle,xSPM.XYZ,xSPM.Z,xSPM.M,colours(c,:),xSPM.title)
        spm_orthviews('AddColourBar',handle,c)
        
        % equalize axes
        if c == 1
            global st; % retrieve handle for later
        else % set all colors to same max
            st.vols{1}.blobs{c}.max =  st.vols{1}.blobs{1}.max;
        end
    end
    spm_orthviews('Redraw')
    
    
    %% reorient: rotate!
    % initialize reorient context menu
    spm_ov_reorient('context_init', 1)
    
    % manipulate edits in context GUI for reorientation
    for iAxis = 1:3
        st.vols{1}.reorient.e(3+iAxis).String = num2str(rotationParameters(iAxis));
    end
    
    % callback to updated orientation w/ current values in context GUI
    spm_ov_reorient('reorient', 1)
    
    %% Adjust view and zoom
    
    spm_orthviews('Reposition',centre)
    spm_orthviews('Zoom', zoom_mm)
    
    %% Adjust view to slice with most activation
    
    [idxSlice, centre] = spifi_find_most_activated_slices(...
        details.glm.tcon, xSPM.u);
    spm_orthviews('Reposition',centre);
  
    %% delete all blobs
    %spm_orthviews('RemoveBlobs',handle)
    
    spm_print(fullfile(details.representation.fig_spm_subject.pathSave, ...
        sprintf('%s_spm_blobs_%s.ps', details.subjectId, contrastCombination)));
end