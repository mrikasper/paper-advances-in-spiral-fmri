% id, FOVm, FOVp, FOVs, dxm, dxp, dxs, Rp, Rs, nIl, maxG, maxSR, maxTacq, in/out
% OR (# added echoes), type, scheme3D, nPlanesPerShot
% InOut/OutIn-Trajektorien
% For 2D trajectories, set FOVs and dxs to zero.
% scheme3D: 'etagere', 'stack' (of spirals),'yarnball','shells','cones',
% 'blippedCones','gaSpiral' (leave empty or 'stack' for 2D traj)


%% 2D Spiral fMRI

% only 1 ms slower than SR 200, but no PNS
40, 0.23, 0.23, 0.00, 0.8e-3, 0.8e-3, 0, 4, 1, 1, 31e-3, 160, 0, out, minTime, stack, 1
41, 0.23, 0.23, 0.00, 1.5e-3, 1.5e-3, 0, 4, 1, 1, 31e-3, 160, 0, inout, minTime, stack, 1

% Using new gradient system, parallel mode
10, 0.23, 0.23, 0.00, 0.8e-3, 0.8e-3, 0, 4, 1, 1, 200e-3, 600, 0, out, minTime, stack, 1
11, 0.23, 0.23, 0.00, 1.5e-3, 1.5e-3, 0, 4, 1, 1, 200e-3, 600, 0, inout, minTime, stack, 1

% Using new gradient system, serial mode
30, 0.23, 0.23, 0.00, 0.8e-3, 0.8e-3, 0, 4, 1, 1, 100e-3, 1200, 0, out, minTime, stack, 1
31, 0.23, 0.23, 0.00, 1.5e-3, 1.5e-3, 0, 4, 1, 1, 100e-3, 1200, 0, inout, minTime, stack, 1


%% Competitive EPI (Same k-space area as spiral), old and new slew
1, 0.23, 0.23, 0, 0.8e-3, 0.8e-3, 0, 4, 1, 1, 31e-3, 200, 0, 0, EPI, stack, 1
2, 0.23, 0.23, 0, 0.9e-3, 0.9e-3, 0, 4, 1, 1, 31e-3, 200, 0, 0, EPI, stack, 1
3, 0.23, 0.23, 0, 0.9e-3, 0.9e-3, 0, 4, 1, 1, 31e-3, 160, 0, 0, EPI, stack, 1
4, 0.23, 0.23, 0, 2.05e-3, 2.05e-3, 0, 1, 1, 1, 31e-3, 160, 0, out, minTime, stack, 1
5, 0.23, 0.23, 0, 2.5e-3, 2.5e-3, 0, 1, 1, 1, 31e-3, 160, 0, 0, EPI, stack, 1

%% old full slew spirals
20, 0.23, 0.23, 0.00, 0.8e-3, 0.8e-3, 0, 4, 1, 1, 31e-3, 200, 0, out, minTime, stack, 1
21, 0.23, 0.23, 0.00, 1.5e-3, 1.5e-3, 0, 4, 1, 1, 31e-3, 200, 0, inout, minTime, stack, 1

%% High Res Flow Trajs
51, 0.204, 0.204, 0.00, 0.1e-3, 0.1e-3, 0, 512, 1, 1, 70e-3, 200, 0, out, minTime, stack, 1
52, 0.051, 0.051, 0.00, 0.1e-3, 0.1e-3, 0, 128, 1, 1, 70e-3, 200, 0, out, minTime, stack, 1
53, 0.204, 0.204, 0.00, 0.08862e-3, 0.08862e-3, 0, 512, 1, 1, 70e-3, 200, 0, out, minTime, stack, 1
54, 0.204, 0.204, 0.00, 0.1063e-3, 0.1063e-3, 0, 512, 1, 1, 70e-3, 200, 0, out, minTime, stack, 1
55, 0.204, 0.204, 0.00, 0.1063e-3, 0.1063e-3, 0, 512, 1, 1, 300e-3, 200, 0, out, minTime, stack, 1
56, 0.204, 0.204, 0.00, 0.1063e-3, 0.1063e-3, 0, 512, 1, 1, 200e-3, 600, 0, out, minTime, stack, 1
57, 0.23, 0.230, 0.00, 0.7e-3, 0.7e-3, 0, 4, 1, 1, 200e-3, 600, 0, out, minTime, stack, 1
58, 0.23, 0.230, 0.00, 0.5e-3, 0.5e-3, 0, 4, 1, 1, 200e-3, 600, 0, out, minTime, stack, 1
59, 0.23, 0.230, 0.00, 0.5e-3, 0.5e-3, 0, 1, 1, 1, 200e-3, 600, 0, out, minTime, stack, 1
