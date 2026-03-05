%% =============================
addpath /storage/work/krf5429/spm 

% --- Define type and region ---
type = 'full';    
ROI  = 'MFG'; 

%% =============================

% General scanning specifications
n_scans = 157;
n_regions = 2;
n_inputs = 2;

% Hypothesis connectivity indices
% A: target, source
A_idxs = [2,1;
          1,2;
          1,1;
          2,2];

% B: input, target, source
B_idxs = [2,2,1;
          2,1,2;
          2,1,1];

% C: region, input
C_idxs = [1,1;
          2,1];

% Initialize DCM structure
DCM = struct();

% Start to populate with things that don't change between subjects and
% sessions
DCM.Y.dt = 2;       % TR (s)
DCM.Y.X0 = [];         % no confounds
DCM.v    = n_scans;    % number of scans
DCM.n    = n_regions;  % number of regions

% Initialize connectivity matrices
DCM.a = zeros(n_regions);                    % intrinsic
DCM.b = zeros(n_regions,n_regions,n_inputs); % modulatory
DCM.c = zeros(n_regions,n_inputs);           % driving

% Populate a-matrix (target, source)
for i = 1:size(A_idxs,1)
    target = A_idxs(i,1);
    source = A_idxs(i,2);
    DCM.a(target, source) = 1;
end

% Populate b-matrix (input, target, source)
for i = 1:size(B_idxs,1)
    input  = B_idxs(i,1);
    target = B_idxs(i,2);
    source = B_idxs(i,3);
    DCM.b(target, source, input) = 1;
end

% Populate c-matrix (region, input)
for i = 1:size(C_idxs,1)
    region = C_idxs(i,1);
    input  = C_idxs(i,2);
    DCM.c(region, input) = 1;
end

% TE and delays (SPM standard delay is TR/2)
DCM.TE = 0.04;
DCM.delays(1:n_regions) = 1;

% Task-based BOLD options
DCM.options.nonlinear  = 0;       % bilinear
DCM.options.two_state  = 0;       % single-state
DCM.options.stochastic = 0;       % deterministic
DCM.options.centre     = 0;       % use U exactly as given
DCM.options.induced    = 0;       % no induced responses
DCM.options.analysis   = 'time';  % time-domain BOLD DCM
DCM.options.nograph    = 1;       % display output plots
DCM.options.maxit      = 512;     % maximum VL iterations

% Define region and input names
region_names = {ROI, 'PCC'};              % regions
input_names  = {'GNG', 'Inhibition'};   % experimental inputs

% --- Build the filename dynamically ---
filename = sprintf('../Data/%s_%s_%s_roi.mat', sub, ROI, type);

% --- Load the file ---
load(filename);

% --- Message ---
msg = sprintf('DCM for %s, %s, and ROI %s', sub, ROI, type);
disp(msg)

% --- Time series data
Y = struct2array(Y);
DCM.Y.y  = Y;  

% --- xY structure
for r = 1:n_regions
    DCM.xY(r).name  = region_names{r};     % region name
    DCM.xY(r).u     = Y(:,r);              % time series
    DCM.xY(r).DT    = DCM.Y.dt;            % TR
    DCM.xY(r).XYZ   = [0 0 0];             % dummy coordinates
    DCM.xY(r).def   = 'manual';            % required
    DCM.xY(r).Sname = 'Session_1';         % dummy session
end

% --- Experimental inputs
DCM.U.u    = U;
DCM.U.name = input_names;                  % input names
DCM.U.dt   = DCM.Y.dt;

% Save and estimate
output   = sprintf('../Output/DCM_%s_%s_%s_roi.mat', sub, ROI, type);

% --- Save the DCM structure ---
save(output, 'DCM');

% --- Estimate the DCM ---
spm_dcm_estimate(output);

% Load fitted DCM back in
load(output);

% Save necessary variables for R
Cp = full(Cp);   

Ep_A = Ep.A;
Ep_B = Ep.B;
Ep_C = Ep.C;
DCM_y = DCM.y;
transit = full(Ep.transit);
decay = full(Ep.decay);
epsilon = full(Ep.epsilon);


save(sprintf('../Output/DCM_%s_%s_%s_roi_out.mat', sub, ROI, type), 'Ep_A', 'Ep_B', 'Ep_C', 'Cp', 'DCM_y', ...
    "transit","decay","epsilon");


