function run_discovery(N, quick, armFilter)
%RUN_DISCOVERY Phase A: compare the two new mechanisms.
%
%   run_discovery(512)          full discovery at 512^2 (2 arms x 5 x 3 = 30 runs).
%   run_discovery(512, true)    smoke: 1 seed per level.
%   run_discovery(512, false, 'mechanical_displacement')  single arm, no handoff.
%
%   Outputs: results/discovery_trajectories.csv (or results/arm_<name>.csv)
%   and results/part3_handoff.json (full run only).
if nargin < 1, N = 512; end
if nargin < 2, quick = false; end
if nargin < 3, armFilter = ''; end

p3 = fileparts(fileparts(mfilename('fullpath')));
addpath(p3); addpath(fullfile(p3, 'config'));

align.discovery(N, quick, armFilter);
end
