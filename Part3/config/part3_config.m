function cfg = part3_config()
%PART3_CONFIG Single source of truth for every Part 3 (alignment) parameter.
%
%   cfg = PART3_CONFIG() returns the Part 2 simulation config verbatim (so the
%   simulator layer is byte-identical to Part 2's), repointed to Part3's own
%   results folder, plus two additions:
%
%     1. The two new mechanism knobs (push-displacement and the growth-rate
%        field), stored as FLAT cfg.sim fields so the +align override loop
%        needs no nesting special-casing.
%     2. cfg.p3, the sweep / sensitivity / handoff settings for Phase A+B.
%
%   The measurement side (cfg.measure) is still Part 1's default_config(),
%   untouched: alignment is only meaningful if the same code scores real and
%   simulated lesions. Reuse-via-addpath, never copy.

thisDir = fileparts(mfilename('fullpath'));
cfg.root   = fileparts(thisDir);                    % Part3/
part2Root  = fullfile(fileparts(cfg.root), 'Part2');
part1Root  = fullfile(fileparts(cfg.root), 'Part1');

% sim_config() must be reachable before it is called.
if ~exist('sim_config', 'file')
    addpath(part2Root);
    addpath(fullfile(part2Root, 'config'));
end
cfg = sim_config();

% Repoint results/figures to Part3 (cfg.measure keeps Part 1's results path).
cfg.root = fileparts(thisDir);
cfg.paths.results = fullfile(cfg.root, 'results');
cfg.paths.figures = fullfile(cfg.root, 'results', 'figures');
if ~exist(cfg.paths.results, 'dir'), mkdir(cfg.paths.results); end
if ~exist(cfg.paths.figures, 'dir'), mkdir(cfg.paths.figures); end

% ---- New Part 3 mechanisms (flat, no nesting) ---------------------------
cfg.sim.displace     = false;   % arm A: push-displacement on division
cfg.sim.pushRate     = 1.0;     % prob a surrounded cell pushes (0..1)
cfg.sim.heteroOn     = false;   % arm B: spatial growth-rate heterogeneity
cfg.sim.heteroLambda = 120;     % correlation length (sites); swept
cfg.sim.heteroAmp    = 0.5;     % growth-rate amplitude

% ---- Part 3 alignment / sweep / sensitivity ----------------------------
cfg.p3 = struct();
cfg.p3.primary    = {'compactness', 'solidity', 'convexity'};
cfg.p3.targetsJson = fullfile(part1Root, 'results', 'calibration_targets_otsu.json');
cfg.p3.numWorkers = 4;
cfg.p3.seeds      = 3;
cfg.p3.levels     = 5;
cfg.p3.sobol.N     = 256;       % Saltelli base sample size (per matrix row)
cfg.p3.sobol.Nrun  = 384;       % lattice size for Sobol'/phase-diagram runs
cfg.p3.phase.nPts  = 16;        % phase-diagram grid points per axis
end
