function cfg = sim_config()
%SIM_CONFIG Single source of truth for every Part 2 simulation parameter.
%
%   cfg = SIM_CONFIG() returns a struct consumed by +ca, +phen, +bridge and
%   the scripts under Part2/scripts. Two design rules mirror Part 1:
%
%   1. The MEASUREMENT side is not reimplemented. cfg.measure is Part 1's
%      default_config() verbatim, so every simulated mask is scored by the
%      exact same code that scored the real lesions. This is what makes the
%      Part 3 alignment meaningful: a copy would silently drift.
%
%   2. No numeric constant is hard-coded in a step function; every value a
%      reviewer might question lives here.

% ---- Paths -------------------------------------------------------------
thisDir     = fileparts(mfilename('fullpath'));
cfg.root    = fileparts(thisDir);                       % Part2/
part1Root   = fullfile(fileparts(cfg.root), 'Part1');

cfg.paths.results  = fullfile(cfg.root, 'results');
cfg.paths.figures  = fullfile(cfg.root, 'results', 'figures');
if ~exist(cfg.paths.results, 'dir'), mkdir(cfg.paths.results); end
if ~exist(cfg.paths.figures, 'dir'), mkdir(cfg.paths.figures); end

% ---- Measurement bridge ------------------------------------------------
% Part 1's config must be reachable before default_config() is called.
if ~exist('default_config', 'file')
    addpath(part1Root);
    addpath(fullfile(part1Root, 'config'));
end
cfg.measure = default_config();          % the unmodified Part 1 pipeline
cfg.measure.paths.results = fullfile(part1Root, 'results'); % keep targets readable

% ---- Lattice -----------------------------------------------------------
% Production 1024x1024 so that a tumour grown to D_eq ~600 sites spans 1.56
% decades of box-counting scale, matching the median real lesion (574 px ->
% 1.555 decades). 512x512 is the CI/dev size: a D_eq ~400 tumour spans only
% 1.40 decades, i.e. the real-lesion p25.
cfg.sim.N               = 1024;    % lattice side (sites)
cfg.sim.nc              = 128;     % coarse nutrient-grid side (N/nc = 8 sites/cell)
cfg.sim.connectivity    = 8;       % 4 = von Neumann, 8 = Moore
cfg.sim.diagWeight      = 1/sqrt(2); % isotropy correction for Moore diagonals

% ---- Growth / proliferation -------------------------------------------
% T_cell is FIXED so that alpha and T_cell (two knobs that both set the
% effective division rate) never trade against each other: that is a
% structural non-identifiability and is removed here, not in Part 3.
cfg.sim.T_cell          = 4;       % min interval between divisions, FIXED
cfg.sim.alpha           = 0.5;     % base division probability per step
cfg.sim.beta            = 0.95;    % contact inhibition in [0,1]; 1 = strong
cfg.sim.initRadius      = 8;       % initial cluster radius (sites)

% ---- Migration ---------------------------------------------------------
cfg.sim.P_move          = 0.02;    % base migration probability per step
cfg.sim.A_i             = 0.90;    % adhesion threshold in [0,1]; 1 = sticky
cfg.sim.chi             = 0.0;     % chemotactic bias strength (0 = off)

% ---- Apoptosis / necrosis ---------------------------------------------
cfg.sim.delta           = 1e-4;    % basal death probability per cell per step
cfg.sim.pStarve         = 0.0;     % extra death prob for C < C_min (necrosis)

% ---- Nutrient field (hybrid reaction-diffusion) ------------------------
% Steady-state equation D*laplacian(C) = k*rho*C, i.e. L^2 laplacian(C) = rho*C
% with L = sqrt(D/k) the diffusion length in fine sites. rho is the local
% occupied-cell fraction on the coarse grid. Boundary is Dirichlet C = C0.
cfg.sim.nutrientCoupling = false;  % true enables the full hybrid model
cfg.sim.nutrientInterval = 10;     % solve the field every N steps
cfg.sim.L               = 60;      % diffusion length (fine sites); swept in arm B
cfg.sim.C0              = 1.0;     % far-field nutrient (Dirichlet)
cfg.sim.C_min           = 0.2;     % below this: necrosis (extra death)
cfg.sim.C_q             = 0.3;     % below this: proliferation quiescent

% ---- Clonal heterogeneity (mechanism arm D) ----------------------------
cfg.sim.subclone.frac       = 0.0; % fraction of cells in the advantaged clone
cfg.sim.subclone.alphaMult   = 3.0; % advantaged clone's alpha multiplier

% ---- Run control -------------------------------------------------------
cfg.sim.maxSteps        = 6000;    % hard cap; stopEqDiameter usually fires first
cfg.sim.stopEqDiameter  = 600;     % stop growing when dEq reaches this (sites)
cfg.sim.seed            = 1;
cfg.sim.snapshotInterval = 50;     % metric history every N steps
cfg.sim.maskSaveInterval = 250;    % full-mask snapshots every N steps (figures)

% ---- Boundary-noise degradation for comparability with real masks -------
% Real masks carry Otsu segmentation noise; simulated masks do not. Each mask
% is therefore also measured after a first-order +/-1 px stochastic boundary
% perturbation, so Part 3 can choose. Main conclusions use between-group
% differences, which the noise mostly cancels.
cfg.degrade.frac        = 0.04;    % fraction of boundary pixels perturbed

% ---- Comparison sweep --------------------------------------------------
cfg.comparison.numWorkers = 4;     % parpool size: 6 workers at 1024^2 exceeds
                                   % this machine's 8 GB RAM (OOM kill); 4 is safe.
cfg.comparison.useParfor  = true;
cfg.comparison.levels     = 5;
cfg.comparison.seeds      = 3;
end
