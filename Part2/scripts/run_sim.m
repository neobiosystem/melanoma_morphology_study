function out = run_sim(phenotype, N)
%RUN_SIM Grow one virtual tumour, score it, and print the descriptors.
%
%   out = RUN_SIM(phenotype, N) applies the benign or melanoma preset, grows
%   a tumour to the target equivalent diameter, measures it clean and
%   degraded with the Part 1 bridge, and prints the headline descriptors.
%
%   Defaults: phenotype = 'benign', N = 512 (the CI/dev size). Use N = 1024
%   for production runs whose D_f scale span matches the real-lesion median.

if nargin < 1 || isempty(phenotype), phenotype = 'benign'; end
if nargin < 2 || isempty(N), N = 512; end

p2 = fileparts(fileparts(mfilename('fullpath')));
addpath(p2); addpath(fullfile(p2, 'config')); addpath(fullfile(p2, 'tests'));
p1 = fullfile(fileparts(p2), 'Part1');
addpath(p1); addpath(fullfile(p1, 'config')); addpath(fullfile(p1, 'tests'));

cfg = sim_config();
cfg.sim.N  = N;
cfg.sim.nc = N / 8;                       % keep 8 sites per coarse cell

% Target equivalent diameter: 600 for the production 1024 lattice (span 1.56
% decades, matches the real median), 0.74*N for the dev size (span ~1.40,
% matches the real p25). Stop BEFORE any border contact so the QC gate never
% fires on the simulation's own frame.
if N >= 1024
    cfg.sim.stopEqDiameter = 600;
else
    cfg.sim.stopEqDiameter = round(0.74 * N);
end

% Apply the phenotype preset onto the simulation parameters.
P = phen.presets();
prm = P.(phenotype);
for f = fieldnames(prm)'
    cfg.sim.(f{1}) = prm.(f{1});
end

cfg.sim.seed = 42;

t0 = tic;
S = ca.run_simulation(cfg);
growTime = toc(t0);

mask = ca.to_mask(S);
outClean = bridge.measure(mask, cfg);
outDeg   = bridge.measure(bridge.degrade(mask, cfg), cfg);

fprintf('\n=== %s tumour, N=%d, %d steps, %d cells (grew %.1f s) ===\n', ...
    phenotype, N, S.step, S.nCells, growTime);
fprintf('accepted=%d (%s)  DfFlag=%s  DfSpan=%.2f dec\n', ...
    outClean.accepted, outClean.rejectReason, char(outClean.DfFlag), ...
    outClean.DfSpanDecades);
fprintf('%-20s %8s %8s\n', 'descriptor', 'clean', 'degraded');
for d = cfg.measure.stats.descriptors
    nm = d{1};
    fprintf('%-20s %8.3f %8.3f\n', nm, outClean.(nm), outDeg.(nm));
end
fprintf('strippedAreaFrac=%.4f  holeAreaFrac=%.4f\n', ...
    outClean.preStrip.strippedAreaFraction, outClean.preStrip.holeAreaFraction);

out = struct('cfg', cfg, 'S', S, 'mask', mask, ...
             'clean', outClean, 'degraded', outDeg);
end
