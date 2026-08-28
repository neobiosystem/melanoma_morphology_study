function T = phase_diagram(armName, p1, p2, cfg0, nPts)
%PHASE_DIAGRAM Phase B3: 2D grid of the two dominant knobs, response = dist.
%
%   T = PHASE_DIAGRAM(armName, p1, p2, cfg0, nPts) builds an nPts x nPts grid
%   over knobs p1 and p2 (the two dominant from OAT), runs one seed per cell
%   (a phase map, not a Monte Carlo average), and returns a table with columns
%   <p1>, <p2>, dist, and the clean primary descriptors.
%
%   The real benign/melanoma points are NOT overlaid here; that is a figure
%   concern (make_figures reads cfg0.p3.targetsJson directly).

if nargin < 5 || isempty(nPts), nPts = cfg0.p3.phase.nPts; end

cfg = cfg0;
cfg.sim.N  = cfg0.p3.sobol.Nrun;
cfg.sim.nc = cfg.sim.N / 8;
cfg.sim.stopEqDiameter = round(0.74 * cfg.sim.N);

K = align.knob_levels(armName, cfg);
v1 = linspace(min(K.(p1)), max(K.(p1)), nPts);
v2 = linspace(min(K.(p2)), max(K.(p2)), nPts);

runs = struct('arm', {}, 'level', {}, 'seed', {}, 'overrides', {});
nRun = 0;
for a = 1:nPts
    for b = 1:nPts
        nRun = nRun + 1;
        runs(nRun).arm = armName; %#ok<AGROW>
        runs(nRun).level = nRun;
        runs(nRun).seed = 1;
        runs(nRun).overrides = struct(char(p1), v1(a), char(p2), v2(b));
    end
end

descs = cfg.measure.stats.descriptors;
rows = cell(nRun, 1);
useParfor = cfg.comparison.useParfor && nRun > 1;
if useParfor && isempty(gcp('nocreate'))
    try, parpool(cfg.p3.numWorkers); catch, useParfor = false; end
end
if useParfor
    parfor i = 1:nRun
        rows{i} = align.run_one(runs(i), cfg, descs);
    end
else
    for i = 1:nRun
        rows{i} = align.run_one(runs(i), cfg, descs);
    end
end

T = vertcat(rows{:});
T.dist = align.target_distance(T, cfg);
T.(p1) = repelem(v1(:), nPts, 1);
T.(p2) = repmat(v2(:), nPts, 1);
end
