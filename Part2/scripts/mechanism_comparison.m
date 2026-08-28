function mechanism_comparison(N, quick, armFilter)
%MECHANISM_COMPARISON Sweep four candidate mechanisms and score their morphology.
%
%   MECHANISM_COMPARISON(N, quick) grows tumours under each of four mechanism
%   arms (contact-inhibition loss, nutrient-driven fingering, enhanced
%   migration, clonal heterogeneity), sweeping the arm's parameter from
%   benign to malignant across `levels` values and `seeds` seeds. Every mask
%   is scored clean and degraded by the Part 1 bridge; results go to
%   results/mechanism_trajectories.csv and results/part2_handoff.json.
%
%   MECHANISM_COMPARISON(N, quick, armFilter) runs only the arm(s) named by
%   armFilter (a string or cell array of arm names) and writes each to
%   results/arm_<name>.csv instead of the combined file, with no handoff.
%   This chunking exists because a full 60-run 1024 sweep exceeds the 10-minute
%   background-task budget on this machine; run scripts/assemble_results.m
%   after the per-arm chunks to rebuild the combined CSV and handoff JSON.
%
%   The success criterion is whether any arm's trajectory in
%   (compactness, solidity, convexity) moves in the benign -> melanoma
%   direction; absolute proximity to the real medians is also reported, but
%   the clean CA baseline does not sit on the real benign point (Eden
%   roughness vs real large-scale irregularity + segmentation noise), so the
%   primary reading is the between-level DELTA, not absolute match.
%
%   Defaults: N = 1024 (production), quick = false. Pass N = 512 and
%   quick = true for a fast smoke run.

if nargin < 1 || isempty(N), N = 1024; end
if nargin < 2 || isempty(quick), quick = false; end
if nargin < 3, armFilter = ''; end
if ischar(armFilter) || isstring(armFilter)
    armFilter = cellstr(armFilter);
end

p2 = fileparts(fileparts(mfilename('fullpath')));
p1 = fullfile(fileparts(p2), 'Part1');
paths = {p2, fullfile(p2, 'config'), fullfile(p2, 'tests'), ...
         p1, fullfile(p1, 'config'), fullfile(p1, 'tests')};
addpath(paths{:});

cfg0 = sim_config();
cfg0.sim.N  = N;
cfg0.sim.nc = N / 8;
if N >= 1024
    cfg0.sim.stopEqDiameter = 600;
else
    cfg0.sim.stopEqDiameter = round(0.74 * N);
end

% ---- Arm definitions ----------------------------------------------------
% Each arm holds ONE mechanism at benign and sweeps its parameter(s) toward
% malignant; every other parameter stays at the benign base (COMP.ARM_BASE)
% so that the shape change is attributable to that mechanism alone.
arms = struct();

arms(1).name   = 'contact_inhibition';
arms(1).levels = { struct('beta', 0.95), struct('beta', 0.70), ...
                   struct('beta', 0.40), struct('beta', 0.20), ...
                   struct('beta', 0.05) };

arms(2).name   = 'nutrient_fingering';
arms(2).levels = { struct('L', 20), struct('L', 60), struct('L', 150), ...
                   struct('L', 300), struct('L', 600) };

arms(3).name   = 'enhanced_migration';
arms(3).levels = { struct('P_move', 0.02, 'A_i', 0.90), ...
                   struct('P_move', 0.06, 'A_i', 0.70), ...
                   struct('P_move', 0.12, 'A_i', 0.50), ...
                   struct('P_move', 0.20, 'A_i', 0.25), ...
                   struct('P_move', 0.35, 'A_i', 0.05) };

arms(4).name   = 'clonal_heterogeneity';
arms(4).levels = { struct('subclone_frac', 0.00), struct('subclone_frac', 0.05), ...
                   struct('subclone_frac', 0.10), struct('subclone_frac', 0.20), ...
                   struct('subclone_frac', 0.40) };

if ~isempty(armFilter)
    keep = false(size(arms));
    for a = 1:numel(arms)
        keep(a) = any(strcmp(armFilter, arms(a).name));
    end
    arms = arms(keep);
    if isempty(arms)
        error('armFilter did not match any arm name.');
    end
end

nLevels = cfg0.comparison.levels;
nSeeds  = cfg0.comparison.seeds;
if quick, nSeeds = 1; end

% ---- Flatten into a run list -------------------------------------------
runs = struct('arm', {}, 'level', {}, 'seed', {}, 'overrides', {});
for a = 1:numel(arms)
    for lv = 1:nLevels
        for s = 1:nSeeds
            r.arm       = arms(a).name;
            r.level     = lv;
            r.seed      = s;
            r.overrides = arms(a).levels{lv};
            runs(end+1) = r; %#ok<AGROW>
        end
    end
end

nRuns = numel(runs);
fprintf('Mechanism comparison: %d arms x %d levels x %d seeds = %d runs (N=%d)\n', ...
    numel(arms), nLevels, nSeeds, nRuns, N);

% Descriptors to collect (clean + degraded), from the shared Part 1 list.
descs = cfg0.measure.stats.descriptors;

useParfor = cfg0.comparison.useParfor && ~quick && nRuns > 1;
if useParfor && isempty(gcp('nocreate'))
    try
        parpool(cfg0.comparison.numWorkers);
    catch
        useParfor = false;   % fall back to serial
    end
end

t0 = tic;
rows = cell(nRuns, 1);

if useParfor
    parfor i = 1:nRuns
        addpath(paths{:});
        rows{i} = comp.run_one(runs(i), cfg0, descs);
    end
else
    for i = 1:nRuns
        rows{i} = comp.run_one(runs(i), cfg0, descs);
        fprintf('[%2d/%d] %-20s L%d seed%d  (%.0f s elapsed)\n', ...
            i, nRuns, runs(i).arm, runs(i).level, runs(i).seed, toc(t0));
    end
end

T = vertcat(rows{:});

% Write path depends on whether this is a filtered (per-arm chunk) run.
if isempty(armFilter)
    csvName = 'mechanism_trajectories.csv';
else
    csvName = ['arm_' strjoin(armFilter, '_') '.csv'];
end
csvPath = fullfile(cfg0.paths.results, csvName);
writetable(T, csvPath);

% Handoff is only meaningful once every arm is present (full run).
if isempty(armFilter)
    comp.write_handoff(T, N, cfg0.sim.stopEqDiameter, p1, ...
        fullfile(cfg0.paths.results, 'part2_handoff.json'));
end

fprintf('\nWrote %d rows to %s (%.0f s total)\n', nRuns, csvName, toc(t0));
if isempty(armFilter)
    fprintf('Wrote part2_handoff.json\n');
end
end
