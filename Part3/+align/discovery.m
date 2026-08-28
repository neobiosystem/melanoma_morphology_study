function discovery(N, quick, armFilter)
%DISCOVERY Phase A: compare the two new mechanisms (push vs growth field).
%
%   DISCOVERY(N, quick) sweeps both Part 3 mechanism arms across `levels`
%   levels and `seeds` seeds at lattice size N, scores every mask clean and
%   degraded, and writes results/discovery_trajectories.csv plus the per-run
%   distance-to-melanoma column. It also writes results/part3_handoff.json.
%
%   DISCOVERY(N, quick, armFilter) runs only the named arm(s), writing
%   results/arm_<name>.csv with no handoff, mirroring Part 2's chunking.
%
%   The success criterion is whether either arm drives
%   (compactness, solidity, convexity) toward the melanoma medians and drops
%   the 3-descriptor distance below Part 2's ~1.15-1.17.

if nargin < 1 || isempty(N), N = 512; end
if nargin < 2 || isempty(quick), quick = false; end
if nargin < 3, armFilter = ''; end
if ischar(armFilter) || isstring(armFilter), armFilter = cellstr(armFilter); end
% cellstr('') yields {''} (non-empty) -- normalise so '' and {} both mean "no filter".
if numel(armFilter) == 1 && isempty(char(armFilter{1})), armFilter = {}; end

p3 = fileparts(fileparts(mfilename('fullpath')));
p2 = fullfile(fileparts(p3), 'Part2');
p1 = fullfile(fileparts(p3), 'Part1');
paths = {p3, fullfile(p3, 'config'), fullfile(p3, 'tests'), ...
         p2, fullfile(p2, 'config'), fullfile(p2, 'tests'), ...
         p1, fullfile(p1, 'config'), fullfile(p1, 'tests')};
addpath(paths{:});

cfg0 = part3_config();
cfg0.sim.N  = N;
cfg0.sim.nc = N / 8;
if N >= 1024
    cfg0.sim.stopEqDiameter = 600;
else
    cfg0.sim.stopEqDiameter = round(0.74 * N);
end

arms = align.arm_defs(cfg0);
if ~isempty(armFilter)
    keep = false(size(arms));
    for a = 1:numel(arms), keep(a) = any(strcmp(armFilter, arms(a).name)); end
    arms = arms(keep);
    if isempty(arms), error('armFilter did not match any arm name.'); end
end

nLevels = cfg0.p3.levels;
nSeeds  = cfg0.p3.seeds;
if quick, nSeeds = 1; end

runs = struct('arm', {}, 'level', {}, 'seed', {}, 'overrides', {});
for a = 1:numel(arms)
    for lv = 1:nLevels
        for s = 1:nSeeds
            r.arm = arms(a).name; r.level = lv; r.seed = s;
            r.overrides = arms(a).levels{lv};
            runs(end+1) = r; %#ok<AGROW>
        end
    end
end

nRuns = numel(runs);
fprintf('Discovery: %d arms x %d levels x %d seeds = %d runs (N=%d)\n', ...
    numel(arms), nLevels, nSeeds, nRuns, N);

descs = cfg0.measure.stats.descriptors;

useParfor = cfg0.comparison.useParfor && ~quick && nRuns > 1;
if useParfor && isempty(gcp('nocreate'))
    try, parpool(cfg0.p3.numWorkers); catch, useParfor = false; end
end

t0 = tic;
rows = cell(nRuns, 1);
if useParfor
    parfor i = 1:nRuns
        addpath(paths{:});
        rows{i} = align.run_one(runs(i), cfg0, descs);
    end
else
    for i = 1:nRuns
        rows{i} = align.run_one(runs(i), cfg0, descs);
        fprintf('[%2d/%d] %-24s L%d seed%d (%.0f s)\n', i, nRuns, runs(i).arm, ...
            runs(i).level, runs(i).seed, toc(t0));
    end
end

T = vertcat(rows{:});
T.dist = align.target_distance(T, cfg0);

if isempty(armFilter)
    csvName = 'discovery_trajectories.csv';
else
    csvName = ['arm_' strjoin(armFilter, '_') '.csv'];
end
csvPath = fullfile(cfg0.paths.results, csvName);
writetable(T, csvPath);

if isempty(armFilter)
    align.write_handoff(T, cfg0, fullfile(cfg0.paths.results, 'part3_handoff.json'));
end

fprintf('\nWrote %d rows to %s (%.0f s total)\n', nRuns, csvName, toc(t0));
if isempty(armFilter), fprintf('Wrote part3_handoff.json\n'); end
end
