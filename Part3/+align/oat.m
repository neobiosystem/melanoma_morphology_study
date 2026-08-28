function [T, ranking] = oat(armName, cfg0, N)
%OAT Phase B1: one-at-a-time screening of a winning arm's knobs.
%
%   [T, ranking] = OAT(armName, cfg0, N) sweeps each free knob of armName one
%   at a time (the others at the arm base), scores dist-to-melanoma per run,
%   and returns the run table T plus a ranking table (knob, maxAbsDelta)
%   sorted by the largest |delta dist| from baseline. The top knobs are the
%   two carried into the reduced Sobol' and phase diagram.

if nargin < 3 || isempty(N), N = 512; end

cfg = cfg0;
cfg.sim.N = N; cfg.sim.nc = N/8;
if N >= 1024, cfg.sim.stopEqDiameter = 600; else, cfg.sim.stopEqDiameter = round(0.74*N); end

K = align.knob_levels(armName, cfg);
knobs = fieldnames(K);
nSeeds = cfg.p3.seeds;

runs = struct('arm', {}, 'knob', {}, 'level', {}, 'seed', {}, 'overrides', {});
for s = 1:nSeeds
    r.arm = armName; r.knob = 'baseline'; r.level = 0; r.seed = s; r.overrides = struct();
    runs(end+1) = r; %#ok<AGROW>
end
for ki = 1:numel(knobs)
    nm = knobs{ki};
    for lv = 1:numel(K.(nm))
        for s = 1:nSeeds
            r.arm = armName; r.knob = nm; r.level = lv; r.seed = s;
            r.overrides = struct(nm, K.(nm)(lv));
            runs(end+1) = r; %#ok<AGROW>
        end
    end
end

nRuns = numel(runs);
descs = cfg.measure.stats.descriptors;
rows = cell(nRuns, 1);

useParfor = cfg.comparison.useParfor && nRuns > 1;
if useParfor && isempty(gcp('nocreate'))
    try, parpool(cfg.p3.numWorkers); catch, useParfor = false; end
end
if useParfor
    parfor i = 1:nRuns
        rows{i} = align.run_one(runs(i), cfg, descs);
    end
else
    for i = 1:nRuns
        rows{i} = align.run_one(runs(i), cfg, descs);
    end
end

T = vertcat(rows{:});
T.dist = align.target_distance(T, cfg);
T.knob  = string({runs.knob}');
T.level = [runs.level]';
T.seed  = [runs.seed]';

% Ranking: per knob, max over levels of |mean-level dist - baseline dist|.
base = mean(T.dist(strcmp(T.knob, 'baseline')), 'omitnan');
ranking = table('Size', [0 2], 'VariableTypes', {'string', 'double'}, ...
    'VariableNames', {'knob', 'maxAbsDelta'});
for ki = 1:numel(knobs)
    nm = knobs{ki};
    idx = strcmp(T.knob, nm);
    if ~any(idx), continue; end
    lvls = unique(T.level(idx));
    maxAbs = 0;
    for lv = lvls'
        d = mean(T.dist(idx & T.level == lv), 'omitnan');
        maxAbs = max(maxAbs, abs(d - base));
    end
    ranking = [ranking; {string(nm), maxAbs}]; %#ok<AGROW>
end
ranking = sortrows(ranking, 'maxAbsDelta', 'descend');
fprintf('OAT baseline dist = %.3f\n', base);
disp(ranking);
end
