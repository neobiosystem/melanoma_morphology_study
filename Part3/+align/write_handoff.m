function H = write_handoff(T, cfg0, outPath)
%WRITE_HANDOFF Write part3_handoff.json from a discovery result table.
%
%   H = WRITE_HANDOFF(T, cfg0, outPath) computes, for each arm, the run
%   closest (in z-like distance across the three primary descriptors) to the
%   real invasive median, and writes the handoff struct to outPath as JSON.
%   The struct records the winning mechanism and the alignment error, as the
%   Part 4 input.

targets = jsondecode(fileread(cfg0.p3.targetsJson));
primary = cfg0.p3.primary;

H = struct();
H.source = 'Part3 mechanism + alignment';
H.N = cfg0.sim.N;
H.stopEqDiameter = cfg0.sim.stopEqDiameter;
H.generated = char(datetime('now'));
H.primaryDescriptors = {primary};
H.benignTargets = struct();
H.melanomaTargets = struct();
for d = primary
    nm = d{1};
    H.benignTargets.(nm)   = targets.descriptors.(nm).benign.median;
    H.melanomaTargets.(nm) = targets.descriptors.(nm).invasive.median;
end

arms = categories(categorical(T.arm));
H.best = struct();
for a = 1:numel(arms)
    nm = arms{a};
    idx = strcmp(T.arm, nm) & (T.accepted == 1);
    if ~any(idx), continue; end
    sub = T(idx, :);
    dist = align.target_distance(sub, cfg0);
    [dmin, bestIdx] = min(dist);
    row = sub(bestIdx, :);
    H.best.(nm) = struct( ...
        'level', row.level, 'seed', row.seed, ...
        'distanceToMelanoma', dmin, ...
        'clean_compactness', row.clean_compactness, ...
        'clean_solidity', row.clean_solidity, ...
        'clean_convexity', row.clean_convexity, ...
        'degraded_compactness', row.degraded_compactness, ...
        'degraded_solidity', row.degraded_solidity, ...
        'degraded_convexity', row.degraded_convexity);
end

% Winner = the arm with the smallest distance (the one that goes to Phase B).
dminAll = inf; winner = '';
for a = 1:numel(arms)
    nm = arms{a};
    if isfield(H.best, nm) && H.best.(nm).distanceToMelanoma < dminAll
        dminAll = H.best.(nm).distanceToMelanoma;
        winner = nm;
    end
end
if ~isempty(winner)
    H.winner = winner;
    H.winnerDistance = dminAll;
end

fid = fopen(outPath, 'w');
fwrite(fid, jsonencode(H));
fclose(fid);
end
