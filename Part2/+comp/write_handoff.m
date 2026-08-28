function H = write_handoff(T, N, stopEqDiameter, p1, outPath)
%WRITE_HANDOFF Write part2_handoff.json from a comparison result table.
%
%   H = WRITE_HANDOFF(T, N, stopEqDiameter, p1, outPath) computes, for each
%   arm, the run closest (in z-like distance across the three primary
%   descriptors) to the real invasive median, and writes the handoff struct to
%   outPath as JSON. p1 is the Part 1 root (for the calibration-targets JSON).

targets = jsondecode(fileread(fullfile(p1, 'results', 'calibration_targets_otsu.json')));
primary = {'compactness', 'solidity', 'convexity'};

H = struct();
H.source = 'Part2 mechanism comparison';
H.N = N;
H.stopEqDiameter = stopEqDiameter;
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
    dist = zeros(height(sub), 1);
    for d = primary
        dn = d{1};
        v   = sub.(['clean_' dn]);
        med = targets.descriptors.(dn).invasive.median;
        iqr = diff(targets.descriptors.(dn).invasive.iqr);
        dist = dist + ((v - med) / iqr).^2;
    end
    dist = sqrt(dist);
    [~, bestIdx] = min(dist);
    row = sub(bestIdx, :);
    H.best.(nm) = struct( ...
        'level', row.level, 'seed', row.seed, ...
        'distanceToMelanoma', dist(bestIdx), ...
        'clean_compactness', row.clean_compactness, ...
        'clean_solidity', row.clean_solidity, ...
        'clean_convexity', row.clean_convexity, ...
        'degraded_compactness', row.degraded_compactness, ...
        'degraded_solidity', row.degraded_solidity, ...
        'degraded_convexity', row.degraded_convexity);
end

fid = fopen(outPath, 'w');
fwrite(fid, jsonencode(H));
fclose(fid);
end
