function out = discriminative_power(values, groups, cfg)
%DISCRIMINATIVE_POWER AUC and ordinal concordance for one descriptor.
%
%   out = DISCRIMINATIVE_POWER(values, groups, cfg) returns:
%       auc, aucCI   - benign versus (in-situ + invasive), with a bootstrap CI
%       cIndex       - concordance across the three ordered grades
%       balancedAuc  - AUC re-estimated on grade-balanced subsamples
%
%   Raw classification accuracy is deliberately not reported. With nevi
%   outnumbering in-situ melanomas roughly ten to one, a classifier that
%   calls everything benign scores about 88%, so accuracy would flatter
%   every descriptor equally and inform nothing.
%
%   balancedAuc guards the same imbalance from the other side: it repeatedly
%   subsamples every grade down to the smallest one, so a descriptor that
%   only looks strong because nevi dominate is exposed.
%
%   The direction of every descriptor is normalised so AUC >= 0.5 and the
%   reported value always means "separates malignant from benign this well",
%   never "separates it backwards".

values = values(:);
groups = categorical(groups(:));

ok = isfinite(values) & ~isundefined(groups);
v  = values(ok);
g  = groups(ok);

out = struct('auc', NaN, 'aucCI', [NaN NaN], 'cIndex', NaN, ...
             'balancedAuc', NaN, 'balancedAucSd', NaN, 'direction', 1, 'n', numel(v));

if numel(v) < cfg.stats.minNForStats
    return
end

isMalignant = g ~= cfg.data.gradeOrder{1};
if all(isMalignant) || ~any(isMalignant)
    return
end

% Orient so a higher score always means more malignant.
direction = 1;
if median(v(isMalignant)) < median(v(~isMalignant))
    direction = -1;
end
out.direction = direction;
score = direction * v;

[~, ~, ~, auc] = perfcurve(double(isMalignant), score, 1);
out.auc = auc;

rng(cfg.stats.seed);
try
    [~, ~, ~, aucBoot] = perfcurve(double(isMalignant), score, 1, ...
        'NBoot', cfg.stats.numBootstrap, 'XVals', 'all');
    if size(aucBoot, 1) >= 1 && size(aucBoot, 2) >= 3
        out.aucCI = aucBoot(1, 2:3);
    end
catch
    % Bootstrap CI is optional; the point estimate still stands.
end

% ---- Ordinal concordance across all three grades ----------------------
out.cIndex = ordinal_c_index(score, g, cfg.data.gradeOrder);

% ---- Balanced-subsample sensitivity -----------------------------------
cats   = cfg.data.gradeOrder;
counts = arrayfun(@(k) nnz(g == cats{k}), 1:numel(cats));
nMin   = min(counts(counts > 0));

if nMin >= cfg.stats.minNForStats
    reps = zeros(cfg.stats.numBalancedReps, 1);
    for r = 1:cfg.stats.numBalancedReps
        idx = [];
        for k = 1:numel(cats)
            ck = find(g == cats{k});
            if numel(ck) >= nMin
                idx = [idx; ck(randperm(numel(ck), nMin))]; %#ok<AGROW>
            end
        end
        mSub = g(idx) ~= cats{1};
        if any(mSub) && ~all(mSub)
            [~, ~, ~, a] = perfcurve(double(mSub), score(idx), 1);
            reps(r) = a;
        else
            reps(r) = NaN;
        end
    end
    out.balancedAuc   = mean(reps, 'omitnan');
    out.balancedAucSd = std(reps, 'omitnan');
end
end

% -------------------------------------------------------------------------
function c = ordinal_c_index(score, g, order)
%ORDINAL_C_INDEX Share of cross-grade pairs ranked in the correct order.
rank = nan(size(g));
for k = 1:numel(order)
    rank(g == order{k}) = k;
end

valid = isfinite(rank) & isfinite(score);
rank  = rank(valid);
score = score(valid);

concordant = 0;
tied       = 0;
total      = 0;

for i = 1:numel(rank)-1
    higher = rank(i+1:end) > rank(i);
    lower  = rank(i+1:end) < rank(i);

    sj = score(i+1:end);
    si = score(i);

    concordant = concordant + nnz(higher & sj > si) + nnz(lower & sj < si);
    tied       = tied       + nnz((higher | lower) & sj == si);
    total      = total      + nnz(higher | lower);
end

if total == 0
    c = NaN;
else
    c = (concordant + 0.5 * tied) / total;
end
end
