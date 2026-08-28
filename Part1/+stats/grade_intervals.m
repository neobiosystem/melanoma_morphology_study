function T = grade_intervals(values, groups, cfg)
%GRADE_INTERVALS Per-grade descriptive intervals for one descriptor.
%
%   T = GRADE_INTERVALS(values, groups, cfg) returns one row per grade with
%   n, median, IQR, the 2.5-97.5 percentile reference interval, and mean/SD.
%
%   This is the headline deliverable of Part 1: the morphological interval
%   of each severity grade, and the form in which Part 3 consumes the
%   calibration targets.
%
%   Median and IQR lead rather than mean and SD because these distributions
%   are skewed and the rare grades are small, where a mean is outlier
%   dominated. Mean and SD are still reported so the numbers can be compared
%   against prior literature that quotes them.
%
%   Reference intervals are suppressed below cfg.stats.minNForRefInt: a
%   2.5-97.5 percentile estimated from a handful of lesions is noise
%   presented as a bound, which is worse than reporting nothing.

values = values(:);
groups = categorical(groups(:));

cats = cfg.data.gradeOrder;
rows = {};

for k = 1:numel(cats)
    sel = groups == cats{k} & isfinite(values);
    v   = values(sel);
    n   = numel(v);

    r.grade = string(cats{k});
    r.n     = n;

    if n < cfg.stats.minNForStats
        [r.median, r.q1, r.q3, r.mean, r.sd, r.refLo, r.refHi] = deal(NaN);
        r.note = "insufficient_n";
        rows{end+1} = r; %#ok<AGROW>
        continue
    end

    q = quantile(v, [0.25 0.5 0.75]);
    r.median = q(2);
    r.q1     = q(1);
    r.q3     = q(3);
    r.mean   = mean(v);
    r.sd     = std(v);

    if n >= cfg.stats.minNForRefInt
        ref = prctile(v, cfg.stats.refIntervalPct);
        r.refLo = ref(1);
        r.refHi = ref(2);
        r.note  = "";
    else
        r.refLo = NaN;
        r.refHi = NaN;
        r.note  = "n_below_reference_interval_threshold";
    end

    rows{end+1} = r; %#ok<AGROW>
end

T = struct2table([rows{:}]);
end
