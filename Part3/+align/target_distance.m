function dist = target_distance(T, cfg0)
%TARGET_DISTANCE Z-like distance to the real invasive medians.
%
%   dist = TARGET_DISTANCE(T, cfg0) returns, for every row of T, the
%   normalised distance
%
%       dist = sqrt( sum_d ((clean_d - inv_median_d) / inv_iqr_d)^2 )
%
%   over the three primary descriptors (compactness, solidity, convexity),
%   reusing the exact formula of Part 2's COMP.WRITE_HANDOFF. This is the
%   single alignment objective: the lower, the closer the simulated shape is
%   to the real malignant signature. cfg0.p3.targetsJson points at the otsu
%   calibration targets.

targets = jsondecode(fileread(cfg0.p3.targetsJson));
primary = cfg0.p3.primary;

dist = zeros(height(T), 1);
for d = primary
    nm = d{1};
    v   = T.(['clean_' nm]);
    med = targets.descriptors.(nm).invasive.median;
    iqr = diff(targets.descriptors.(nm).invasive.iqr);
    dist = dist + ((v - med) / iqr).^2;
end
dist = sqrt(dist);
end
