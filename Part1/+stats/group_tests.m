function out = group_tests(values, groups, cfg)
%GROUP_TESTS Omnibus and post-hoc comparison of one descriptor across grades.
%
%   out = GROUP_TESTS(values, groups, cfg) runs Kruskal-Wallis, reports the
%   eta-squared-H effect size, and follows with Dunn's pairwise test.
%
%   Kruskal-Wallis rather than one-way ANOVA because these descriptors are
%   skewed and bounded (compactness has a hard floor at 1, asymmetry lives
%   in [0,1]), so the normality assumption is not met.
%
%   A p value alone is not enough here: with ~1900 nevi against ~200
%   in-situ, a trivial difference can reach significance. eta^2_H and the
%   pairwise rank-biserial effect sizes are what determine whether a
%   descriptor is actually useful as a calibration target.

values = values(:);
groups = categorical(groups(:));

ok     = isfinite(values) & ~isundefined(groups);
v      = values(ok);
g      = removecats(groups(ok));

out = struct('n', numel(v), 'numGroups', numel(categories(g)), ...
             'H', NaN, 'df', NaN, 'p', NaN, 'etaSqH', NaN, ...
             'postHoc', table());

if out.numGroups < 2 || out.n < cfg.stats.minNForStats
    return
end

[p, tbl] = kruskalwallis(v, g, 'off');

out.p  = p;
out.H  = tbl{2, 5};
out.df = tbl{2, 3};

k = out.numGroups;
N = out.n;
out.etaSqH = (out.H - k + 1) / (N - k);

out.postHoc = stats.dunn_test(v, g, cfg.stats.postHocCorrect);
end
