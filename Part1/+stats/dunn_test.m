function R = dunn_test(values, groups, correction)
%DUNN_TEST Non-parametric pairwise post-hoc following Kruskal-Wallis.
%
%   R = DUNN_TEST(values, groups) returns a table of pairwise comparisons
%   with z statistics, raw and adjusted p values, and the rank-biserial
%   effect size.
%
%   MATLAB's multcompare on a kruskalwallis result performs a Tukey-style
%   comparison, which is not the correct companion to Kruskal-Wallis. Dunn's
%   test is, because it compares mean RANKS using the pooled rank variance
%   from the omnibus test, with the tie correction that dermoscopic
%   descriptor data needs.
%
%   See also STATS.GROUP_TESTS, STATS.ADJUST_P.

if nargin < 3
    correction = 'holm';
end

values = values(:);
groups = categorical(groups(:));

ok     = isfinite(values) & ~isundefined(groups);
values = values(ok);
groups = removecats(groups(ok));

cats = categories(groups);
k    = numel(cats);
N    = numel(values);

R = table();
if k < 2 || N < 3
    return
end

ranks = tiedrank(values);

% Tie correction for the pooled rank variance.
[~, ~, idx] = unique(values);
tieCounts   = accumarray(idx, 1);
tieTerm     = sum(tieCounts.^3 - tieCounts) / (12 * (N - 1));

meanRank = zeros(k, 1);
nGroup   = zeros(k, 1);
for i = 1:k
    sel         = groups == cats{i};
    nGroup(i)   = nnz(sel);
    meanRank(i) = mean(ranks(sel));
end

rows = {};
for i = 1:k-1
    for j = i+1:k
        se = sqrt((N * (N + 1) / 12 - tieTerm) * (1/nGroup(i) + 1/nGroup(j)));
        z  = (meanRank(i) - meanRank(j)) / se;
        p  = 2 * (1 - normcdf(abs(z)));

        % Rank-biserial correlation from the Mann-Whitney U of the pair.
        a = values(groups == cats{i});
        b = values(groups == cats{j});
        U = mann_whitney_u(a, b);
        rb = 1 - 2 * U / (numel(a) * numel(b));

        r.groupA    = string(cats{i});
        r.groupB    = string(cats{j});
        r.nA        = nGroup(i);
        r.nB        = nGroup(j);
        r.medianA   = median(a);
        r.medianB   = median(b);
        r.z         = z;
        r.p         = p;
        r.rankBiserial = rb;
        rows{end+1} = r; %#ok<AGROW>
    end
end

R = struct2table([rows{:}]);
R.pAdjusted = stats.adjust_p(R.p, correction);
end

% -------------------------------------------------------------------------
function U = mann_whitney_u(a, b)
%MANN_WHITNEY_U U statistic for group a against group b.
combined = [a(:); b(:)];
r        = tiedrank(combined);
Ra       = sum(r(1:numel(a)));
U        = Ra - numel(a) * (numel(a) + 1) / 2;
end
