function padj = adjust_p(p, method)
%ADJUST_P Multiple-comparison correction.
%
%   padj = ADJUST_P(p, 'holm') applies Holm-Bonferroni, controlling the
%   family-wise error rate. Used WITHIN a descriptor across its three
%   pairwise post-hoc comparisons, where any false positive would be claimed
%   as a specific between-grade difference.
%
%   padj = ADJUST_P(p, 'bh') applies Benjamini-Hochberg, controlling the
%   false discovery rate. Used ACROSS descriptors on the omnibus tests,
%   where the question is which descriptors are worth pursuing and FWER
%   would be needlessly conservative.

p = p(:);
n = numel(p);
padj = nan(n, 1);

if n == 0
    return
end

[sorted, order] = sort(p);

switch lower(method)
    case 'holm'
        adj = sorted .* (n - (1:n)' + 1);
        adj = cummax(adj);              % enforce monotonicity
    case {'bh', 'fdr'}
        adj = sorted .* n ./ (1:n)';
        adj = flipud(cummin(flipud(adj)));
    otherwise
        error('stats:adjust_p:badMethod', 'Unknown method: %s', method);
end

padj(order) = min(adj, 1);
end
