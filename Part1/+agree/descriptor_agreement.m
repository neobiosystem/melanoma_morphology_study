function out = descriptor_agreement(x, y, groups)
%DESCRIPTOR_AGREEMENT Paired agreement between two measurement routes.
%
%   out = DESCRIPTOR_AGREEMENT(x, y) compares descriptor values measured on
%   the same lesions by the expert route (x) and the automatic route (y).
%   out = DESCRIPTOR_AGREEMENT(x, y, groups) additionally reports the same
%   statistics stratified by grade.
%
%   Fields:
%       ccc              - Lin's concordance correlation coefficient
%       bias, loA        - Bland-Altman mean difference and 95% limits
%       spearman         - rank correlation
%       reliabilityRatio - between-lesion SD divided by within-lesion route
%                          difference SD
%
%   reliabilityRatio is the number that matters for Part 3. It says how much
%   of the biological spread across lesions survives the noise introduced by
%   choosing a segmentation route. A ratio near 1 means route choice is as
%   large as the biology, and the descriptor cannot serve as a calibration
%   target without an expert mask.
%
%   CCC is primary rather than Pearson r because r is blind to bias: two
%   routes can correlate perfectly while one reads systematically high.

if nargin < 3
    groups = [];
end

out = struct();

x = x(:);
y = y(:);
ok = isfinite(x) & isfinite(y);
x  = x(ok);
y  = y(ok);

out.n = numel(x);
if out.n < 3
    [out.ccc, out.bias, out.spearman, out.reliabilityRatio] = deal(NaN);
    out.loA = [NaN NaN];
    out.byGroup = table();
    return
end

out.ccc      = lin_ccc(x, y);
d            = y - x;
out.bias     = mean(d);
out.sdDiff   = std(d);
out.loA      = out.bias + [-1.96, 1.96] * out.sdDiff;
out.spearman = corr(x, y, 'Type', 'Spearman');

betweenSd = std([x; y]);
if out.sdDiff > 0
    out.reliabilityRatio = betweenSd / out.sdDiff;
else
    out.reliabilityRatio = Inf;
end

% ---- Stratified by grade ----------------------------------------------
% Whether agreement degrades with severity is a finding in its own right,
% not a footnote: if ill-defined infiltrative borders are objectively harder
% to segment, that supports the biological narrative AND warns that an
% automatic-only analysis would confuse segmentation difficulty with
% biology.
out.byGroup = table();
if ~isempty(groups)
    g = groups(ok);
    cats = categories(removecats(categorical(g)));
    rows = {};
    for k = 1:numel(cats)
        sel = categorical(g) == cats{k};
        if nnz(sel) < 3
            continue
        end
        r.grade     = string(cats{k});
        r.n         = nnz(sel);
        r.ccc       = lin_ccc(x(sel), y(sel));
        r.bias      = mean(y(sel) - x(sel));
        r.sdDiff    = std(y(sel) - x(sel));
        r.medianAbs = median(abs(y(sel) - x(sel)));
        rows{end+1} = r; %#ok<AGROW>
    end
    if ~isempty(rows)
        out.byGroup = struct2table([rows{:}]);
    end
end
end

% -------------------------------------------------------------------------
function c = lin_ccc(x, y)
%LIN_CCC Lin's concordance correlation coefficient.
mx = mean(x);  my = mean(y);
vx = var(x, 1); vy = var(y, 1);
sxy = mean((x - mx) .* (y - my));
denom = vx + vy + (mx - my)^2;
if denom == 0
    c = NaN;
else
    c = 2 * sxy / denom;
end
end
