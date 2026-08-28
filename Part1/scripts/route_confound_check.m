function out = route_confound_check()
%ROUTE_CONFOUND_CHECK Is the automatic route's grade separation biological?
%
%   The automatic route separates grades on nearly every descriptor while the
%   expert route separates almost none. There are two incompatible readings:
%
%     (A) Real geometry. Expert masks are polygon-traced, which erases the
%         fine boundary structure the descriptors measure, so the expert
%         route is blind and the automatic route sees the biology.
%
%     (B) Segmentation artefact. Ill-defined infiltrative borders are harder
%         to threshold, so the automatic route fails progressively worse with
%         severity. The "signal" is then a gradient of segmentation quality
%         wearing the costume of biology.
%
%   (B) is the more dangerous reading because it produces exactly the result
%   the study hopes to see. These checks discriminate between them:
%
%     1. Does route agreement (Dice) fall with grade? Monotone decline is
%        the signature of (B) -- tested with Jonckheere-Terpstra.
%     2. Does automatic-route acceptance fall with grade? Differential
%        dropout means the surviving malignant lesions are a selected,
%        easier subset.
%     3. Does the separation survive on well-segmented lesions only? This is
%        the decisive one. Under (A) the effect persists once poorly
%        segmented cases are removed; under (B) it collapses.
%
%   Nothing here is conclusive alone, but (3) failing while (1) shows a
%   monotone decline would be strong evidence against promoting the
%   automatic route.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root, 'config'));
cfg = default_config();

R = readtable(fullfile(cfg.paths.results, 'descriptors.csv'), 'TextType', 'string');
R.grade = categorical(R.grade, cfg.data.gradeOrder, 'Ordinal', true);

out = struct();

%% ---- 1. Acceptance rate by grade -------------------------------------
fprintf('===== 1. route acceptance by grade =====\n');
fprintf('%-10s %8s %10s %10s\n', 'grade', 'n', 'expert', 'otsu');
for g = cfg.data.gradeOrder
    sel = R.grade == g{1};
    fprintf('%-10s %8d %9.1f%% %9.1f%%\n', g{1}, nnz(sel), ...
        100 * nnz(sel & R.expert_accepted) / nnz(sel), ...
        100 * nnz(sel & R.otsu_accepted)   / nnz(sel));
end

%% ---- 2. Agreement by grade --------------------------------------------
both = R(logical(R.expert_accepted) & logical(R.otsu_accepted) & isfinite(R.dice), :);
fprintf('\n===== 2. expert-vs-otsu agreement by grade =====\n');
fprintf('%-10s %6s %10s %10s %12s\n', 'grade', 'n', 'medDice', 'meanDice', 'medAreaRatio');
for g = cfg.data.gradeOrder
    sel = both.grade == g{1};
    if ~any(sel), continue; end
    fprintf('%-10s %6d %10.3f %10.3f %12.3f\n', g{1}, nnz(sel), ...
        median(both.dice(sel)), mean(both.dice(sel)), ...
        median(both.areaRatio(sel), 'omitnan'));
end

% Kruskal-Wallis alongside the trend test. They answer different questions
% and the pilot showed why both are needed: agreement differed sharply by
% grade while the MONOTONE component was weak, because in-situ segmented
% worst of the three. A trend test alone would report "no problem" for a
% pattern that still leaves one grade's descriptors resting on poorer masks.
kwDice = stats.group_tests(both.dice, both.grade, cfg);
fprintf('\nKruskal-Wallis on Dice across grades: p = %.4g, eta2_H = %.4f\n', ...
    kwDice.p, kwDice.etaSqH);
if kwDice.p < cfg.stats.alpha
    fprintf(['  -> agreement differs by grade even if not monotonically;\n' ...
             '     check which grade is worst, not just the direction.\n']);
end
out.kwDice = kwDice;

[jt.p, jt.z, jt.trend] = jonckheere_terpstra(both.dice, both.grade, cfg.data.gradeOrder);
fprintf('\nJonckheere-Terpstra trend in Dice across grades: z = %.3f, p = %.4g (%s)\n', ...
    jt.z, jt.p, jt.trend);
if jt.p < cfg.stats.alpha && jt.z < 0
    fprintf(['  -> agreement DEGRADES monotonically with severity. The automatic\n' ...
             '     route''s grade separation is confounded with segmentation quality.\n']);
else
    fprintf(['  -> no significant monotone decline; segmentation quality is not\n' ...
             '     obviously tracking severity.\n']);
end
out.jt = jt;

%% ---- 2b. Hair contamination of D_f ------------------------------------
% Morphological opening is disabled in the automatic route because it eroded
% the lesion periphery. That leaves hair suppression to DullRazor alone, and
% residual hair is thin, dark and boundary-attached -- the exact structure
% that inflates a box dimension. If D_f rises with the fraction of pixels
% flagged as hair, the automatic route is measuring hair, not lesion.
if ismember('otsu_hairFraction', R.Properties.VariableNames)
    H = R(logical(R.otsu_accepted) & isfinite(R.otsu_hairFraction) & ...
          isfinite(R.otsu_Df), :);
    if height(H) >= 12
        [rho, pRho] = corr(H.otsu_hairFraction, H.otsu_Df, 'Type', 'Spearman');
        fprintf('\n===== 2b. hair fraction vs D_f (automatic route) =====\n');
        fprintf('Spearman rho = %.3f, p = %.4g, n = %d\n', rho, pRho, height(H));
        if pRho < cfg.stats.alpha && rho > 0
            fprintf(['  -> D_f rises with residual hair. Re-enable a small\n' ...
                     '     opening or strengthen hair removal before trusting D_f.\n']);
        else
            fprintf('  -> no significant hair-driven inflation of D_f.\n');
        end

        % Hair burden must also not track grade, or it would masquerade as
        % between-grade signal even without inflating D_f overall.
        kwHair = stats.group_tests(H.otsu_hairFraction, H.grade, cfg);
        fprintf('hair fraction differs by grade: p = %.4g\n', kwHair.p);
        out.hairRho = rho; out.hairP = pRho; out.hairByGrade = kwHair.p;

        % Whether hair actually DRIVES the between-grade D_f difference is a
        % separate question from whether it inflates D_f overall. Restricting
        % to the least hairy lesions answers it: if the grade difference
        % survives where there is almost no hair to confound it, the effect
        % is not a hair artefact.
        cut  = quantile(H.otsu_hairFraction, cfg.qc.lowHairQuantile);
        L    = H(H.otsu_hairFraction <= cut, :);
        fprintf('\nrestricted to the least-hairy %.0f%% (hairFrac <= %.5f, n = %d):\n', ...
            100 * cfg.qc.lowHairQuantile, cut, height(L));
        lc = arrayfun(@(g) nnz(L.grade == g), categories(L.grade));
        fprintf('  per grade: %s\n', strjoin(string(lc), ' / '));

        for nm = ["Df", "compactness", "solidity"]
            c = "otsu_" + nm;
            if ~ismember(c, L.Properties.VariableNames), continue; end
            aAll = summarise(H.(c), H.grade, cfg);
            aLow = summarise(L.(c), L.grade, cfg);
            fprintf('  %-12s all: p=%-10.4g eta2=%.4f  |  low-hair: p=%-10.4g eta2=%.4f\n', ...
                nm, aAll.p, aAll.eta2, aLow.p, aLow.eta2);
        end
    end
end

%% ---- 3. Does separation survive on well-segmented lesions? ------------
% Restricted to lesions where both routes substantially agree, so the
% automatic mask is known to be reasonable. If the effect is geometric it
% should persist here; if it was segmentation failure it should vanish.
thr  = cfg.qc.wellSegmentedDice;
good = both(both.dice >= thr, :);

fprintf('\n===== 3. separation restricted to well-segmented lesions (Dice >= %.2f) =====\n', thr);
fprintf('n = %d of %d paired lesions', height(good), height(both));
gc = arrayfun(@(g) nnz(good.grade == g), categories(good.grade));
fprintf('  (per grade: %s)\n\n', strjoin(string(gc), ' / '));

fprintf('%-20s | %-22s | %-22s\n', '', 'ALL paired', 'well-segmented only');
fprintf('%-20s | %8s %6s %6s | %8s %6s %6s\n', ...
    'descriptor', 'KW p', 'eta2', 'AUC', 'KW p', 'eta2', 'AUC');
fprintf('%s\n', repmat('-', 1, 74));

rows = {};
for d = 1:numel(cfg.stats.descriptors)
    name = cfg.stats.descriptors{d};
    col  = ['otsu_' name];
    if ~ismember(col, both.Properties.VariableNames), continue; end

    a = summarise(both.(col), both.grade, cfg);
    b = summarise(good.(col), good.grade, cfg);

    fprintf('%-20s | %8.4g %6.3f %6.3f | %8.4g %6.3f %6.3f\n', ...
        name, a.p, a.eta2, a.auc, b.p, b.eta2, b.auc);

    rows{end+1} = struct('descriptor', string(name), ...
        'allP', a.p, 'allEta2', a.eta2, 'allAuc', a.auc, ...
        'goodP', b.p, 'goodEta2', b.eta2, 'goodAuc', b.auc); %#ok<AGROW>
end
out.table = struct2table([rows{:}]);

nSurv = nnz(out.table.goodP < cfg.stats.alpha);
fprintf(['\n%d of %d descriptors remain significant among well-segmented\n' ...
         'lesions.\n'], nSurv, height(out.table));

writetable(out.table, fullfile(cfg.paths.results, 'route_confound_check.csv'));
end

% =========================================================================
function s = summarise(v, g, cfg)
ok = isfinite(v);
v = v(ok); g = removecats(g(ok));
s = struct('p', NaN, 'eta2', NaN, 'auc', NaN);
if numel(v) < 12 || numel(categories(g)) < 2
    return
end
t = stats.group_tests(v, g, cfg);
p = stats.discriminative_power(v, g, cfg);
s.p = t.p; s.eta2 = t.etaSqH; s.auc = p.auc;
end

% =========================================================================
function [p, z, trend] = jonckheere_terpstra(x, g, order)
%JONCKHEERE_TERPSTRA Test for a monotone trend across ordered groups.
%
%   Unlike Kruskal-Wallis, which is sensitive to any difference, this is
%   directed: it asks specifically whether values increase (or decrease)
%   along the grade ordering. That is the right question for "does agreement
%   degrade with severity", where an unordered difference would not support
%   the confounding explanation.

x = x(:);
groups = cell(1, numel(order));
for k = 1:numel(order)
    groups{k} = x(g == order{k} & isfinite(x));
end
groups = groups(~cellfun(@isempty, groups));
K = numel(groups);

if K < 2
    [p, z] = deal(NaN); trend = 'undetermined'; return
end

% U statistic: count of concordant pairs across every ordered group pair.
U = 0;
for i = 1:K-1
    for j = i+1:K
        a = groups{i}; b = groups{j};
        U = U + sum(sum(b(:)' > a(:))) + 0.5 * sum(sum(b(:)' == a(:)));
    end
end

n  = cellfun(@numel, groups);
N  = sum(n);
mu = (N^2 - sum(n.^2)) / 4;
sd = sqrt((N^2 * (2*N + 3) - sum(n.^2 .* (2*n + 3))) / 72);

if sd == 0
    [p, z] = deal(NaN); trend = 'undetermined'; return
end

z = (U - mu) / sd;
p = 2 * (1 - normcdf(abs(z)));

if z > 0
    trend = 'increasing with grade';
else
    trend = 'decreasing with grade';
end
end
