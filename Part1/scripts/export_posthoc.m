function T = export_posthoc()
%EXPORT_POSTHOC Tidy per-pair post-hoc table for both routes.
%
%   Writes results/posthoc_all.csv with one row per (route, descriptor,
%   group pair): medians, Holm-corrected p, and rank-biserial effect size.
%
%   Direction claims ("melanoma measures less compact than nevi on expert
%   masks") are pairwise claims, but the omnibus Kruskal-Wallis p in
%   route_comparison.csv is a three-group statistic and cannot support them.
%   Reporting a direction next to an omnibus p invites the reader to assume
%   the direction was tested; this file is what makes it actually testable.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root, 'config'));
cfg = default_config();

routes = {'expert', 'otsu'};
T = table();

for r = 1:numel(routes)
    f = fullfile(cfg.paths.results, ['analysis_' routes{r} '.mat']);
    if ~isfile(f)
        fprintf('missing %s, skipping\n', f);
        continue
    end
    S  = load(f);
    A  = S.A;
    ds = fieldnames(A.tests);

    for i = 1:numel(ds)
        t  = A.tests.(ds{i});
        ph = t.postHoc;
        if isempty(ph)
            continue
        end
        n   = height(ph);
        add = table();
        add.route        = repmat(string(routes{r}), n, 1);
        add.descriptor   = repmat(string(ds{i}),     n, 1);
        add.etaSqH       = repmat(t.etaSqH,          n, 1);
        add.kwP          = repmat(t.p,               n, 1);
        add.groupA       = string(ph.groupA);
        add.groupB       = string(ph.groupB);
        add.medianA      = ph.medianA;
        add.medianB      = ph.medianB;
        add.pAdjusted    = ph.pAdjusted;
        add.rankBiserial = ph.rankBiserial;
        T = [T; add]; %#ok<AGROW>
    end
end

outPath = fullfile(cfg.paths.results, 'posthoc_all.csv');
writetable(T, outPath);
fprintf('Wrote %s (%d rows)\n', outPath, height(T));
end
