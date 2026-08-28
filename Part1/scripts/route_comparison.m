function T = route_comparison()
%ROUTE_COMPARISON Side-by-side between-grade separation for both routes.
%
%   T = ROUTE_COMPARISON() reports, for every descriptor, the Kruskal-Wallis
%   p, effect size and benign-vs-melanoma AUC computed independently on the
%   expert-mask route and on the automatic route.
%
%   This exists to answer one question that the per-route analyses cannot:
%   whether a descriptor's between-grade signal depends on which
%   segmentation produced it. That matters most for D_f. The ISIC expert
%   masks are largely polygon-traced, and a polygon boundary is locally
%   straight by construction, so its box dimension is pulled toward 1
%   regardless of the lesion underneath. If D_f separates grades on the
%   automatic route but not on the expert route, the flat expert result is an
%   annotation artefact rather than a biological null -- and the primary
%   route for D_f specifically should be the automatic one, even though the
%   expert masks remain the better reference for perimeter-derived
%   descriptors.
%
%   Each route is analysed on the lesions that route accepted, so the two
%   columns are not paired. That is deliberate: the question is what each
%   route would conclude on its own, which is how each would be used.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root, 'config'));
cfg = default_config();

R = readtable(fullfile(cfg.paths.results, 'descriptors.csv'), 'TextType', 'string');
R.grade = categorical(R.grade, cfg.data.gradeOrder, 'Ordinal', true);

descriptors = cfg.stats.descriptors;
rows = {};

for d = 1:numel(descriptors)
    name = descriptors{d};

    % Pre-seed every field so struct2table still works when a route produced
    % too few accepted lesions to analyse; a missing route should leave a
    % visible NaN rather than abort the whole comparison.
    r = struct('descriptor', string(name));
    for route = {'expert', 'otsu'}
        for f = {'n','medBen','medInv','p','eta2','auc'}
            r.([route{1} '_' f{1}]) = NaN;
        end
    end

    for route = {'expert', 'otsu'}
        rt  = route{1};
        col = [rt '_' name];
        acc = [rt '_accepted'];

        if ~all(ismember({col, acc}, R.Properties.VariableNames))
            continue
        end

        S  = R(logical(R.(acc)), :);
        v  = S.(col);
        g  = S.grade;
        ok = isfinite(v);
        v  = v(ok); g = g(ok);

        if numel(v) < 12 || numel(categories(removecats(g))) < 2
            continue
        end

        t = stats.group_tests(v, g, cfg);
        p = stats.discriminative_power(v, g, cfg);

        r.([rt '_n'])      = numel(v);
        r.([rt '_medBen']) = median(v(g == cfg.data.gradeOrder{1}), 'omitnan');
        r.([rt '_medInv']) = median(v(g == cfg.data.gradeOrder{end}), 'omitnan');
        r.([rt '_p'])      = t.p;
        r.([rt '_eta2'])   = t.etaSqH;
        r.([rt '_auc'])    = p.auc;
    end

    rows{end+1} = r; %#ok<AGROW>
end

T = struct2table([rows{:}]);

fprintf('\n%-20s | %-28s | %-28s\n', '', 'EXPERT route', 'OTSU route');
fprintf('%-20s | %5s %8s %7s %6s | %5s %8s %7s %6s\n', ...
    'descriptor', 'n', 'KW p', 'eta2_H', 'AUC', 'n', 'KW p', 'eta2_H', 'AUC');
fprintf('%s\n', repmat('-', 1, 84));
for k = 1:height(T)
    fprintf('%-20s | %5.0f %8.4g %7.4f %6.3f | %5.0f %8.4g %7.4f %6.3f\n', ...
        T.descriptor(k), ...
        T.expert_n(k), T.expert_p(k), T.expert_eta2(k), T.expert_auc(k), ...
        T.otsu_n(k),   T.otsu_p(k),   T.otsu_eta2(k),   T.otsu_auc(k));
end

fprintf('\nbenign -> invasive median shift (sign shows direction):\n');
fprintf('%-20s %12s %12s\n', 'descriptor', 'expert', 'otsu');
for k = 1:height(T)
    fprintf('%-20s %12.4f %12.4f\n', T.descriptor(k), ...
        T.expert_medInv(k) - T.expert_medBen(k), ...
        T.otsu_medInv(k)   - T.otsu_medBen(k));
end

writetable(T, fullfile(cfg.paths.results, 'route_comparison.csv'));
end
