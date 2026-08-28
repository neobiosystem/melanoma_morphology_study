function make_figures(route)
%MAKE_FIGURES Publication figures for Part 1.
%
%   Writes vector PDFs plus PNG previews into results/figures:
%       fig1_descriptor_distributions - per-grade box plots, n annotated
%       fig2_boxcount_fit             - log-log fit with the scaling region
%       fig3_roc                      - ROC per descriptor
%       fig4_bland_altman             - expert vs Otsu agreement
%       fig5_agreement_by_grade       - does agreement degrade with severity
%       fig6_perimeter_scale_stability - compactness gap across smoothing scales
%
%   Figure 6 is the one that defends the headline result: if the between
%   grade compactness gap survives across smoothing scales it is geometric,
%   and if it vanishes once sigma exceeds a few pixels it was a rasterisation
%   artefact.

if nargin < 1
    route = 'expert';
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root, 'config'));
cfg = default_config();

if ~isfolder(cfg.paths.figures)
    mkdir(cfg.paths.figures);
end

R = readtable(fullfile(cfg.paths.results, 'descriptors.csv'), 'TextType', 'string');
R.grade = categorical(R.grade, cfg.data.gradeOrder, 'Ordinal', true);
A = R(logical(R.([route '_accepted'])), :);

labels = cfg.data.gradeLabels;
descriptors = cfg.stats.descriptors;

% ---- Fig 1: distributions ---------------------------------------------
f = figure('Visible', 'off', 'Position', [0 0 1500 850]);
tiledlayout(2, ceil(numel(descriptors)/2), 'Padding', 'compact');
for d = 1:numel(descriptors)
    col = [route '_' descriptors{d}];
    if ~ismember(col, A.Properties.VariableNames), continue; end
    nexttile;
    % boxchart rather than boxplot: boxplot does not lay out correctly
    % inside a tiledlayout (MATLAB warns), and these are publication figures.
    boxchart(A.grade, A.(col));
    xticklabels(labels);       % clinical names, not the internal category keys
    title(descriptors{d}, 'Interpreter', 'none');
    grid on;

    % Several descriptors are heavy-tailed (compactness reaches ~14 on a few
    % ragged lesions). Autoscaling to the extremes squashes every box into an
    % unreadable sliver, hiding the between-grade differences the figure
    % exists to show. Clip the view to the bulk; the outliers are still in
    % the data and in the statistics, just not dictating the axis.
    v  = A.(col); v = v(isfinite(v));
    if ~isempty(v)
        lo = quantile(v, 0.01); hi = quantile(v, 0.99);
        if hi > lo
            pad = 0.05 * (hi - lo);
            ylim([lo - pad, hi + pad]);
        end
    end

    counts = arrayfun(@(g) nnz(A.grade == g), categories(A.grade));
    xlabel(sprintf('n = %s', strjoin(string(counts), ' / ')));
end
sgtitle(sprintf('Morphological descriptors by grade (%s route, n = %d)', ...
    route, height(A)));
save_fig(f, cfg, ['fig1_descriptor_distributions_' route]);

% ---- Fig 2: box-count fit for one representative lesion ---------------
f = figure('Visible', 'off', 'Position', [0 0 900 400]);
example_boxcount_figure(cfg, A, route);
save_fig(f, cfg, ['fig2_boxcount_fit_' route]);

% ---- Fig 3: ROC --------------------------------------------------------
f = figure('Visible', 'off', 'Position', [0 0 700 620]);
hold on;
isMal = A.grade ~= cfg.data.gradeOrder{1};
legendEntries = {};
for d = 1:numel(descriptors)
    col = [route '_' descriptors{d}];
    if ~ismember(col, A.Properties.VariableNames), continue; end
    v = A.(col);
    ok = isfinite(v);
    if nnz(ok) < 10, continue; end
    dir = 1;
    if median(v(ok & isMal)) < median(v(ok & ~isMal)), dir = -1; end
    [x, y, ~, auc] = perfcurve(double(isMal(ok)), dir * v(ok), 1);
    plot(x, y, 'LineWidth', 1.4);
    legendEntries{end+1} = sprintf('%s (AUC %.3f)', descriptors{d}, auc); %#ok<AGROW>
end
plot([0 1], [0 1], 'k--');
legendEntries{end+1} = 'chance';
xlabel('1 - specificity'); ylabel('sensitivity');
title('Benign vs melanoma (in situ + invasive)');
legend(legendEntries, 'Location', 'southeast', 'Interpreter', 'none');
grid on; axis square;
save_fig(f, cfg, ['fig3_roc_' route]);

% ---- Fig 4 and 5: dual-route agreement ---------------------------------
if ismember('dice', R.Properties.VariableNames) && any(isfinite(R.dice))
    both = R(logical(R.expert_accepted) & logical(R.otsu_accepted), :);

    f = figure('Visible', 'off', 'Position', [0 0 1400 800]);
    tiledlayout(2, 3, 'Padding', 'compact');
    for d = 1:min(6, numel(descriptors))
        ec = ['expert_' descriptors{d}];
        oc = ['otsu_'   descriptors{d}];
        if ~all(ismember({ec, oc}, both.Properties.VariableNames)), continue; end
        x = both.(ec); y = both.(oc);
        ok = isfinite(x) & isfinite(y);
        if nnz(ok) < 5, continue; end

        nexttile;
        m = (x(ok) + y(ok)) / 2;
        df = y(ok) - x(ok);
        scatter(m, df, 6, 'filled', 'MarkerFaceAlpha', 0.3);
        yline(mean(df), 'r-', 'LineWidth', 1.2);
        yline(mean(df) + 1.96*std(df), 'r--');
        yline(mean(df) - 1.96*std(df), 'r--');
        title(descriptors{d}, 'Interpreter', 'none');
        xlabel('mean of routes'); ylabel('otsu - expert'); grid on;
    end
    sgtitle('Bland-Altman: automatic vs expert route');
    save_fig(f, cfg, 'fig4_bland_altman');

    f = figure('Visible', 'off', 'Position', [0 0 900 420]);
    tiledlayout(1, 2, 'Padding', 'compact');
    nexttile;
    boxchart(both.grade, both.dice);
    xticklabels(labels);
    ylabel('Dice (expert vs Otsu)'); grid on;
    title('Segmentation agreement by grade');
    nexttile;
    boxchart(both.grade, both.areaRatio);
    xticklabels(labels);
    yline(1, 'k--');
    ylabel('area ratio (Otsu / expert)'); grid on;
    title('Over/under-segmentation by grade');
    save_fig(f, cfg, 'fig5_agreement_by_grade');
end

% ---- Fig 6: perimeter-bias scale stability ----------------------------
f = figure('Visible', 'off', 'Position', [0 0 800 450]);
rawCol = [route '_compactnessRaw'];
smCol  = [route '_compactness'];
if all(ismember({rawCol, smCol}, A.Properties.VariableNames))
    g = categories(A.grade);
    rawMed = arrayfun(@(k) median(A.(rawCol)(A.grade == g{k}), 'omitnan'), 1:numel(g));
    smMed  = arrayfun(@(k) median(A.(smCol)(A.grade == g{k}), 'omitnan'), 1:numel(g));
    bar([rawMed(:), smMed(:)]);
    set(gca, 'XTickLabel', labels);
    legend({'raw raster perimeter', 'fixed-scale smoothed'}, 'Location', 'best');
    ylabel('median compactness'); grid on;
    title({'Perimeter-bias check', ...
           'a gap that survives smoothing is geometric, not rasterisation'});
end
save_fig(f, cfg, ['fig6_perimeter_scale_stability_' route]);

fprintf('Figures written to %s\n', cfg.paths.figures);
end

% =========================================================================
function example_boxcount_figure(cfg, A, route)
%EXAMPLE_BOXCOUNT_FIGURE Recompute one lesion's curve to show the fit.
idx = find(A.([route '_DfFlag']) == "ok", 1);
if isempty(idx)
    idx = 1;
end
% The mask must come from the route being plotted, or an "otsu" figure would
% actually display an expert-mask fit under an Otsu label.
if strcmp(route, 'otsu')
    imgPath = fullfile(cfg.paths.images, char(A.isic_id(idx) + '.jpg'));
    if ~isfile(imgPath)
        return
    end
    mask = seg.otsu_route(imread(imgPath), cfg);
else
    maskPath = fullfile(cfg.paths.masks, char(A.isic_id(idx) + cfg.data.maskSuffix));
    if ~isfile(maskPath)
        return
    end
    mask = seg.expert_route(maskPath);
end
out  = morph.process_lesion(mask, cfg);
b    = out.boxcount;
if isempty(b) || ~isfield(b, 'epsAll') || isempty(b.epsAll)
    return
end

x = log(1 ./ b.epsAll);
y = log(b.countsAll);

subplot(1,2,1);
imshow(mask); title(sprintf('%s (%s)', A.isic_id(idx), A.grade(idx)));

subplot(1,2,2);
plot(x, y, 'o', 'MarkerSize', 4); hold on;
if ~isempty(b.fitIdx)
    xf = x(b.fitIdx); yf = y(b.fitIdx);
    plot(xf, yf, 'r.', 'MarkerSize', 14);
    pf = polyfit(xf, yf, 1);
    plot(xf, polyval(pf, xf), 'r-', 'LineWidth', 1.5);
end
xlabel('log(1/\epsilon)'); ylabel('log N(\epsilon)');
title(sprintf('D_f = %.4f (R^2 = %.4f)', out.Df, out.DfR2));
legend({'all scales', 'scaling region', 'fit'}, 'Location', 'northwest');
grid on;
end

% =========================================================================
function save_fig(f, cfg, name)
exportgraphics(f, fullfile(cfg.paths.figures, [name '.pdf']), 'ContentType', 'vector');
exportgraphics(f, fullfile(cfg.paths.figures, [name '.png']), 'Resolution', 150);
close(f);
end
