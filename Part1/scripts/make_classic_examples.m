function T = make_classic_examples(k)
%MAKE_CLASSIC_EXAMPLES Pick and render the k most representative lesions per grade.
%
%   T = MAKE_CLASSIC_EXAMPLES()  picks 3 per grade (9 total).
%   T = MAKE_CLASSIC_EXAMPLES(k) picks k per grade.
%
%   Writes into classical_examples/:
%       <grade>_<rank>_<id>_panel.png   four-panel exemplar sheet
%       <grade>_<rank>_<id>_image.jpg   source dermoscopic image
%       <grade>_<rank>_<id>_mask.png    primary-route (Otsu) mask
%       manifest.csv                    values, group medians, license, attribution
%       README.md                       selection rule and licensing
%
%   "Typical" is defined operationally, because picking exemplars by eye is
%   how a figure ends up arguing for a conclusion the cohort does not support.
%   A lesion is typical when it sits close to its own grade's median across
%   ALL seven descriptors simultaneously, measured on the primary (automatic)
%   route.
%
%   Three choices in that rule carry the weight:
%
%   * Distance uses the MEDIAN of |z| across descriptors, not the mean or the
%     Euclidean norm. A lesion that is dead-central on six descriptors and
%     unusual on one is still a fair representative; a sum-based distance
%     would disqualify it while admitting a lesion that is mildly off on
%     everything.
%   * z-scores are formed per grade from that grade's own median and IQR, so
%     "typical" means typical OF ITS GRADE. Scoring against the pooled cohort
%     would systematically pick the most benign-looking melanomas.
%   * Candidates are restricted to well-segmented lesions (Dice >= the
%     configured threshold), but the reference median is taken from the FULL
%     accepted grade population. An exemplar is a claim that this is what the
%     measurement saw, so its mask must be right; but it is also a claim to
%     represent the published interval, so it must be scored against that
%     interval rather than against the filtered pool. The two populations are
%     not interchangeable -- the Dice gate keeps smoother lesions, since
%     irregular borders are the ones the routes disagree on.
%
%   Among the most-typical pool the final k are spread out by a greedy
%   max-min rule rather than taken as the top k. The top k by distance are
%   frequently near-duplicates -- same site, same morphology -- which looks
%   like corroboration while carrying the information of a single case.

if nargin < 1
    k = 3;
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root, 'config'));
cfg = default_config();

outDir = cfg.paths.classicExamples;
if ~isfolder(outDir)
    mkdir(outDir);
end

R = readtable(fullfile(cfg.paths.results, 'descriptors.csv'), 'TextType', 'string');
R.grade = categorical(R.grade, cfg.data.gradeOrder, 'Ordinal', true);

% readtable gives logicals back as numeric and blanks as <missing>; both break
% logical indexing and string comparison respectively.
for v = ["otsu_accepted", "expert_accepted"]
    R.(v) = logical(R.(v));
end
for v = ["otsu_DfFlag", "otsu_flags", "expert_flags"]
    if ismember(v, R.Properties.VariableNames)
        R.(v)(ismissing(R.(v))) = "";
    end
end

% Attribution is required for reuse and is not carried in descriptors.csv.
R.attribution = strings(height(R), 1);
if isfile(cfg.paths.metadata)
    M = readtable(cfg.paths.metadata, 'TextType', 'string');
    if ismember('attribution', M.Properties.VariableNames)
        M.attribution(ismissing(M.attribution)) = "";
        [tf, loc] = ismember(R.isic_id, M.isic_id);
        R.attribution(tf) = M.attribution(loc(tf));
    end
end

descriptors = cfg.stats.descriptors;
V = nan(height(R), numel(descriptors));
for d = 1:numel(descriptors)
    V(:, d) = R.(['otsu_' descriptors{d}]);
end

eligible = R.otsu_accepted & R.expert_accepted & ...
           R.otsu_DfFlag == "ok" & ...
           isfinite(R.dice) & R.dice >= cfg.qc.wellSegmentedDice & ...
           all(isfinite(V), 2);

fprintf('Eligible for exemplar selection: %d / %d lesions\n', ...
    nnz(eligible), height(R));

rows = table();
for gi = 1:numel(cfg.data.gradeOrder)
    gname = cfg.data.gradeOrder{gi};
    idx   = find(eligible & R.grade == gname);
    if numel(idx) < k
        fprintf('  %s: only %d eligible, skipping\n', gname, numel(idx));
        continue
    end

    % The reference median comes from the FULL accepted grade population, not
    % from the well-segmented candidates. Those two differ substantially: the
    % Dice gate preferentially keeps smoother lesions, because irregular
    % borders are exactly the ones the two routes disagree on. Scoring against
    % the filtered subset would yield exemplars that are typical of an
    % artefact of the filter while visibly contradicting the published
    % intervals -- invasive compactness would show ~1.35 against a reported
    % grade median of ~1.81. Scoring against the full population instead picks
    % the best available representatives of the numbers actually reported, and
    % any residual gap stays visible in each panel's z-profile.
    pop = R.otsu_accepted & R.grade == gname & all(isfinite(V), 2);
    med = median(V(pop, :), 1);
    sc  = iqr(V(pop, :), 1) / 1.349;   % robust sd
    sc(~isfinite(sc) | sc <= 0) = 1;   % a degenerate descriptor must not dominate

    Z = (V(idx, :) - med) ./ sc;

    dist     = median(abs(Z), 2);
    [~, ord] = sort(dist);
    poolN    = min(numel(idx), max(12, round(0.05 * numel(idx))));
    pool     = ord(1:poolN);
    pick     = greedy_diverse(Z(pool, :), k);

    fprintf('  %s: %d in grade, %d well-segmented candidates, pool %d, picked %d\n', ...
        gname, nnz(pop), numel(idx), poolN, numel(pick));

    for r = 1:numel(pick)
        p   = pool(pick(r));
        sel = idx(p);
        try
            render_example(R(sel, :), cfg, gi, r, descriptors, ...
                Z(p, :), med, outDir);
        catch ME
            fprintf('    render failed for %s: %s\n', R.isic_id(sel), ME.message);
            continue
        end
        rows = [rows; summarise(R(sel, :), cfg, gi, r, descriptors, ...
            Z(p, :), med, dist(p))]; %#ok<AGROW>
    end
end

T = rows;
if ~isempty(T)
    writetable(T, fullfile(outDir, 'manifest.csv'));
end
write_readme(outDir, cfg, k, descriptors);
fprintf('Exemplars written to %s\n', outDir);
end

% =========================================================================
function pick = greedy_diverse(Z, k)
%GREEDY_DIVERSE Spread k choices across a pool already filtered for typicality.
%
%   Z is ordered most-typical-first, so index 1 seeds the set; each further
%   pick maximises its minimum distance to everything already chosen.
k    = min(k, size(Z, 1));
pick = 1;
for j = 2:k
    dmin = inf(size(Z, 1), 1);
    for p = pick
        dmin = min(dmin, vecnorm(Z - Z(p, :), 2, 2));
    end
    dmin(pick) = -inf;
    [~, best]  = max(dmin);
    pick(end+1) = best; %#ok<AGROW>
end
end

% =========================================================================
function row = summarise(r, cfg, gi, rank, descriptors, z, med, dist)
row = table();
row.grade          = string(cfg.data.gradeLabels{gi});
row.gradeKey       = string(cfg.data.gradeOrder{gi});
row.rank           = rank;
row.isic_id        = r.isic_id;
row.typicalityDist = dist;
row.dice           = r.dice;
row.pixels_x       = r.pixels_x;
row.pixels_y       = r.pixels_y;
row.histoConfirmed = r.histoConfirmed;
for d = 1:numel(descriptors)
    row.(descriptors{d})               = r.(['otsu_' descriptors{d}]);
    row.([descriptors{d} '_gradeMed']) = med(d);
    row.([descriptors{d} '_z'])        = z(d);
end
row.license     = r.license;
row.attribution = r.attribution;
end

% =========================================================================
function render_example(r, cfg, gi, rank, descriptors, z, med, outDir)
imgPath  = fullfile(cfg.paths.images, char(r.isic_id + ".jpg"));
maskPath = fullfile(cfg.paths.masks,  char(r.isic_id + cfg.data.maskSuffix));
if ~isfile(imgPath) || ~isfile(maskPath)
    error('missing image or mask');
end

img        = imread(imgPath);
expertMask = seg.expert_route(maskPath);
otsuMask   = seg.otsu_route(img, cfg);

% Measure on the primary route, since that is what the exemplar illustrates.
out = morph.process_lesion(otsuMask, cfg);

stem = sprintf('%s_%02d_%s', cfg.data.gradeOrder{gi}, rank, r.isic_id);

% Self-contained copies so the folder can be used for figures without
% reaching back into the 10 GB raw archive.
copyfile(imgPath, fullfile(outDir, [stem '_image.jpg']));
imwrite(otsuMask, fullfile(outDir, [stem '_mask.png']));

scale = min(1, 900 / max(size(img, 1), size(img, 2)));
imgS  = imresize(img, scale);
expS  = imresize(expertMask, scale, 'nearest');
otsS  = imresize(otsuMask,   scale, 'nearest');

f = figure('Visible', 'off', 'Position', [0 0 1500 420]);
tiledlayout(1, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
imshow(imgS);
title(sprintf('%s  (%s)', r.isic_id, cfg.data.gradeLabels{gi}), ...
    'Interpreter', 'none');

nexttile;
imshow(imgS); hold on;
visboundaries(otsS, 'Color', 'r', 'LineWidth', 0.9);
visboundaries(expS, 'Color', 'g', 'LineWidth', 0.7, 'LineStyle', '--');
title(sprintf('red = automatic (primary), green = expert\nDice %.3f', r.dice));

nexttile;
b = out.boxcount;
if ~isempty(b) && isfield(b, 'epsAll') && ~isempty(b.epsAll)
    x = log(1 ./ b.epsAll);
    y = log(b.countsAll);
    plot(x, y, 'o', 'MarkerSize', 3); hold on;
    if ~isempty(b.fitIdx)
        xf = x(b.fitIdx); yf = y(b.fitIdx);
        plot(xf, yf, 'r.', 'MarkerSize', 12);
        plot(xf, polyval(polyfit(xf, yf, 1), xf), 'r-', 'LineWidth', 1.3);
    end
    grid on;
    xlabel('log(1/\epsilon)'); ylabel('log N(\epsilon)');
    title(sprintf('D_f = %.4f  (R^2 = %.4f)', out.Df, out.DfR2));
end

% The last panel is the evidence for the word "typical": every descriptor
% shown as a robust z against its own grade's median.
nexttile;
barh(z); hold on;
xline(0, 'k-');
xline(-1, 'k:'); xline(1, 'k:');
yticks(1:numel(descriptors));
yticklabels(descriptors);
set(gca, 'TickLabelInterpreter', 'none', 'YDir', 'reverse');
xlabel('robust z vs grade median');
xlim([-3 3]); grid on;
title('deviation from grade median');

sgtitle(sprintf(['%s exemplar %d  |  D_f %.4f (grade median %.4f)  ' ...
                 '|  compactness %.3f (median %.3f)'], ...
    cfg.data.gradeLabels{gi}, rank, ...
    r.otsu_Df, med(strcmp(descriptors, 'Df')), ...
    r.otsu_compactness, med(strcmp(descriptors, 'compactness'))), ...
    'Interpreter', 'none');

exportgraphics(f, fullfile(outDir, [stem '_panel.png']), 'Resolution', 130);
close(f);
end

% =========================================================================
function write_readme(outDir, cfg, k, descriptors)
fid = fopen(fullfile(outDir, 'README.md'), 'w');
c = onCleanup(@() fclose(fid));
fprintf(fid, '# Classical examples\n\n');
fprintf(fid, ['%d representative lesions per grade, selected automatically ' ...
    'from the Part 1 cohort.\n\n'], k);
fprintf(fid, '## Selection rule\n\n');
fprintf(fid, 'Candidates are accepted on both routes, carry an `ok` ');
fprintf(fid, 'box-count scaling flag, and\nare well segmented ');
fprintf(fid, '(Dice >= %.2f). Every descriptor is turned into a robust\n', ...
    cfg.qc.wellSegmentedDice);
fprintf(fid, 'z-score against *that grade''s* median and IQR:\n\n');
fprintf(fid, '    %s\n\n', strjoin(descriptors, ', '));
fprintf(fid, '**The reference median comes from the full accepted grade ');
fprintf(fid, 'population, not from the\nwell-segmented candidates.** Those ');
fprintf(fid, 'differ: the Dice gate preferentially keeps\nsmoother lesions, ');
fprintf(fid, 'because irregular borders are the ones the two routes disagree\n');
fprintf(fid, 'on. Scoring against the filtered pool would produce exemplars ');
fprintf(fid, 'typical of the\nfilter rather than of the published intervals ');
fprintf(fid, '(invasive compactness would land\nnear 1.35 against a reported ');
fprintf(fid, 'grade median of ~1.81). Any residual gap remains\nvisible in ');
fprintf(fid, 'each panel''s z-profile rather than being hidden by rescaling.\n\n');
fprintf(fid, 'Typicality is the **median** |z| across descriptors, so a lesion ');
fprintf(fid, 'central on most\ndescriptors is not disqualified by one unusual ');
fprintf(fid, 'value. The %d exemplars are then\nspread across the most-typical ', k);
fprintf(fid, 'pool by a greedy max-min rule, because the top\n%d by distance ', k);
fprintf(fid, 'are often near-duplicates.\n\n');
fprintf(fid, 'Measurements shown are from the **automatic (Otsu) route**, ');
fprintf(fid, 'which is the primary\nanalysis route for this study. The expert ');
fprintf(fid, 'mask is drawn dashed for comparison\nonly.\n\n');
fprintf(fid, '## Files\n\n');
fprintf(fid, '| suffix | contents |\n|---|---|\n');
fprintf(fid, '| `_panel.png` | image, dual-route boundaries, box-count fit, z-profile |\n');
fprintf(fid, '| `_image.jpg` | source dermoscopic image |\n');
fprintf(fid, '| `_mask.png` | automatic-route mask actually measured |\n');
fprintf(fid, '| `manifest.csv` | descriptor values, grade medians, z-scores, license |\n\n');
fprintf(fid, '## Licensing\n\n');
fprintf(fid, 'ISIC images are licensed **per image** (a CC-0 / CC-BY / CC-BY-NC ');
fprintf(fid, 'mix). The\n`license` and `attribution` columns of `manifest.csv` ');
fprintf(fid, 'apply to each file here and\nmust travel with any reuse. ');
fprintf(fid, 'CC-BY-NC images cannot be used in a commercial\ncontext.\n\n');
fprintf(fid, 'Regenerate with `make_classic_examples(%d)`.\n', k);
end
