function T = tune_otsu(n)
%TUNE_OTSU Diagnose the automatic route's systematic under-segmentation.
%
%   The automatic masks come out around 55-70% of expert mask area, and the
%   shortfall grows with severity. Under-segmentation that tracks grade is
%   dangerous: it would deform exactly the descriptors under study, in a
%   grade-dependent way.
%
%   Rather than guess, this sweeps the specific stages that can shrink a
%   mask and reports Dice and area ratio for each variant, so the cause is
%   identified before anything is re-run at cohort scale.
%
%   Variants:
%       baseline      - current configuration
%       noActive      - active-contour refinement off (it can shrink)
%       noOpen        - morphological opening off (erodes the periphery)
%       smallOpen     - opening radius halved
%       otsuOnCrop    - threshold computed inside the aperture only, so dark
%                       vignette corners cannot drag the threshold down
%       noActiveNoOpen- both suspects removed together

if nargin < 1
    n = 40;
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root, 'config'));
cfg = default_config();

M = io.build_manifest(cfg, true);
rng(cfg.stats.seed);
idx = [];
for g = cfg.data.gradeOrder
    gi  = find(M.grade == g{1});
    gi  = gi(randperm(numel(gi), min(ceil(n/3), numel(gi))));
    idx = [idx; gi]; %#ok<AGROW>
end
M = M(sort(idx), :);
fprintf('Tuning on %d lesions.\n\n', height(M));

variants = build_variants(cfg);
names    = fieldnames(variants);

acc = struct();
for v = 1:numel(names)
    acc.(names{v}) = struct('dice', [], 'areaRatio', [], 'grade', {{}}, 'fail', 0);
end

for k = 1:height(M)
    img = imread(char(M.imagePath(k)));
    exp = seg.expert_route(char(M.maskPath(k)));
    if ~any(exp(:)), continue; end

    for v = 1:numel(names)
        try
            m = seg.otsu_route(img, variants.(names{v}));
            if ~any(m(:))
                acc.(names{v}).fail = acc.(names{v}).fail + 1;
                continue
            end
            ov = agree.overlap_metrics(exp, m);
            acc.(names{v}).dice(end+1)      = ov.dice;
            acc.(names{v}).areaRatio(end+1) = ov.areaRatio;
            acc.(names{v}).grade{end+1}     = char(M.grade(k));
        catch
            acc.(names{v}).fail = acc.(names{v}).fail + 1;
        end
    end

    if mod(k, 10) == 0
        fprintf('  %d / %d\n', k, height(M));
    end
end

% ---- Report ------------------------------------------------------------
fprintf('\n%-16s %8s %10s %10s %12s %8s\n', ...
    'variant', 'n', 'medDice', 'meanDice', 'medAreaRat', 'fails');
fprintf('%s\n', repmat('-', 1, 70));

rows = {};
for v = 1:numel(names)
    a = acc.(names{v});
    if isempty(a.dice), continue; end
    fprintf('%-16s %8d %10.3f %10.3f %12.3f %8d\n', names{v}, numel(a.dice), ...
        median(a.dice), mean(a.dice), median(a.areaRatio), a.fail);
    rows{end+1} = struct('variant', string(names{v}), 'n', numel(a.dice), ...
        'medDice', median(a.dice), 'meanDice', mean(a.dice), ...
        'medAreaRatio', median(a.areaRatio), 'fails', a.fail); %#ok<AGROW>
end
T = struct2table([rows{:}]);

% ---- Grade-stratified bias, for EVERY variant --------------------------
% This, not Dice, is the selection criterion. A route that under-segments
% uniformly shifts all grades together and largely cancels in a between-grade
% comparison; one that under-segments differentially deforms the comparison
% itself. So what matters is the SPREAD of area ratio across grades, and a
% variant with worse average Dice but flatter bias is the better instrument.
fprintf('\n%-20s %24s %24s %10s\n', '', 'median Dice by grade', ...
    'median areaRatio by grade', 'bias');
fprintf('%-20s %7s %7s %7s %8s %7s %7s %10s\n', 'variant', ...
    'ben', 'situ', 'inv', 'ben', 'situ', 'inv', 'spread');
fprintf('%s\n', repmat('-', 1, 84));

for v = 1:numel(names)
    a = acc.(names{v});
    if isempty(a.dice), continue; end
    dm = nan(1, 3); am = nan(1, 3);
    for k = 1:numel(cfg.data.gradeOrder)
        sel = strcmp(a.grade, cfg.data.gradeOrder{k});
        if ~any(sel), continue; end
        dm(k) = median(a.dice(sel));
        am(k) = median(a.areaRatio(sel));
    end
    fprintf('%-20s %7.3f %7.3f %7.3f %8.3f %7.3f %7.3f %10.3f\n', ...
        names{v}, dm, am, max(am) - min(am));
end

fprintf(['\nSelection rule: prefer the smallest bias spread among variants\n' ...
         'with acceptable Dice. Uniform under-segmentation largely cancels in a\n' ...
         'between-grade comparison; grade-dependent under-segmentation does not.\n']);

writetable(T, fullfile(cfg.paths.results, 'otsu_tuning.csv'));
end

% =========================================================================
function variants = build_variants(cfg)
%BUILD_VARIANTS Round 2: combinations of the two fixes that round 1 isolated.
%
%   Round 1 found two independent improvements that were never tested
%   together. Aperture-restricted thresholding gave the largest gain in MEAN
%   Dice (0.538 -> 0.623) against little change in the median, which is the
%   signature of rescuing catastrophic failures rather than nudging typical
%   cases -- consistent with dark vignette corners dragging the threshold.
%   Dropping the morphological opening gave the best area ratio (0.666 ->
%   0.753), consistent with the opening eroding the lesion periphery.
%
%   These address different failure modes, so the combination is the
%   candidate worth testing rather than either alone.
variants = struct();

variants.baseline = cfg;

c = cfg; c.seg.thresholdInsideAperture = true;
variants.crop = c;

c = cfg; c.seg.thresholdInsideAperture = true; c.seg.openRadiusRatio = 0;
variants.cropNoOpen = c;

c = cfg; c.seg.thresholdInsideAperture = true; c.seg.useActiveContour = false;
variants.cropNoActive = c;

c = cfg; c.seg.thresholdInsideAperture = true;
c.seg.openRadiusRatio = 0; c.seg.useActiveContour = false;
variants.cropNoOpenNoActive = c;

% Opening exists to strip hairs; with it disabled, hair suppression has to
% carry that load alone. Tested explicitly so the tradeoff is visible rather
% than assumed away.
c = cfg; c.seg.thresholdInsideAperture = true;
c.seg.openRadiusRatio = cfg.seg.openRadiusRatio / 3;
variants.cropTinyOpen = c;
end
