function [R, cfg] = run_analysis(varargin)
%RUN_ANALYSIS Measure every cohort lesion by both segmentation routes.
%
%   R = RUN_ANALYSIS() processes the full manifest and writes
%   results/descriptors.csv and results/qc_report.csv.
%
%   Name-value options:
%       'Limit'      - process only the first N lesions (0 = all, default)
%       'Routes'     - {'expert'}, {'otsu'} or both (default both)
%       'UseParfor'  - override cfg.run.useParfor
%       'Stratified' - with Limit, sample evenly across grades (default true)
%
%   Both routes run through the identical MORPH.PROCESS_LESION, differing
%   only in where the mask comes from. Any measured disagreement is
%   therefore attributable to segmentation and not to the measurement code.
%
%   Each lesion is wrapped in try/catch so a single unreadable file cannot
%   abort a run of several thousand.

p = inputParser;
p.addParameter('Limit', 0, @isnumeric);
p.addParameter('Routes', {'expert', 'otsu'}, @iscell);
p.addParameter('UseParfor', [], @(x) isempty(x) || islogical(x));
p.addParameter('Stratified', true, @islogical);
p.parse(varargin{:});
opt = p.Results;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(fullfile(root, 'config'));
cfg = default_config();

if ~isempty(opt.UseParfor)
    cfg.run.useParfor = opt.UseParfor;
end

doExpert = any(strcmp(opt.Routes, 'expert'));
doOtsu   = any(strcmp(opt.Routes, 'otsu'));

% ---- Cohort ------------------------------------------------------------
fprintf('Building manifest...\n');
M = io.build_manifest(cfg, true);

if opt.Limit > 0 && opt.Limit < height(M)
    rng(cfg.stats.seed);
    if opt.Stratified
        idx = [];
        perGrade = ceil(opt.Limit / numel(cfg.data.gradeOrder));
        for g = cfg.data.gradeOrder
            gi = find(M.grade == g{1});
            gi = gi(randperm(numel(gi), min(perGrade, numel(gi))));
            idx = [idx; gi]; %#ok<AGROW>
        end
        M = M(sort(idx), :);
    else
        M = M(1:opt.Limit, :);
    end
    fprintf('Limited to %d lesions.\n', height(M));
end

n = height(M);
fprintf('Processing %d lesions (expert=%d, otsu=%d)...\n', n, doExpert, doOtsu);

imagePaths = cellstr(M.imagePath);
maskPaths  = cellstr(M.maskPath);
results    = cell(n, 1);

tStart = tic;

if cfg.run.useParfor
    pool = gcp('nocreate');
    if isempty(pool)
        % Clamp to what the cluster profile actually permits. Requesting more
        % than NumWorkers is a hard error, not a silent downgrade, so a
        % hardcoded count makes the run machine-specific for no benefit.
        clust   = parcluster('local');
        nWanted = min(cfg.run.numWorkers, clust.NumWorkers);
        pool    = parpool(clust, nWanted);
    end
    fprintf('Using %d workers.\n', pool.NumWorkers);

    parfor k = 1:n
        % Workers do not inherit the client path, so the packages must be
        % made visible inside the loop.
        addpath(root); addpath(fullfile(root, 'config'));
        results{k} = process_one(imagePaths{k}, maskPaths{k}, cfg, doExpert, doOtsu);
    end
else
    for k = 1:n
        results{k} = process_one(imagePaths{k}, maskPaths{k}, cfg, doExpert, doOtsu);
        if mod(k, 100) == 0
            fprintf('  %d / %d  (%.1f min elapsed)\n', k, n, toc(tStart)/60);
        end
    end
end

elapsed = toc(tStart);
fprintf('Done in %.1f min (%.2f s/lesion).\n', elapsed/60, elapsed/n);

% ---- Assemble ----------------------------------------------------------
R = struct2table([results{:}]);
R = [M(:, {'isic_id','grade','diagnosis_3','histoConfirmed', ...
           'pixels_x','pixels_y','license'}), R];

if ~isfolder(cfg.paths.results)
    mkdir(cfg.paths.results);
end

writetable(R, fullfile(cfg.paths.results, 'descriptors.csv'));
save(fullfile(cfg.paths.results, 'run_config.mat'), 'cfg', 'elapsed');

% ---- QC report ---------------------------------------------------------
qcVars = intersect({'isic_id','grade','expert_accepted','expert_rejectReason', ...
    'expert_flags','expert_DfFlag','expert_DfR2','expert_borderContactFraction', ...
    'otsu_accepted','otsu_rejectReason','otsu_flags','dice','jaccard','areaRatio'}, ...
    R.Properties.VariableNames, 'stable');
writetable(R(:, qcVars), fullfile(cfg.paths.results, 'qc_report.csv'));

print_summary(R, doExpert, doOtsu, cfg);
end

% =========================================================================
function row = process_one(imagePath, maskPath, cfg, doExpert, doOtsu)
%PROCESS_ONE Measure a single lesion by both routes, never throwing.

row = blank_row();

% ---- Expert route ------------------------------------------------------
expertMask = [];
if doExpert
    try
        [expertMask, info] = seg.expert_route(maskPath);
        if info.valid
            out = morph.process_lesion(expertMask, cfg);
            row = merge_fields(row, flatten(out, 'expert_'));
        else
            row.expert_rejectReason = string(info.reason);
        end
    catch ME
        row.expert_rejectReason = "error:" + string(ME.identifier);
    end
end

% ---- Otsu route --------------------------------------------------------
otsuMask = [];
if doOtsu
    try
        img = imread(imagePath);
        [otsuMask, info] = seg.otsu_route(img, cfg);
        if info.valid
            out = morph.process_lesion(otsuMask, cfg);
            row = merge_fields(row, flatten(out, 'otsu_'));
            row.otsu_hairFraction = info.hairFraction;
            row.otsu_usedActiveContour = info.usedActiveContour;
        else
            row.otsu_rejectReason = string(info.reason);
        end
    catch ME
        row.otsu_rejectReason = "error:" + string(ME.identifier);
    end
end

% ---- Route overlap -----------------------------------------------------
if doExpert && doOtsu && ~isempty(expertMask) && ~isempty(otsuMask)
    try
        ov = agree.overlap_metrics(expertMask, otsuMask);
        row.dice          = ov.dice;
        row.jaccard       = ov.jaccard;
        row.areaRatio     = ov.areaRatio;
        row.hausdorffNorm = ov.hausdorffNorm;
    catch
        % Leave overlap as NaN.
    end
end
end

% =========================================================================
function f = flatten(out, prefix)
%FLATTEN Select the scalar descriptor fields into a prefixed struct.

keep = {'accepted','area','eqDiameter','eqDiameterNative','frameFraction', ...
        'borderContactFraction','downsampled','numBoundaryPixels', ...
        'Df','DfStdErr','DfR2','DfSpanDecades','DfNumPoints','DfCrossCheck', ...
        'perimeterRaw','perimeterSmoothed','compactness','compactnessRaw', ...
        'circularity','radialCv','asymmetry','asymmetryMajor','asymmetryMinor', ...
        'asymmetryCentral','elongation','solidity','convexity','perimeterDecaySlope'};

f = struct();
for c = keep
    f.([prefix c{1}]) = out.(c{1});
end
f.([prefix 'rejectReason']) = string(out.rejectReason);
f.([prefix 'DfFlag'])       = string(out.DfFlag);
f.([prefix 'flags'])        = string(strjoin(out.flags, '|'));
end

% =========================================================================
function row = blank_row()
row = struct();
for prefix = {'expert_', 'otsu_'}
    pre = prefix{1};
    row.([pre 'accepted']) = false;
    for c = {'area','eqDiameter','eqDiameterNative','frameFraction', ...
             'borderContactFraction','numBoundaryPixels', ...
             'Df','DfStdErr','DfR2','DfSpanDecades','DfNumPoints','DfCrossCheck', ...
             'perimeterRaw','perimeterSmoothed','compactness','compactnessRaw', ...
             'circularity','radialCv','asymmetry','asymmetryMajor','asymmetryMinor', ...
             'asymmetryCentral','elongation','solidity','convexity','perimeterDecaySlope'}
        row.([pre c{1}]) = NaN;
    end
    row.([pre 'downsampled'])  = false;
    row.([pre 'rejectReason']) = "";
    row.([pre 'DfFlag'])       = "";
    row.([pre 'flags'])        = "";
end
row.otsu_hairFraction      = NaN;
row.otsu_usedActiveContour = false;
row.dice          = NaN;
row.jaccard       = NaN;
row.areaRatio     = NaN;
row.hausdorffNorm = NaN;
end

% =========================================================================
function a = merge_fields(a, b)
f = fieldnames(b);
for k = 1:numel(f)
    a.(f{k}) = b.(f{k});
end
end

% =========================================================================
function print_summary(R, doExpert, doOtsu, cfg)
fprintf('\n===== cohort =====\n');
summary(R.grade);

if doExpert
    fprintf('\nexpert route accepted: %d / %d\n', nnz(R.expert_accepted), height(R));
    rr = R.expert_rejectReason(R.expert_rejectReason ~= "");
    if ~isempty(rr)
        disp(groupcounts(table(rr), 'rr'));
    end
end
if doOtsu
    fprintf('otsu route accepted:   %d / %d\n', nnz(R.otsu_accepted), height(R));
end
if doExpert && doOtsu
    d = R.dice(isfinite(R.dice));
    if ~isempty(d)
        fprintf('\nDice expert-vs-otsu: median %.3f, IQR [%.3f %.3f]\n', ...
            median(d), quantile(d,0.25), quantile(d,0.75));
    end
end

% Both routes are summarised, not just the primary one. In a dual-route
% design the whole point is to see whether the two tell the same story, and
% printing only the expert route hides the comparison that justifies running
% the second route at all.
if doExpert
    fprintf('\n===== descriptor medians by grade (expert route, accepted) =====\n');
    A = R(R.expert_accepted, :);
    if ~isempty(A)
        disp(groupsummary(A, 'grade', 'median', ...
            {'expert_Df','expert_compactness','expert_circularity', ...
             'expert_asymmetry','expert_solidity'}));
    end
end

if doOtsu
    fprintf('\n===== descriptor medians by grade (otsu route, accepted) =====\n');
    B = R(R.otsu_accepted, :);
    if ~isempty(B)
        disp(groupsummary(B, 'grade', 'median', ...
            {'otsu_Df','otsu_compactness','otsu_circularity', ...
             'otsu_asymmetry','otsu_solidity'}));
    end
end
end
