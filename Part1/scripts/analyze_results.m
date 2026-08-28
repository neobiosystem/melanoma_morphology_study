function A = analyze_results(route)
%ANALYZE_RESULTS Per-grade intervals, group tests and discriminative power.
%
%   A = ANALYZE_RESULTS() analyses the expert route (the primary analysis).
%   A = ANALYZE_RESULTS('otsu') analyses the automatic route instead.
%
%   Reads results/descriptors.csv, produces the per-grade interval tables
%   that are Part 1's headline deliverable, and writes
%   results/calibration_targets.json for Part 3 to load directly.
%
%   Omnibus p values across descriptors are FDR-corrected together, and the
%   pairwise post-hocs within each descriptor are Holm-corrected, so the
%   two levels of multiplicity are handled separately rather than ignored.

if nargin < 1
    route = 'expert';
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root, 'config'));
cfg = default_config();

R = readtable(fullfile(cfg.paths.results, 'descriptors.csv'), 'TextType', 'string');
R.grade = categorical(R.grade, cfg.data.gradeOrder, 'Ordinal', true);

acceptedCol = [route '_accepted'];
R = R(logical(R.(acceptedCol)), :);
fprintf('Analysing %d accepted lesions (%s route).\n\n', height(R), route);

descriptors = cfg.stats.descriptors;
A = struct('route', route, 'intervals', struct(), 'tests', struct(), ...
           'power', struct(), 'n', height(R));

omnibusP = nan(numel(descriptors), 1);

for d = 1:numel(descriptors)
    name = descriptors{d};
    col  = [route '_' name];
    if ~ismember(col, R.Properties.VariableNames)
        continue
    end

    v = R.(col);
    g = R.grade;

    A.intervals.(name) = stats.grade_intervals(v, g, cfg);
    A.tests.(name)     = stats.group_tests(v, g, cfg);
    A.power.(name)     = stats.discriminative_power(v, g, cfg);
    omnibusP(d)        = A.tests.(name).p;
end

% FDR across descriptors.
omnibusAdj = stats.adjust_p(omnibusP, cfg.stats.omnibusCorrect);

% ---- Report ------------------------------------------------------------
fprintf('%-22s %8s %8s %9s %8s %8s %10s\n', ...
    'descriptor','KW p','p(FDR)','eta2_H','AUC','balAUC','C-index');
fprintf('%s\n', repmat('-', 1, 82));
for d = 1:numel(descriptors)
    name = descriptors{d};
    if ~isfield(A.tests, name), continue; end
    t = A.tests.(name);
    p = A.power.(name);
    fprintf('%-22s %8.4g %8.4g %9.4f %8.3f %8.3f %10.3f\n', ...
        name, t.p, omnibusAdj(d), t.etaSqH, p.auc, p.balancedAuc, p.cIndex);
end

fprintf('\n===== per-grade intervals =====\n');
for d = 1:numel(descriptors)
    name = descriptors{d};
    if ~isfield(A.intervals, name), continue; end
    fprintf('\n--- %s ---\n', name);
    disp(A.intervals.(name));
end

fprintf('\n===== pairwise post-hoc (Dunn, Holm-corrected) =====\n');
for d = 1:numel(descriptors)
    name = descriptors{d};
    if ~isfield(A.tests, name) || isempty(A.tests.(name).postHoc), continue; end
    ph = A.tests.(name).postHoc;
    sig = ph(ph.pAdjusted < cfg.stats.alpha, :);
    if ~isempty(sig)
        fprintf('\n--- %s ---\n', name);
        disp(sig(:, {'groupA','groupB','medianA','medianB','z','pAdjusted','rankBiserial'}));
    end
end

A.omnibusAdjusted = omnibusAdj;
save(fullfile(cfg.paths.results, ['analysis_' route '.mat']), 'A');

export_calibration_targets(A, cfg, route);
end

% =========================================================================
function export_calibration_targets(A, cfg, route)
%EXPORT_CALIBRATION_TARGETS Write the Part 3 target file.
%
%   Machine-readable so Part 3 loads the intervals directly instead of
%   transcribing numbers, which is where quiet errors creep into a
%   multi-stage study.

targets = struct();
if strcmp(route, 'expert')
    maskSource = 'expert ground-truth masks (polygon-traced)';
else
    maskSource = 'automatic Otsu segmentation of image pixels';
end
targets.source        = ['Part1 ISIC 2018 Task1-2 Training, ' maskSource];
targets.route         = route;
targets.numLesions    = A.n;
targets.gradeOrder    = cfg.data.gradeOrder;
targets.descriptors   = struct();

names = fieldnames(A.intervals);
for k = 1:numel(names)
    T = A.intervals.(names{k});
    entry = struct();
    for r = 1:height(T)
        gname = matlab.lang.makeValidName(char(T.grade(r)));
        entry.(gname) = struct( ...
            'n',      T.n(r), ...
            'median', T.median(r), ...
            'iqr',    [T.q1(r), T.q3(r)], ...
            'refInterval', [T.refLo(r), T.refHi(r)], ...
            'mean',   T.mean(r), ...
            'sd',     T.sd(r));
    end
    if isfield(A.power, names{k})
        entry.auc     = A.power.(names{k}).auc;
        entry.cIndex  = A.power.(names{k}).cIndex;
    end
    targets.descriptors.(names{k}) = entry;
end

% Route-specific filename: both routes are exported, and a shared name would
% mean whichever ran last silently became "the" calibration target.
outPath = fullfile(cfg.paths.results, ['calibration_targets_' route '.json']);
fid = fopen(outPath, 'w');
fprintf(fid, '%s', jsonencode(targets, 'PrettyPrint', true));
fclose(fid);
fprintf('\nWrote %s\n', outPath);
end
