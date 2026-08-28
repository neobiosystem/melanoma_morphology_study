function R = smoke_masks(nPerGrade)
%SMOKE_MASKS Pre-flight the measurement chain on real expert masks.
%
%   R = SMOKE_MASKS(nPerGrade) runs the expert route and the full descriptor
%   pipeline over a small stratified sample of real ISIC ground-truth masks
%   and prints the per-grade descriptor values, QC outcomes and timing.
%
%   Purpose is to catch problems that synthetic shapes cannot: real lesions
%   touch frame borders, come in wildly different resolutions, and produce
%   box-count curves whose scaling region may or may not exist. Run this
%   before committing to the full cohort.

if nargin < 1
    nPerGrade = 8;
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(fullfile(root, 'config'));
cfg = default_config();

T = readtable(cfg.paths.metadata, 'TextType', 'string');

% Map to the three grades.
grade = strings(height(T), 1);
for k = 1:size(cfg.data.gradeMap, 1)
    grade(T.diagnosis_3 == cfg.data.gradeMap{k, 1}) = cfg.data.gradeMap{k, 2};
end
T = T(grade ~= "", :);
grade = grade(grade ~= "");

rng(cfg.stats.seed);

rows = {};
for g = 1:numel(cfg.data.gradeOrder)
    thisGrade = cfg.data.gradeOrder{g};
    idx = find(grade == thisGrade);
    idx = idx(randperm(numel(idx), min(nPerGrade, numel(idx))));

    for k = idx(:)'
        maskPath = fullfile(cfg.paths.masks, ...
                            char(T.isic_id(k) + cfg.data.maskSuffix));
        if ~isfile(maskPath)
            continue
        end

        tic;
        [mask, mInfo] = seg.expert_route(maskPath);
        if ~mInfo.valid
            fprintf('%-14s %-9s SKIP %s\n', T.isic_id(k), thisGrade, mInfo.reason);
            continue
        end

        out = morph.process_lesion(mask, cfg);
        elapsed = toc;

        r.isic_id     = T.isic_id(k);
        r.grade       = string(thisGrade);
        r.accepted    = out.accepted;
        r.reason      = string(out.rejectReason);
        r.flags       = string(strjoin(out.flags, '|'));
        r.pixels      = sprintf('%dx%d', T.pixels_y(k), T.pixels_x(k));
        r.eqDiameter  = out.eqDiameter;
        r.Df          = out.Df;
        r.DfR2        = out.DfR2;
        r.DfFlag      = string(out.DfFlag);
        r.DfSpan      = out.DfSpanDecades;
        r.compactness = out.compactness;
        r.circularity = out.circularity;
        r.asymmetry   = out.asymmetry;
        r.solidity    = out.solidity;
        r.seconds     = elapsed;
        rows{end+1} = r; %#ok<AGROW>
    end
end

R = struct2table([rows{:}]);

if height(R) <= 40
    fprintf('\n===== per-lesion =====\n');
    disp(R(:, {'isic_id','grade','accepted','reason','pixels','eqDiameter', ...
               'Df','DfR2','DfFlag','compactness','circularity','asymmetry','seconds'}));
end

fprintf('\n===== QC summary =====\n');
fprintf('accepted %d / %d\n', nnz(R.accepted), height(R));
if any(~R.accepted)
    disp(groupcounts(R(~R.accepted, :), 'reason'));
end
fprintf('boxcount flags:\n');
disp(groupcounts(R, 'DfFlag'));
nonEmptyFlags = R.flags(R.flags ~= "");
if ~isempty(nonEmptyFlags)
    fprintf('lesion flags raised: %s\n', strjoin(cellstr(unique(nonEmptyFlags)), ', '));
end

fprintf('\n===== descriptor medians by grade (accepted only) =====\n');
A = R(R.accepted, :);
if ~isempty(A)
    disp(groupsummary(A, 'grade', 'median', ...
        {'Df','compactness','circularity','asymmetry','solidity'}));
end

% Does any descriptor actually separate the grades in this sample?
fprintf('\n===== separation check (Kruskal-Wallis across 3 grades) =====\n');
descriptors = {'Df','compactness','circularity','asymmetry','solidity'};
g = categorical(A.grade);
for d = descriptors
    v = A.(d{1});
    ok = isfinite(v);
    if nnz(ok) < 6 || numel(unique(g(ok))) < 2
        continue
    end
    p = kruskalwallis(v(ok), g(ok), 'off');
    q = quantile(v(ok), [0.25 0.5 0.75]);
    fprintf('  %-12s median=%7.4f IQR=[%7.4f %7.4f]  p=%.4g%s\n', ...
        d{1}, q(2), q(1), q(3), p, ternary(p < 0.05, '  *', ''));
end

fprintf('\nmean %.2f s/lesion -> projected %.1f min for 2223 lesions single-threaded\n', ...
    mean(R.seconds), mean(R.seconds) * 2223 / 60);
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
