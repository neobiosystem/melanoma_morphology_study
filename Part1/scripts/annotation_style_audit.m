function S = annotation_style_audit(nPerGrade)
%ANNOTATION_STYLE_AUDIT Quantify ISIC mask annotation style and its confound risk.
%
%   The ISIC 2018 ground-truth masks are not produced by one method. Some
%   are coarse polygon tracings whose boundary is mostly straight runs;
%   others are pixel-level flood fills that follow the true border. Because
%   the box-counting dimension responds to exactly that fine structure, a
%   polygon-traced lesion returns D_f ~ 1 no matter how irregular the real
%   lesion is.
%
%   That attenuation is tolerable if annotation style is independent of
%   diagnosis. It is FATAL if style correlates with grade, because then any
%   between-grade difference in D_f could be an artefact of who annotated
%   which lesion and with what tool. This script measures the style
%   distribution and tests it against grade.
%
%   Style proxy: the fraction of boundary vertices lying on a straight run
%   (turn angle ~ 0). A hand-drawn polygon is dominated by these; an
%   organically traced border is not.

if nargin < 1
    nPerGrade = 120;
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root, 'config'));
cfg = default_config();

T = readtable(cfg.paths.metadata, 'TextType', 'string');
grade = strings(height(T), 1);
for k = 1:size(cfg.data.gradeMap, 1)
    grade(T.diagnosis_3 == cfg.data.gradeMap{k, 1}) = cfg.data.gradeMap{k, 2};
end
T = T(grade ~= "", :);  grade = grade(grade ~= "");

rng(cfg.stats.seed);
rows = {};
for g = cfg.data.gradeOrder
    idx = find(grade == g{1});
    idx = idx(randperm(numel(idx), min(nPerGrade, numel(idx))));
    for k = idx(:)'
        p = fullfile(cfg.paths.masks, char(T.isic_id(k) + cfg.data.maskSuffix));
        if ~isfile(p), continue; end
        mask = imread(p) > 127;
        mask = bwareafilt(logical(mask), 1);
        if ~any(mask(:)), continue; end
        mask = imfill(mask, 'holes');

        dEq = morph.equivalent_diameter(mask);
        if dEq < cfg.qc.minEqDiameterPx, continue; end
        C   = morph.boundary_trace(mask);
        Cs  = morph.smooth_contour(C, round(cfg.smooth.primaryRatio * dEq));

        r.isic_id   = T.isic_id(k);
        r.grade     = string(g{1});
        r.dEq       = dEq;
        r.collinear = collinear_fraction(C);
        r.roughness = morph.contour_perimeter(C) / morph.contour_perimeter(Cs);
        rows{end+1} = r; %#ok<AGROW>
    end
end

S = struct2table([rows{:}]);
S.isPolygon = S.collinear > 0.85;   % clearly a traced polygon

fprintf('\n===== annotation style, n = %d =====\n', height(S));
fprintf('collinear-vertex fraction: median %.3f, IQR [%.3f %.3f]\n', ...
    median(S.collinear), quantile(S.collinear,0.25), quantile(S.collinear,0.75));
fprintf('polygon-traced masks (collinear > 0.85): %d / %d = %.1f%%\n\n', ...
    nnz(S.isPolygon), height(S), 100*mean(S.isPolygon));

fprintf('--- polygon-traced share BY GRADE (the confound test) ---\n');
G = groupsummary(S, 'grade', {'mean','median'}, {'isPolygon','collinear','roughness'});
disp(G);

% Chi-square on style vs grade.
tbl = crosstab(categorical(S.grade), S.isPolygon);
[~, chi2, pChi] = crosstab(categorical(S.grade), S.isPolygon);
fprintf('contingency (rows = %s):\n', strjoin(cfg.data.gradeOrder, ', '));
disp(tbl);
fprintf('chi-square = %.3f, p = %.4f\n', chi2, pChi);
if pChi < 0.05
    fprintf(['*** STYLE IS CONFOUNDED WITH GRADE. D_f differences from ' ...
             'expert masks would be uninterpretable. ***\n']);
else
    fprintf(['Style appears independent of grade: D_f is attenuated ' ...
             'toward 1 but not differentially biased.\n']);
end

fprintf('\n--- roughness (raw/smoothed perimeter) by style ---\n');
fprintf('polygon-traced : median %.4f\n', median(S.roughness(S.isPolygon)));
fprintf('organic        : median %.4f\n', median(S.roughness(~S.isPolygon)));
p = ranksum(S.roughness(S.isPolygon), S.roughness(~S.isPolygon));
fprintf('rank-sum p = %.3g\n', p);
end

function frac = collinear_fraction(C)
if size(C,1) < 5, frac = NaN; return; end
v1 = C(2:end-1,:) - C(1:end-2,:);
v2 = C(3:end,:)   - C(2:end-1,:);
cosang = sum(v1.*v2,2) ./ max(vecnorm(v1,2,2).*vecnorm(v2,2,2), eps);
frac = mean(cosang > 0.999);
end
