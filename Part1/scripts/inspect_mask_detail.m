function inspect_mask_detail()
%INSPECT_MASK_DETAIL Test whether expert masks preserve fine border structure.
%
%   D_f came out ~1.02 for all three grades, far below the 1.1-1.5 range the
%   melanoma border literature reports. Either the measurement is wrong or
%   the ISIC ground-truth masks are smooth annotations that do not carry
%   pixel-level border detail. This script distinguishes the two.
%
%   The discriminating quantity is how much perimeter survives smoothing. A
%   genuinely rough border sheds perimeter quickly; an already-smooth
%   polygon barely changes. It also writes a zoomed boundary crop so the
%   geometry can be seen directly.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root, 'config')); addpath(fullfile(root, 'tests'));
cfg = default_config();

T = readtable(cfg.paths.metadata, 'TextType', 'string');
grade = strings(height(T), 1);
for k = 1:size(cfg.data.gradeMap, 1)
    grade(T.diagnosis_3 == cfg.data.gradeMap{k, 1}) = cfg.data.gradeMap{k, 2};
end
T = T(grade ~= "", :);  grade = grade(grade ~= "");

rng(cfg.stats.seed);
sel = [];
for g = cfg.data.gradeOrder
    idx = find(grade == g{1});
    sel = [sel; idx(randperm(numel(idx), 12))]; %#ok<AGROW>
end

fprintf('%-14s %-9s %8s %8s %9s %9s %9s\n', ...
    'isic_id','grade','dEq','rawP/smP','decaySlp','vertices','collinear');

decay = []; ratio = []; collin = [];
for k = sel(:)'
    p = fullfile(cfg.paths.masks, char(T.isic_id(k) + cfg.data.maskSuffix));
    if ~isfile(p), continue; end
    mask = imread(p) > 127;
    mask = imfill(bwareafilt(logical(mask), 1), 'holes');
    if ~any(mask(:)), continue; end

    dEq = morph.equivalent_diameter(mask);
    C   = morph.boundary_trace(mask);
    Cs  = morph.smooth_contour(C, round(cfg.smooth.primaryRatio * dEq));
    rP  = morph.contour_perimeter(C);
    sP  = morph.contour_perimeter(Cs);
    ms  = morph.multi_scale(mask, dEq, cfg);

    % Collinearity: share of contour vertices whose turn angle is ~0, i.e.
    % lying on a straight run. A hand-drawn polygon is dominated by these.
    frac = collinear_fraction(C);

    fprintf('%-14s %-9s %8.1f %8.4f %9.4f %9d %9.3f\n', ...
        T.isic_id(k), grade(k), dEq, rP/sP, ms.perimeterDecaySlope, ...
        size(C,1), frac);

    decay(end+1)  = ms.perimeterDecaySlope; %#ok<AGROW>
    ratio(end+1)  = rP/sP;                  %#ok<AGROW>
    collin(end+1) = frac;                   %#ok<AGROW>
end

fprintf('\nmedian raw/smoothed perimeter ratio : %.4f\n', median(ratio));
fprintf('median perimeter-decay slope        : %.4f\n', median(decay));
fprintf('median collinear vertex fraction    : %.3f\n', median(collin));

fprintf(['\nReference: a rasterised CIRCLE (perfectly smooth) and a KOCH ' ...
         'curve (maximally rough)\n']);
d = shapes.disk(500);
Cd = morph.boundary_trace(d);
dq = morph.equivalent_diameter(d);
msd = morph.multi_scale(d, dq, cfg);
fprintf('  circle : raw/sm=%.4f  decay=%.4f  collinear=%.3f\n', ...
    morph.contour_perimeter(Cd)/morph.contour_perimeter( ...
        morph.smooth_contour(Cd, round(cfg.smooth.primaryRatio*dq))), ...
    msd.perimeterDecaySlope, collinear_fraction(Cd));

% Visual: zoom on a boundary segment so the geometry is directly visible.
p = fullfile(cfg.paths.masks, char(T.isic_id(sel(1)) + cfg.data.maskSuffix));
mask = imread(p) > 127;
C = morph.boundary_trace(bwareafilt(logical(mask),1));
mid = round(size(C,1)/2);
seg = C(max(1,mid-150):min(size(C,1),mid+150), :);

f = figure('Visible','off','Position',[0 0 1400 500]);
subplot(1,2,1); imshow(mask); title('full expert mask');
subplot(1,2,2);
plot(seg(:,2), seg(:,1), 'k.-', 'MarkerSize', 8); axis equal; grid on;
title('300 consecutive boundary pixels (zoom)');
outPng = fullfile(cfg.paths.qcFigures, 'mask_boundary_zoom.png');
exportgraphics(f, outPng, 'Resolution', 110);
close(f);
fprintf('\nWrote %s\n', outPng);
end

function frac = collinear_fraction(C)
if size(C,1) < 5, frac = NaN; return; end
v1 = C(2:end-1,:) - C(1:end-2,:);
v2 = C(3:end,:)   - C(2:end-1,:);
cosang = sum(v1.*v2,2) ./ max(vecnorm(v1,2,2).*vecnorm(v2,2,2), eps);
frac = mean(cosang > 0.999);
end
