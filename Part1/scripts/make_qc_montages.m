function make_qc_montages(n)
%MAKE_QC_MONTAGES Per-lesion visual QC sheets for manual spot-checking.
%
%   MAKE_QC_MONTAGES(n) renders n lesions as five-panel sheets:
%   original image, expert mask, Otsu mask, extracted boundary overlay, and
%   the box-counting log-log fit with its scaling region marked.
%
%   Numbers alone cannot tell you that a mask is following a hair rather
%   than the lesion, or that the Otsu route locked onto a vignette. These
%   sheets are how that gets caught. Lesions are sampled across grades and
%   deliberately biased toward flagged cases, since those are the ones worth
%   a human's attention.

if nargin < 1
    n = 30;
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root, 'config'));
cfg = default_config();

if ~isfolder(cfg.paths.qcFigures)
    mkdir(cfg.paths.qcFigures);
end

R = readtable(fullfile(cfg.paths.results, 'descriptors.csv'), 'TextType', 'string');
R.grade = categorical(R.grade, cfg.data.gradeOrder, 'Ordinal', true);

% readtable turns empty CSV cells into <missing>, which propagates through
% comparisons and concatenation instead of behaving like "". Normalise once
% here so downstream string handling is straightforward.
for v = ["expert_flags", "expert_DfFlag", "expert_rejectReason"]
    if ismember(v, R.Properties.VariableNames)
        R.(v)(ismissing(R.(v))) = "";
    end
end

% Prefer flagged or low-agreement lesions, then fill out with normal ones.
suspicious = R.expert_flags ~= "" | R.expert_DfFlag ~= "ok";
if ismember('dice', R.Properties.VariableNames)
    suspicious = suspicious | (isfinite(R.dice) & R.dice < cfg.qc.minRouteDice);
end

rng(cfg.stats.seed);
idxSusp = find(suspicious);
idxOk   = find(~suspicious);
idxSusp = idxSusp(randperm(numel(idxSusp), min(round(n*0.6), numel(idxSusp))));
idxOk   = idxOk(randperm(numel(idxOk), min(n - numel(idxSusp), numel(idxOk))));
sel     = [idxSusp; idxOk];

for k = sel(:)'
    try
        render_one(R(k, :), cfg);
    catch ME
        fprintf('  skipped %s: %s\n', R.isic_id(k), ME.message);
    end
end

fprintf('QC montages written to %s\n', cfg.paths.qcFigures);
end

% =========================================================================
function render_one(row, cfg)
imgPath  = fullfile(cfg.paths.images, char(row.isic_id + ".jpg"));
maskPath = fullfile(cfg.paths.masks,  char(row.isic_id + cfg.data.maskSuffix));

if ~isfile(imgPath) || ~isfile(maskPath)
    return
end

img        = imread(imgPath);
expertMask = seg.expert_route(maskPath);
otsuMask   = seg.otsu_route(img, cfg);
out        = morph.process_lesion(expertMask, cfg);

% Downsample the display copy; these sheets are for looking at, not measuring.
scale = min(1, 900 / max(size(img, 1), size(img, 2)));
imgS  = imresize(img, scale);
expS  = imresize(expertMask, scale, 'nearest');
otsS  = imresize(otsuMask,   scale, 'nearest');

f = figure('Visible', 'off', 'Position', [0 0 1700 420]);
tiledlayout(1, 5, 'Padding', 'compact', 'TileSpacing', 'compact');

gradeStr = string(row.grade);      % categorical will not concatenate with strings
idStr    = string(row.isic_id);

nexttile; imshow(imgS); title(sprintf('%s  (%s)', idStr, gradeStr));
nexttile; imshow(expS); title('expert mask');
nexttile; imshow(otsS); title(sprintf('Otsu mask (Dice %.3f)', row.dice));

nexttile;
imshow(imgS); hold on;
visboundaries(expS, 'Color', 'g', 'LineWidth', 0.8);
visboundaries(otsS, 'Color', 'r', 'LineWidth', 0.8);
title('green = expert, red = Otsu');

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
    xlabel('log(1/\epsilon)'); ylabel('log N');
    title(sprintf('D_f=%.4f R^2=%.4f\n%s', out.Df, out.DfR2, out.DfFlag));
end

flagText = row.expert_flags;
if ~ismissing(flagText) && flagText ~= ""
    sgtitle(sprintf('flags: %s', flagText), 'Interpreter', 'none');
end

exportgraphics(f, fullfile(cfg.paths.qcFigures, ...
    char("qc_" + gradeStr + "_" + idStr + ".png")), 'Resolution', 110);
close(f);
end
