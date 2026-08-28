function out = multi_scale(mask, dEq, cfg)
%MULTI_SCALE Convex-deficiency and perimeter-decay irregularity measures.
%
%   out = MULTI_SCALE(mask, dEq, cfg) returns a struct with fields:
%       solidity            - area / convex-hull area, 1 for a convex shape
%       convexity           - hull perimeter / contour perimeter, 1 if convex
%       sigmaPx             - smoothing scales used, in pixels
%       perimeters          - perimeter at each smoothing scale
%       perimeterRatio      - perimeters normalised by the unsmoothed value
%       perimeterDecaySlope - slope of log P against log(sigma/dEq)
%
%   The decay slope is the multi-scale irregularity descriptor: a boundary
%   carrying fine detail sheds perimeter quickly as the smoothing window
%   grows, giving a steeper negative slope, whereas a smooth boundary barely
%   changes. It is also the evidence that the perimeter-bias correction in
%   MORPH.SMOOTH_CONTOUR is doing real work -- if a between-grade difference
%   in compactness survives across smoothing scales it is geometric, and if
%   it collapses as soon as sigma exceeds a few pixels it was a rasterisation
%   artefact all along.
%
%   Scales are fractions of dEq rather than fixed pixel counts, so the curve
%   is comparable across the cohort's wide resolution spread.
%
%   See also MORPH.SMOOTH_CONTOUR, MORPH.COMPACTNESS.

out = struct('solidity', NaN, 'convexity', NaN, 'sigmaPx', [], ...
             'perimeters', [], 'perimeterRatio', [], ...
             'perimeterDecaySlope', NaN);

C = morph.boundary_trace(mask);
if size(C, 1) < 4 || ~isfinite(dEq) || dEq <= 0
    return
end

% ---- Convex deficiency -------------------------------------------------
% Numerator and denominator must come from the SAME representation. Pairing
% a pixel count (nnz) with the polygon area of the traced contour mixes two
% different area measures: the pixel count includes the full area of every
% boundary pixel while the contour runs through their centres, so the ratio
% creeps above 1 for convex shapes -- which is impossible by definition and
% was showing up in the reference intervals.
try
    hullIdx  = convhull(C(:, 1), C(:, 2));
    hullC    = C(hullIdx, :);
    hullArea = morph.contour_area(hullC);
catch
    return
end

contourArea = morph.contour_area(C);
if hullArea > 0
    out.solidity = contourArea / hullArea;
end

% ---- Perimeter decay across smoothing scales ---------------------------
ratios  = cfg.smooth.sigmaRatios(:)';
sigmaPx = ratios * dEq;

perims = nan(size(sigmaPx));
for k = 1:numel(sigmaPx)
    Cs         = morph.smooth_contour(C, round(sigmaPx(k)));
    perims(k)  = morph.contour_perimeter(Cs);
end

out.sigmaPx    = sigmaPx;
out.perimeters = perims;

% Convexity uses the perimeter at the primary smoothing scale so that both
% terms of the ratio describe comparably smooth geometry. That means the hull
% must be taken from the SMOOTHED contour as well: pairing the raw contour's
% hull perimeter with a smoothed contour perimeter compares a rough outline
% against a smooth one, and since smoothing shortens the contour the ratio
% can exceed 1 -- impossible, because the convex hull is the shortest closed
% curve enclosing a set.
primaryWin = round(cfg.smooth.primaryRatio * dEq);
Cprimary   = morph.smooth_contour(C, primaryWin);
Pprimary   = morph.contour_perimeter(Cprimary);
if Pprimary > 0
    try
        hIdx = convhull(Cprimary(:, 1), Cprimary(:, 2));
        out.convexity = morph.contour_perimeter(Cprimary(hIdx, :)) / Pprimary;
    catch
        % Degenerate (collinear) smoothed contour; leave convexity NaN.
    end
end

if perims(1) > 0
    out.perimeterRatio = perims / perims(1);
end

% Fit the decay over the interior scales; sigma = 0 has no logarithm.
useIdx = find(ratios > 0 & isfinite(perims) & perims > 0);
if numel(useIdx) >= 3
    x = log(ratios(useIdx));
    y = log(perims(useIdx));
    p = polyfit(x, y, 1);
    out.perimeterDecaySlope = p(1);
end
end
