function out = fractal_boxcount(P, dEq, cfg)
%FRACTAL_BOXCOUNT Box-counting fractal dimension of a boundary point set.
%
%   out = FRACTAL_BOXCOUNT(P, dEq, cfg) covers the N-by-2 boundary point set
%   P with square grids of side eps, counts the occupied boxes N(eps), and
%   fits the slope of log N(eps) against log(1/eps) over an automatically
%   detected scaling region. dEq is the lesion's equivalent diameter.
%
%   P must be the RAW boundary point set from MORPH.BOUNDARY_PIXELS. Feeding
%   a smoothed contour here would erase the very irregularity being measured.
%
%   Three implementation choices carry most of the correctness:
%
%   1. Counting primitive. A sparse ACCUMARRAY 2-D histogram counts distinct
%      occupied cells in O(n) and composes cleanly with the offset loop.
%
%   2. Grid origin. Box counts oscillate with grid placement, and that
%      oscillation biases the fitted slope. The box dimension is defined
%      through the INFIMUM over covers, so N(eps) is minimised over a
%      q-by-q lattice of sub-box origin shifts. Averaging over origins is a
%      common shortcut but estimates the mean cover, which is a different
%      quantity; min is the principled choice.
%
%   3. Scale bounds. eps is log-spaced -- a linear sequence would cluster
%      points at small eps where digital noise dominates and would skew the
%      least-squares slope. Bounds are tied to dEq: the lower bound also
%      respects a hard pixel floor, below which a digital boundary is
%      locally just a straight line and D_f collapses toward 1; the upper
%      bound stops before N(eps) saturates at the lesion scale.
%
%   See also MORPH.SCALING_REGION, MORPH.BOUNDARY_PIXELS.

out = struct('Df', NaN, 'DfStdErr', NaN, 'R2', NaN, ...
             'epsAll', [], 'countsAll', [], 'epsRange', [NaN NaN], ...
             'numPoints', 0, 'spanDecades', NaN, 'localSlopes', [], ...
             'fitIdx', [], 'numOffsets', cfg.boxcount.numOffsets, ...
             'method', 'local_slope_plateau', 'flag', 'insufficient_range', ...
             'DfCrossCheck', NaN);

if isempty(P) || size(P, 1) < 8 || ~isfinite(dEq) || dEq <= 0
    return
end

epsMin = max(cfg.boxcount.minBoxPx, cfg.boxcount.minBoxRatio * dEq);
epsMax = cfg.boxcount.maxBoxRatio * dEq;

% Need a real span to fit a slope through.
if epsMax <= epsMin * 1.5
    return
end

epsVals = logspace(log10(epsMin), log10(epsMax), cfg.boxcount.numScales)';

% Shift to a non-negative origin so grid indices stay >= 1 under any offset.
P0 = double(P) - min(double(P), [], 1);

q      = cfg.boxcount.numOffsets;
counts = zeros(numel(epsVals), 1);

for k = 1:numel(epsVals)
    e    = epsVals(k);
    best = Inf;

    for a = 0:(q - 1)
        for b = 0:(q - 1)
            % Shifting the points by +d is equivalent to shifting the grid
            % origin by -d, and keeps every index positive.
            dr = a * e / q;
            dc = b * e / q;

            i = floor((P0(:, 1) + dr) / e) + 1;
            j = floor((P0(:, 2) + dc) / e) + 1;

            occupied = accumarray([i, j], 1, [], [], [], true);
            best     = min(best, nnz(occupied));
        end
    end

    counts(k) = best;
end

fit = morph.scaling_region(epsVals, counts, cfg);

out.Df           = fit.Df;
out.DfStdErr     = fit.DfStdErr;
out.R2           = fit.R2;
out.epsAll       = epsVals;
out.countsAll    = counts;
out.epsRange     = fit.epsRange;
out.numPoints    = fit.numPoints;
out.spanDecades  = fit.spanDecades;
out.localSlopes  = fit.localSlopes;
out.fitIdx       = fit.fitIdx;
out.method       = fit.method;
out.flag         = fit.flag;
out.DfCrossCheck = fit.DfCrossCheck;
end
