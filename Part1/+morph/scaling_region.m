function out = scaling_region(epsVals, counts, cfg)
%SCALING_REGION Locate the genuine power-law range of a box-counting curve.
%
%   out = SCALING_REGION(epsVals, counts, cfg) finds the sub-range of box
%   sizes over which log N(eps) is genuinely linear in log(1/eps), and fits
%   the fractal dimension there.
%
%   epsVals must be ascending. The returned struct carries the fitted slope,
%   its standard error, R^2, the selected index range, the local-slope curve
%   and a QC flag.
%
%   METHOD -- local-slope plateau, not R^2 maximisation. A power law means
%   "locally constant slope in log-log", which is exactly what the local
%   slope measures. Maximising R^2 instead is a trap: R^2 rises trivially as
%   a window shrinks, so unconstrained R^2 search latches onto two- or
%   three-point windows and reports a confident slope fitted to noise. R^2
%   maximisation is therefore used only as an independent cross-check,
%   constrained to the same minimum width; if the two methods disagree by
%   more than cfg.scaling.crossCheckTol the lesion is flagged 'unstable'.
%
%   Truncation at either end biases the slope, which is why the caller bounds
%   eps by the raster floor below and the lesion scale above, and why this
%   function then narrows further within those bounds.
%
%   See also MORPH.FRACTAL_BOXCOUNT.

x = log(1 ./ epsVals(:));
y = log(counts(:));

out = struct('Df', NaN, 'DfStdErr', NaN, 'R2', NaN, ...
             'fitIdx', [], 'epsRange', [NaN NaN], 'numPoints', 0, ...
             'spanDecades', NaN, 'localSlopes', [], ...
             'method', 'local_slope_plateau', 'flag', 'insufficient_range', ...
             'DfCrossCheck', NaN);

% Drop scales where nothing was counted; log(0) is not usable.
valid = isfinite(x) & isfinite(y) & counts(:) > 0;
x = x(valid);  y = y(valid);
epsValid = epsVals(valid);
numScales = numel(x);

win = cfg.scaling.slopeWindow;
w2  = floor(win / 2);

if numScales < max(cfg.scaling.minPoints, win)
    % Too few usable scales to claim a scaling region; fit everything so a
    % number still exists, but keep the flag honest.
    if numScales >= 3
        [slope, ~, r2, se] = local_ols(x, y);
        out.Df = slope;  out.R2 = r2;  out.DfStdErr = se;
        out.fitIdx = 1:numScales;
        out.epsRange = [epsValid(1), epsValid(end)];
        out.numPoints = numScales;
        out.spanDecades = log10(epsValid(end) / epsValid(1));
    end
    return
end

% ---- Local slope at every interior scale -------------------------------
centres = (1 + w2) : (numScales - w2);
localSlopes = nan(numScales, 1);
for k = centres
    idx = (k - w2) : (k + w2);
    localSlopes(k) = local_ols(x(idx), y(idx));
end
out.localSlopes = localSlopes;

% ---- Longest plateau of near-constant local slope ----------------------
bestSupport = [];
bestScore   = [-Inf, -Inf];   % [numPoints, spanDecades]

for a = centres
    for b = a : centres(end)
        run = localSlopes(a:b);
        med = median(run);
        if any(abs(run - med) > cfg.scaling.slopeTol)
            % Do not break: the run median shifts as b grows, so a window
            % that violates the tolerance here can satisfy it once extended.
            continue
        end

        support = (a - w2) : (b + w2);
        nPts    = numel(support);
        span    = log10(epsValid(support(end)) / epsValid(support(1)));

        if nPts < cfg.scaling.minPoints || span < cfg.scaling.minDecades
            continue
        end

        score = [nPts, span];
        if score(1) > bestScore(1) || (score(1) == bestScore(1) && score(2) > bestScore(2))
            bestScore   = score;
            bestSupport = support;
        end
    end
end

if isempty(bestSupport)
    % No plateau wide enough. Fall back to a full-range fit so a number
    % still exists, but flag it distinctly from 'insufficient_range': this
    % lesion HAS enough scales, they simply never settle onto a power law,
    % so the slope carries end-truncation bias and must not be pooled
    % silently with well-behaved lesions.
    [slope, ~, r2, se] = local_ols(x, y);
    out.Df = slope;  out.R2 = r2;  out.DfStdErr = se;
    out.fitIdx = 1:numScales;
    out.epsRange = [epsValid(1), epsValid(end)];
    out.numPoints = numScales;
    out.spanDecades = log10(epsValid(end) / epsValid(1));
    out.flag = 'no_plateau';
    return
end

[slope, ~, r2, se] = local_ols(x(bestSupport), y(bestSupport));
out.Df          = slope;
out.R2          = r2;
out.DfStdErr    = se;
out.fitIdx      = bestSupport;
out.epsRange    = [epsValid(bestSupport(1)), epsValid(bestSupport(end))];
out.numPoints   = numel(bestSupport);
out.spanDecades = log10(out.epsRange(2) / out.epsRange(1));
out.flag        = 'ok';

% ---- Independent cross-check: best-R^2 window of comparable width ------
bestR2    = -Inf;
crossSlope = NaN;
for a = 1 : numScales
    for b = (a + cfg.scaling.minPoints - 1) : numScales
        span = log10(epsValid(b) / epsValid(a));
        if span < cfg.scaling.minDecades
            continue
        end
        [s, ~, r2c] = local_ols(x(a:b), y(a:b));
        if r2c > bestR2
            bestR2     = r2c;
            crossSlope = s;
        end
    end
end
out.DfCrossCheck = crossSlope;

if isfinite(crossSlope) && abs(crossSlope - out.Df) > cfg.scaling.crossCheckTol
    out.flag = 'unstable';
end
end

% -------------------------------------------------------------------------
function [slope, intercept, r2, seSlope] = local_ols(x, y)
%LOCAL_OLS Least-squares line fit returning slope, R^2 and slope std error.

n  = numel(x);
mx = mean(x);
my = mean(y);
sxx = sum((x - mx).^2);

if n < 2 || sxx == 0
    slope = NaN;  intercept = NaN;  r2 = NaN;  seSlope = NaN;
    return
end

slope     = sum((x - mx) .* (y - my)) / sxx;
intercept = my - slope * mx;

resid = y - (slope * x + intercept);
ssRes = sum(resid.^2);
ssTot = sum((y - my).^2);

if ssTot == 0
    r2 = NaN;
else
    r2 = 1 - ssRes / ssTot;
end

if n > 2
    seSlope = sqrt(ssRes / (n - 2) / sxx);
else
    seSlope = NaN;
end
end
