function out = process_lesion(mask, cfg)
%PROCESS_LESION Compute every morphological descriptor for one lesion mask.
%
%   out = PROCESS_LESION(mask, cfg) applies the quality-control gates and
%   returns a flat struct holding all descriptors plus QC bookkeeping.
%
%   This function is the SINGLE entry point for both segmentation routes.
%   The expert-mask route and the Otsu route differ only in where the mask
%   comes from; everything downstream is byte-identical. That is deliberate:
%   if the two routes used even slightly different boundary extraction or
%   smoothing, the mismatch would show up in the agreement analysis
%   disguised as genuine segmentation disagreement.
%
%   Key fields:
%       accepted      - false when a hard QC gate rejected the lesion
%       rejectReason  - '' when accepted, otherwise the failing gate
%       flags         - cell array of non-fatal warnings
%       Df, R2, ...   - descriptors (NaN when not computable)
%
%   See also MORPH.FRACTAL_BOXCOUNT, SEG.EXPERT_ROUTE, SEG.OTSU_ROUTE.

out = init_output();

mask = logical(mask);
if ~any(mask(:))
    out.rejectReason = 'empty_mask';
    return
end

% ---- Connected components ---------------------------------------------
cc = bwconncomp(mask, 8);
if cc.NumObjects > 1
    sizes = cellfun(@numel, cc.PixelIdxList);
    [sorted, order] = sort(sizes, 'descend');
    secondFrac = sorted(2) / sorted(1);

    if secondFrac > cfg.qc.rejectSecondCcFrac
        out.rejectReason = 'multiple_components';
        out.secondCcFraction = secondFrac;
        return
    end

    if secondFrac > cfg.qc.flagSecondCcFrac
        out.flags{end+1} = 'minor_second_component';
    end

    keep = false(size(mask));
    keep(cc.PixelIdxList{order(1)}) = true;
    mask = keep;
    out.secondCcFraction = secondFrac;
end

mask = imfill(mask, 'holes');

% ---- Framing gates -----------------------------------------------------
borderFrac = border_contact_fraction(mask);
out.borderContactFraction = borderFrac;
if cfg.qc.rejectBorderTouch && borderFrac > cfg.qc.maxBorderContactFrac
    % Substantial truncation: perimeter and area no longer describe the
    % lesion. A light graze is kept and flagged instead, because rejecting
    % on any contact at all would preferentially discard large lesions.
    out.rejectReason = 'truncated_at_border';
    return
elseif borderFrac > 0
    out.flags{end+1} = 'grazes_border';
end

frameFraction = nnz(mask) / numel(mask);
out.frameFraction = frameFraction;
if frameFraction > cfg.qc.maxFrameFraction
    out.rejectReason = 'fills_frame';
    return
end

dEq = morph.equivalent_diameter(mask);
out.eqDiameterNative = dEq;

if dEq < cfg.qc.minEqDiameterPx
    out.rejectReason = 'lesion_too_small';
    return
end
if dEq > cfg.qc.maxEqDiameterPx
    out.rejectReason = 'lesion_too_large';
    return
end

% ---- Compute cap -------------------------------------------------------
% Purely a cost control. Masks are otherwise left at native resolution:
% resampling every lesion to a common size would erase the fine boundary
% detail under study, and would do so worst on the highest-resolution
% images, which carry the most signal.
if dEq > cfg.qc.computeCapPx
    scale = cfg.qc.computeCapPx / dEq;
    mask  = imresize(mask, scale, 'nearest');
    mask  = imfill(logical(mask), 'holes');
    dEq   = morph.equivalent_diameter(mask);
    out.downsampled = true;
    out.flags{end+1} = 'compute_cap_downsampled';

    if ~any(mask(:)) || dEq < cfg.qc.minEqDiameterPx
        out.rejectReason = 'lesion_too_small_after_cap';
        return
    end
end

out.eqDiameter = dEq;
out.area       = nnz(mask);

% ---- Fractal dimension on the RAW boundary ----------------------------
P = morph.boundary_pixels(mask);
out.numBoundaryPixels = size(P, 1);

box = morph.fractal_boxcount(P, dEq, cfg);
out.Df             = box.Df;
out.DfStdErr       = box.DfStdErr;
out.DfR2           = box.R2;
out.DfSpanDecades  = box.spanDecades;
out.DfNumPoints    = box.numPoints;
out.DfFlag         = box.flag;
out.DfCrossCheck   = box.DfCrossCheck;
out.boxcount       = box;

switch box.flag
    case 'insufficient_range'
        out.flags{end+1} = 'boxcount_insufficient_range';
    case 'no_plateau'
        out.flags{end+1} = 'boxcount_no_plateau';
    case 'unstable'
        out.flags{end+1} = 'boxcount_unstable';
end
if isfinite(box.R2) && box.R2 < cfg.qc.minBoxcountR2
    out.flags{end+1} = 'boxcount_low_r2';
end
if isfinite(box.DfStdErr) && box.DfStdErr > cfg.qc.maxSlopeStdErr
    out.flags{end+1} = 'boxcount_high_stderr';
end

% ---- Perimeter-based descriptors on the SMOOTHED contour --------------
C = morph.boundary_trace(mask);
primaryWin = round(cfg.smooth.primaryRatio * dEq);
Csmooth    = morph.smooth_contour(C, primaryWin);

out.perimeterRaw      = morph.contour_perimeter(C);
out.perimeterSmoothed = morph.contour_perimeter(Csmooth);
out.compactness       = morph.compactness(Csmooth);
out.compactnessRaw    = morph.compactness(C);   % kept only to quantify the bias

% ---- Radial and symmetry descriptors on the raw mask ------------------
props = regionprops(mask, 'Centroid');
centroidRC = [props(1).Centroid(2), props(1).Centroid(1)];   % [row col]
out.centroidRC = centroidRC;

[out.circularity, out.radialCv] = morph.circularity(P, centroidRC);

asym            = morph.asymmetry_index(mask);
out.asymmetry   = asym.mean;
out.asymmetryMajor = asym.major;
out.asymmetryMinor = asym.minor;
out.asymmetryCentral = asym.central;
out.elongation  = asym.elongation;

ms                       = morph.multi_scale(mask, dEq, cfg);
out.solidity             = ms.solidity;
out.convexity            = ms.convexity;
out.perimeterDecaySlope  = ms.perimeterDecaySlope;
out.perimeterCurve       = ms.perimeters;
out.perimeterRatioCurve  = ms.perimeterRatio;
out.sigmaPx              = ms.sigmaPx;

% ---- Numerical plausibility -------------------------------------------
if isfinite(out.compactness) && (out.compactness < 1 || out.compactness > 100)
    out.flags{end+1} = 'implausible_compactness';
end
if isfinite(out.asymmetry) && (out.asymmetry < 0 || out.asymmetry > 1)
    out.flags{end+1} = 'implausible_asymmetry';
end

out.accepted = true;
end

% -------------------------------------------------------------------------
function frac = border_contact_fraction(mask)
%BORDER_CONTACT_FRACTION Share of the lesion boundary lying on the frame edge.
edge = bwperim(mask, 4);
total = nnz(edge);
if total == 0
    frac = 0;
    return
end
onFrame = nnz(edge(1, :)) + nnz(edge(end, :)) + ...
          nnz(edge(:, 1)) + nnz(edge(:, end));
frac = onFrame / total;
end

% -------------------------------------------------------------------------
function out = init_output()
out = struct( ...
    'accepted', false, 'rejectReason', '', 'flags', {{}}, ...
    'area', NaN, 'eqDiameter', NaN, 'eqDiameterNative', NaN, ...
    'frameFraction', NaN, 'secondCcFraction', 0, 'downsampled', false, ...
    'borderContactFraction', NaN, ...
    'numBoundaryPixels', NaN, 'centroidRC', [NaN NaN], ...
    'Df', NaN, 'DfStdErr', NaN, 'DfR2', NaN, 'DfSpanDecades', NaN, ...
    'DfNumPoints', NaN, 'DfFlag', '', 'DfCrossCheck', NaN, 'boxcount', struct(), ...
    'perimeterRaw', NaN, 'perimeterSmoothed', NaN, ...
    'compactness', NaN, 'compactnessRaw', NaN, ...
    'circularity', NaN, 'radialCv', NaN, ...
    'asymmetry', NaN, 'asymmetryMajor', NaN, 'asymmetryMinor', NaN, ...
    'asymmetryCentral', NaN, 'elongation', NaN, ...
    'solidity', NaN, 'convexity', NaN, 'perimeterDecaySlope', NaN, ...
    'perimeterCurve', [], 'perimeterRatioCurve', [], 'sigmaPx', []);
end
