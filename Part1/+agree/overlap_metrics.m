function out = overlap_metrics(maskA, maskB)
%OVERLAP_METRICS Spatial agreement between two segmentations of one lesion.
%
%   out = OVERLAP_METRICS(maskA, maskB) returns Dice, Jaccard, the area
%   ratio and a Hausdorff distance normalised by equivalent diameter.
%
%   maskA is the reference (expert) and maskB the comparison (Otsu). The
%   area ratio is signed in the sense that values below 1 mean the automatic
%   route under-segmented, which is the direction that matters: systematic
%   under-segmentation would shrink measured lesions and shift every
%   size-dependent descriptor at once.
%
%   Hausdorff is normalised because the cohort spans a ~10x resolution
%   range, and an unnormalised distance in pixels would simply rank images
%   by how big they are.

out = struct('dice', NaN, 'jaccard', NaN, 'areaRatio', NaN, ...
             'hausdorffNorm', NaN, 'valid', false);

maskA = logical(maskA);
maskB = logical(maskB);

if ~isequal(size(maskA), size(maskB))
    return
end

nA = nnz(maskA);
nB = nnz(maskB);
if nA == 0 || nB == 0
    return
end

inter = nnz(maskA & maskB);
union = nnz(maskA | maskB);

out.dice      = 2 * inter / (nA + nB);
out.jaccard   = inter / union;
out.areaRatio = nB / nA;

% Symmetric Hausdorff over boundary pixels, normalised by the reference dEq.
try
    pa = morph.boundary_pixels(maskA);
    pb = morph.boundary_pixels(maskB);
    if ~isempty(pa) && ~isempty(pb)
        dAB = max(min(pdist2(pa, pb), [], 2));
        dBA = max(min(pdist2(pb, pa), [], 2));
        dEq = morph.equivalent_diameter(maskA);
        if dEq > 0
            out.hausdorffNorm = max(dAB, dBA) / dEq;
        end
    end
catch
    % pdist2 can exhaust memory on very large boundaries; leave as NaN.
end

out.valid = true;
end
