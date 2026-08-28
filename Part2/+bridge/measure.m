function out = measure(mask, cfg)
%MEASURE Score a simulated mask with Part 1's pipeline, keeping QC side-effects.
%
%   out = MEASURE(mask, cfg) calls MORPH.PROCESS_LESION(mask, cfg.measure) --
%   the unmodified Part 1 entry point -- but first records what that call is
%   about to silently hide:
%
%     * secondary connected components (satellite detachment) that Part 1
%       strips to the largest component, and the area fraction removed;
%     * holes (necrotic cores) that imfill is about to fill.
%
%   These are stored in out.preStrip. Without them the migration mechanism's
%   most distinctive phenotype -- cells detaching from the aggregate -- would
%   be measured as "no change", because the QC gate erases it before scoring.

mask = logical(mask);

pre = struct();
cc = bwconncomp(mask, 8);
pre.nComponents = cc.NumObjects;
sizes = cellfun(@numel, cc.PixelIdxList);
sizes = sort(sizes, 'descend');

pre.totalArea    = nnz(mask);
pre.largestArea  = sizes(1);
if numel(sizes) > 1
    pre.secondCcFraction     = sizes(2) / sizes(1);
    pre.strippedAreaFraction = (pre.totalArea - pre.largestArea) / pre.totalArea;
else
    pre.secondCcFraction     = 0;
    pre.strippedAreaFraction = 0;
end

filled = imfill(mask, 'holes');
pre.holeAreaFraction = (nnz(filled) - nnz(mask)) / max(1, nnz(mask));

out = morph.process_lesion(mask, cfg.measure);
out.preStrip = pre;
end
