function [mask, info] = expert_route(maskPath)
%EXPERT_ROUTE Load and validate an expert ground-truth segmentation mask.
%
%   [mask, info] = EXPERT_ROUTE(maskPath) reads a challenge ground-truth PNG
%   and returns a logical mask plus a small info struct.
%
%   This is the PRIMARY analysis route. Expert masks are the reference
%   standard and the descriptors measured from them become Part 3's
%   calibration targets; the Otsu route exists to quantify how much
%   descriptor error a mask-less pipeline would incur, not to replace this.
%
%   Challenge masks are 0/255 binary PNGs, so a mid-range threshold is
%   sufficient and robust to any incidental anti-aliasing.
%
%   See also SEG.OTSU_ROUTE, MORPH.PROCESS_LESION.

info = struct('valid', false, 'reason', '', 'fillFraction', NaN);

raw = imread(maskPath);
if size(raw, 3) > 1
    raw = raw(:, :, 1);
end

mask = raw > 127;

info.fillFraction = nnz(mask) / numel(mask);

if ~any(mask(:))
    info.reason = 'empty_expert_mask';
    return
end

info.valid = true;
end
