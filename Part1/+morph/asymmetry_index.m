function out = asymmetry_index(mask)
%ASYMMETRY_INDEX Principal-axis and central reflection asymmetry of a mask.
%
%   out = ASYMMETRY_INDEX(mask) returns a struct with fields:
%       major      - non-overlap fraction when folded about the major axis
%       minor      - non-overlap fraction when folded about the minor axis
%       mean       - (major + minor)/2, the reported asymmetry index
%       max        - max(major, minor)
%       central    - non-overlap fraction under 180 deg point reflection
%       elongation - (a - b)/(a + b) from the principal axis lengths
%
%   All non-overlap fractions lie in [0, 1] and are 0 for a perfectly
%   symmetric shape.
%
%   'mean' is the primary index and corresponds to the A of the clinical
%   ABCD rule, matching the source design document. 'central' is a
%   robustness backup: point reflection needs no axis orientation, so it is
%   immune to the orientation ambiguity that appears when a lesion is nearly
%   circular and the principal axis is poorly determined.
%
%   elongation is reported alongside deliberately. It is second-moment shape
%   only, so a symmetric ellipse has elongation > 0 while its asymmetry is
%   ~0 -- the pair is independent, and the test suite asserts exactly that.

out = struct('major', NaN, 'minor', NaN, 'mean', NaN, 'max', NaN, ...
             'central', NaN, 'elongation', NaN);

area = nnz(mask);
if area == 0
    return
end

props = regionprops(mask, 'Centroid', 'Orientation', ...
                    'MajorAxisLength', 'MinorAxisLength');
if isempty(props)
    return
end

% regionprops may return several components; the caller is expected to have
% reduced to the largest, but guard anyway.
if numel(props) > 1
    stats = regionprops(mask, 'Area');
    [~, keep] = max([stats.Area]);
    props = props(keep);
end

a = props.MajorAxisLength;
b = props.MinorAxisLength;
if (a + b) > 0
    out.elongation = (a - b) / (a + b);
end

% ---- Re-centre and align the lesion to its principal axes --------------
% A square canvas large enough for any rotation, with the centroid placed
% at the exact centre so that flips are about the true principal axes.
[h, w]   = size(mask);
canvasSz = 2 * ceil(hypot(h, w) / 2) + 4;
centre   = canvasSz / 2;

centroidXY = props.Centroid;          % [x y] == [col row]
rowShift   = round(centre - centroidXY(2));
colShift   = round(centre - centroidXY(1));

canvas = false(canvasSz, canvasSz);
rows   = (1:h) + rowShift;
cols   = (1:w) + colShift;
inR    = rows >= 1 & rows <= canvasSz;
inC    = cols >= 1 & cols <= canvasSz;
canvas(rows(inR), cols(inC)) = mask(inR, inC);

% Orientation is measured counter-clockwise from the x axis, so rotating by
% -orientation lays the major axis along the rows.
aligned = imrotate(canvas, -props.Orientation, 'nearest', 'crop');

alignedArea = nnz(aligned);
if alignedArea == 0
    return
end

out.major   = nnz(xor(aligned, flipud(aligned))) / (2 * alignedArea);
out.minor   = nnz(xor(aligned, fliplr(aligned))) / (2 * alignedArea);
out.mean    = (out.major + out.minor) / 2;
out.max     = max(out.major, out.minor);

% Point reflection about the centroid; no alignment needed.
centredArea = nnz(canvas);
if centredArea > 0
    out.central = nnz(xor(canvas, rot90(canvas, 2))) / (2 * centredArea);
end
end
