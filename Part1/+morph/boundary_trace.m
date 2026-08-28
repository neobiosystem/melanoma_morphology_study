function C = boundary_trace(mask)
%BOUNDARY_TRACE Ordered closed contour of the largest component of a mask.
%
%   C = BOUNDARY_TRACE(mask) returns an N-by-2 array of [row, col] contour
%   vertices in traversal order, with the first vertex repeated as the last
%   so the polygon is explicitly closed.
%
%   Unlike MORPH.BOUNDARY_PIXELS (an unordered point set for box counting),
%   this ordered contour is what the perimeter-based descriptors need,
%   because they must smooth along the curve and integrate segment lengths.
%
%   Holes are ignored: a lesion is treated as a simply-connected region, so
%   only the outer boundary of the largest connected component is traced.
%
%   See also MORPH.SMOOTH_CONTOUR, MORPH.BOUNDARY_PIXELS.

B = bwboundaries(mask, 8, 'noholes');

if isempty(B)
    C = zeros(0, 2);
    return
end

% bwboundaries returns one contour per component; keep the largest lesion.
lengths       = cellfun(@(b) size(b, 1), B);
[~, largest]  = max(lengths);
C             = B{largest};

% Guarantee explicit closure regardless of what bwboundaries returned.
if ~isequal(C(1, :), C(end, :))
    C = [C; C(1, :)];
end
end
