function P = contour_perimeter(C)
%CONTOUR_PERIMETER Polyline length of a closed contour.
%
%   P = CONTOUR_PERIMETER(C) sums the Euclidean lengths of the segments of
%   the closed contour C (N-by-2, first vertex repeated last).
%
%   Pair this with MORPH.CONTOUR_AREA on the SAME contour. Mixing a smoothed
%   perimeter with a raster pixel-count area would make compactness
%   internally inconsistent, since the two quantities would describe
%   different geometries.

if size(C, 1) < 2
    P = 0;
    return
end

if ~isequal(C(1, :), C(end, :))
    C = [C; C(1, :)];
end

d = diff(C, 1, 1);
P = sum(hypot(d(:, 1), d(:, 2)));
end
