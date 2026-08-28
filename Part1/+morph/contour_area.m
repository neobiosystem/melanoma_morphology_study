function A = contour_area(C)
%CONTOUR_AREA Polygon area enclosed by a closed contour (shoelace formula).
%
%   A = CONTOUR_AREA(C) returns the area enclosed by the closed contour C
%   (N-by-2, first vertex repeated last). The result is unsigned, so
%   traversal direction does not matter.
%
%   Use together with MORPH.CONTOUR_PERIMETER so that compactness is
%   computed from a single consistent geometry.

if size(C, 1) < 3
    A = 0;
    return
end

% polyarea expects an open vertex list.
if isequal(C(1, :), C(end, :))
    C(end, :) = [];
end

A = polyarea(C(:, 1), C(:, 2));
end
