function K = compactness(C)
%COMPACTNESS Perimeter-to-area efficiency of a closed contour, P^2/(4*pi*A).
%
%   K = COMPACTNESS(C) returns 1 for a perfect circle and grows without
%   bound as the contour becomes more convoluted.
%
%   Pass a SMOOTHED contour (see MORPH.SMOOTH_CONTOUR). Computing this on a
%   raw raster boundary would import the staircase perimeter bias, which is
%   larger for irregular boundaries than for smooth ones and would therefore
%   fake a between-grade difference. Perimeter and area are both taken from
%   the contour passed in, so the two come from one consistent geometry.
%
%   See also MORPH.CIRCULARITY, MORPH.SMOOTH_CONTOUR.

P = morph.contour_perimeter(C);
A = morph.contour_area(C);

if A <= 0
    K = NaN;
    return
end

K = P^2 / (4 * pi * A);
end
