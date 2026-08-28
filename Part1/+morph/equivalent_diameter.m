function dEq = equivalent_diameter(mask)
%EQUIVALENT_DIAMETER Diameter of the circle with the same area as the mask.
%
%   dEq = EQUIVALENT_DIAMETER(mask) returns 2*sqrt(A/pi) in pixels, where A
%   is the number of true pixels in the binary MASK.
%
%   dEq is the reference length for every scale-dependent setting in Part 1:
%   box-counting bounds and contour-smoothing windows are expressed as
%   fractions of it so that measurements stay comparable across the ~10x
%   spread of image resolutions in the ISIC cohort.

area = nnz(mask);
dEq  = 2 * sqrt(area / pi);
end
