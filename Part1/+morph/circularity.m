function [circ, radialCv] = circularity(P, centroidRC)
%CIRCULARITY Radial-deviation circularity, 1 - std(r)/mean(r).
%
%   [circ, radialCv] = CIRCULARITY(P, centroidRC) measures how far a
%   boundary departs from an equal-area circle, using the radii r from the
%   region centroid to each boundary pixel in P (N-by-2 [row col]).
%   centroidRC is the region centroid as [row col].
%
%   circ is 1 for a perfect circle and falls as the boundary becomes more
%   radially anisotropic. radialCv is the underlying coefficient of
%   variation std(r)/mean(r), returned for diagnostics.
%
%   WHY NOT 4*pi*A/P^2 -- the source design document lists both compactness
%   and circularity as descriptors, but that naive circularity is exactly
%   1/compactness, so the pair would carry zero independent information and
%   would inflate the apparent descriptor count. This radial definition
%   matches the document's own wording ("deviation from an equal-area
%   circle") while measuring genuinely different geometry: compactness is a
%   global perimeter-to-area efficiency, circularity is radial anisotropy.
%   An ellipse has low circularity but only mildly raised compactness; an
%   angularly uniform but crenellated shape is the reverse.
%
%   Radii are taken from the RAW boundary because a mean-and-spread
%   statistic over many radii is insensitive to per-pixel staircase jitter,
%   unlike a path-length integral.
%
%   See also MORPH.COMPACTNESS.

if isempty(P)
    circ = NaN;  radialCv = NaN;
    return
end

dr = double(P(:, 1)) - centroidRC(1);
dc = double(P(:, 2)) - centroidRC(2);
r  = hypot(dr, dc);

meanR = mean(r);
if meanR <= 0
    circ = NaN;  radialCv = NaN;
    return
end

radialCv = std(r) / meanR;
circ     = 1 - radialCv;
end
