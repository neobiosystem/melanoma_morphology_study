function P = boundary_pixels(mask)
%BOUNDARY_PIXELS Thin 8-connected boundary point set of a binary mask.
%
%   P = BOUNDARY_PIXELS(mask) returns an N-by-2 array of [row, col] pixel
%   coordinates lying on the lesion boundary.
%
%   This is the point set fed to box counting, and it is deliberately the
%   RAW, UNSMOOTHED boundary. Smoothing before box counting would destroy
%   precisely the fine-scale irregularity that the fractal dimension is
%   meant to quantify, so no caller may smooth on this path. Perimeter-based
%   descriptors take the separate MORPH.BOUNDARY_TRACE route instead.
%
%   bwperim is called with connectivity 4, which marks a pixel as perimeter
%   when any of its 4-neighbours is background. That yields the THINNER,
%   8-connected boundary curve; connectivity 8 would produce a thicker
%   4-connected band and inflate box counts at small scales.
%
%   See also MORPH.FRACTAL_BOXCOUNT, MORPH.BOUNDARY_TRACE.

edgeMask = bwperim(mask, 4);
[r, c]   = find(edgeMask);
P        = [r, c];
end
