function [imgOut, hairMask] = remove_hair(img, cfg)
%REMOVE_HAIR Suppress dark hair shafts in a dermoscopic image (DullRazor).
%
%   [imgOut, hairMask] = REMOVE_HAIR(img, cfg) detects thin dark curvilinear
%   structures and replaces them with locally interpolated skin.
%
%   Hair is the dominant artefact in dermoscopy and it attaches spurious
%   filaments to the lesion boundary. Because the fractal dimension responds
%   precisely to fine boundary structure, uncorrected hair inflates D_f in a
%   way that mimics infiltrative growth -- the exact signal the study is
%   trying to measure. It therefore has to be removed before thresholding.
%
%   Detection is a bank of morphological closings with line structuring
%   elements at several orientations: closing with a line longer than the
%   hair width removes the shaft, so the closed-minus-original residual is
%   large exactly on hair.
%
%   Filling uses the closed image values rather than a Laplace inpaint. That
%   is the original DullRazor approximation, it is what the closing already
%   computed, and it avoids a per-image PDE solve that would dominate
%   runtime at native dermoscopic resolutions.
%
%   See also SEG.OTSU_ROUTE.

if size(img, 3) == 3
    gray = rgb2gray(img);
else
    gray = img;
end

minSide = min(size(gray, 1), size(gray, 2));
seLen   = max(5, round(cfg.seg.hairSeLengthRatio * minSide));

angles  = linspace(0, 180, cfg.seg.hairNumOrient + 1);
angles(end) = [];

% Maximum response across orientations: hair at any angle is closed away.
closedMax = gray;
for a = angles
    closedMax = max(closedMax, imclose(gray, strel('line', seLen, a)));
end

residual = double(closedMax) - double(gray);
hairMask = residual > cfg.seg.hairThreshold;

% Hair is elongated; drop isolated speckle that is more likely to be
% pigment network or bubbles than a shaft.
hairMask = bwareaopen(hairMask, max(10, round(seLen / 2)));
hairMask = imdilate(hairMask, strel('disk', 2));

imgOut = img;
if ~any(hairMask(:))
    return
end

for c = 1:size(img, 3)
    channel   = img(:, :, c);
    closedCh  = imclose(channel, strel('disk', max(3, round(seLen / 3))));
    channel(hairMask) = closedCh(hairMask);
    imgOut(:, :, c) = channel;
end
end
