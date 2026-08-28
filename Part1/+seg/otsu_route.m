function [mask, info] = otsu_route(img, cfg)
%OTSU_ROUTE Fully automatic lesion segmentation, the reproducibility arm.
%
%   [mask, info] = OTSU_ROUTE(img, cfg) runs the colour-space transform,
%   hair suppression, Otsu thresholding, morphological cleanup and optional
%   active-contour refinement described in the Part 1 design document.
%
%   This route exists to answer a specific question: how much of the
%   between-grade descriptor signal survives when no expert mask is
%   available? It is a sensitivity arm, not the primary analysis, because
%   segmentation quality plausibly degrades with severity -- ill-defined
%   infiltrative borders are objectively harder to threshold -- and an
%   Otsu-only comparison could therefore confuse segmentation difficulty
%   with biology.
%
%   Two dermoscopy-specific traps are handled explicitly:
%
%   * Vignetting. Most dermatoscopes image through a circular aperture, so
%     the frame corners are dark. Otsu happily labels those corners as
%     'lesion'. Components touching the image border are therefore discarded
%     when an interior component survives.
%
%   * Polarity. The lesion is darker than surrounding skin, so the class of
%     interest is BELOW the Otsu threshold, not above it.
%
%   See also SEG.REMOVE_HAIR, SEG.EXPERT_ROUTE.

info = struct('valid', false, 'reason', '', 'threshold', NaN, ...
              'hairFraction', NaN, 'usedActiveContour', false, ...
              'droppedBorderComponents', false);

if size(img, 3) == 3
    imgRgb = img;
else
    imgRgb = repmat(img, 1, 1, 3);
end

% ---- Hair suppression --------------------------------------------------
[imgClean, hairMask] = seg.remove_hair(imgRgb, cfg);
info.hairFraction = nnz(hairMask) / numel(hairMask);

% ---- Colour-space transform -------------------------------------------
switch cfg.seg.channel
    case 'Lab_L'
        lab  = rgb2lab(imgClean);
        chan = lab(:, :, 1) / 100;              % luminance, normalised
    case 'HSV_V'
        hsv  = rgb2hsv(imgClean);
        chan = hsv(:, :, 3);
    case 'HSV_H'
        hsv  = rgb2hsv(imgClean);
        chan = hsv(:, :, 1);
    otherwise
        error('seg:otsu_route:badChannel', ...
              'Unknown cfg.seg.channel: %s', cfg.seg.channel);
end

chan = mat2gray(chan);

% ---- Otsu threshold ----------------------------------------------------
% Dark vignette corners belong to neither skin nor lesion, but Otsu counts
% them in the dark class and so pulls the threshold down, which costs the
% lighter lesion periphery. Restricting the threshold's support to the
% aperture removes that pull; the threshold is still applied to the whole
% image so nothing is cropped away.
if cfg.seg.thresholdInsideAperture
    support = aperture_mask(size(chan));
    level   = graythresh(chan(support));
else
    level = graythresh(chan);
end
info.threshold = level;
bw = chan < level;              % lesion is the darker class

% ---- Morphological cleanup --------------------------------------------
minSide  = min(size(bw, 1), size(bw, 2));
closeRad = max(1, round(cfg.seg.closeRadiusRatio * minSide));

% Opening strips hairs and speckle but also erodes the lesion periphery, so
% a ratio of 0 disables it rather than quietly falling back to radius 1.
if cfg.seg.openRadiusRatio > 0
    openRad = max(1, round(cfg.seg.openRadiusRatio * minSide));
    bw = imopen(bw, strel('disk', openRad));
end
bw = imclose(bw, strel('disk', closeRad));   % close bubbles and gaps
bw = imfill(bw, 'holes');

if ~any(bw(:))
    info.reason = 'empty_after_cleanup';
    mask = bw;
    return
end

% ---- Discard vignette corners -----------------------------------------
% Dropping every border-touching component is wrong for a lesion that
% genuinely runs off the frame: if any small speck survives in the interior,
% the speck would be kept and the lesion thrown away. So the interior version
% is only preferred when its largest component is a substantial fraction of
% the largest component overall -- i.e. when what was removed really did look
% like vignette rather than the lesion itself.
interior = imclearborder(bw);
if any(interior(:)) && nnz(interior) < nnz(bw)
    if largest_component_size(interior) >= ...
            cfg.seg.minInteriorRetainRatio * largest_component_size(bw)
        bw = interior;
        info.droppedBorderComponents = true;
    end
end

% ---- Keep the single largest component --------------------------------
cc = bwconncomp(bw, 8);
if cc.NumObjects == 0
    info.reason = 'no_component';
    mask = false(size(bw));
    return
end
sizes = cellfun(@numel, cc.PixelIdxList);
[~, biggest] = max(sizes);

mask = false(size(bw));
mask(cc.PixelIdxList{biggest}) = true;

% ---- Optional active-contour refinement -------------------------------
% Initialised from the Otsu mask, this pulls the boundary onto the local
% intensity gradient and recovers detail that the morphological opening
% rounded off.
if cfg.seg.useActiveContour
    try
        refined = activecontour(chan, mask, cfg.seg.activeContourIter, 'Chan-Vese');
        if any(refined(:))
            refined = imfill(refined, 'holes');
            cc2 = bwconncomp(refined, 8);
            if cc2.NumObjects > 0
                s2 = cellfun(@numel, cc2.PixelIdxList);
                [~, b2] = max(s2);
                candidate = false(size(refined));
                candidate(cc2.PixelIdxList{b2}) = true;

                % Guard against the contour collapsing or exploding.
                ratio = nnz(candidate) / nnz(mask);
                if ratio > 0.5 && ratio < 2
                    mask = candidate;
                    info.usedActiveContour = true;
                end
            end
        end
    catch
        % Refinement is optional; keep the Otsu mask on failure.
    end
end

info.valid = any(mask(:));
if ~info.valid
    info.reason = 'empty_final_mask';
end
end

% =========================================================================
function m = aperture_mask(sz)
%APERTURE_MASK Inscribed ellipse approximating the dermatoscope aperture.
%
%   Deliberately geometric rather than intensity-derived: deciding which
%   pixels are vignette by their darkness would use the very quantity the
%   threshold is about to estimate.
[rr, cc] = ndgrid(1:sz(1), 1:sz(2));
m = ((rr - sz(1)/2) / (0.48 * sz(1))).^2 + ...
    ((cc - sz(2)/2) / (0.48 * sz(2))).^2 <= 1;
end

% =========================================================================
function n = largest_component_size(bw)
%LARGEST_COMPONENT_SIZE Pixel count of the biggest 8-connected component.
cc = bwconncomp(bw, 8);
if cc.NumObjects == 0
    n = 0;
else
    n = max(cellfun(@numel, cc.PixelIdxList));
end
end
