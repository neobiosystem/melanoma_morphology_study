function maskDeg = degrade(mask, cfg)
%DEGRADE First-order Otsu-like segmentation noise on the boundary.
%
%   maskDeg = DEGRADE(mask, cfg) perturbs a fraction of boundary pixels by
%   +/-1 px (stochastic erosion/dilation) and refills holes. This is a
%   deliberately simple model of the segmentation roughness a real Otsu mask
%   carries, not a claim to reproduce Otsu: it exists so every simulated mask
%   can be measured BOTH clean and degraded, and Part 3 can decide which
%   comparison to trust. Main conclusions use between-group differences,
%   which this symmetric noise mostly cancels.

% Off the simulation stream so it never disturbs run determinism.
rng(cfg.measure.stats.seed + 9000, 'twister');

frac = cfg.degrade.frac;
se   = strel('disk', 1, 0);

ero = imerode(mask, se);
dil = imdilate(mask, se);
r   = rand(size(mask));

maskDeg = mask;
maskDeg(ero ~= mask & r < frac)      = true;    % grow some boundary pixels
maskDeg(dil ~= mask & r > 1 - frac)  = false;   % shrink some boundary pixels
maskDeg = imfill(maskDeg, 'holes');
end
