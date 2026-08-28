function S = step_apoptosis(S)
%STEP_APOPTOSIS Basal cell death plus starvation-driven necrosis.
%
%   Basal deaths are sampled as Poisson(delta * nCells) to avoid drawing a
%   full N^2 uniform matrix every step. Starvation deaths hit only the small
%   hypoxic set (C < C_min) with probability pStarve. Both are rare events;
%   the point is to clear nutrient-starved cores, whose imprint on the
%   boundary is invisible after Part 1's imfill -- recorded, not measured.

cfg = S.cfg;
occ = S.occ;

% Starvation / necrosis.
if cfg.sim.pStarve > 0 && ~isempty(S.Cfine)
    hyp = find(occ & S.Cfine < cfg.sim.C_min);
    if ~isempty(hyp)
        die = hyp(rand(numel(hyp), 1) < cfg.sim.pStarve);
        if ~isempty(die)
            occ(die) = false;
            S.totalDeaths = S.totalDeaths + numel(die);
        end
    end
end

% Basal death.
if cfg.sim.delta > 0
    nCells = nnz(occ);
    k = poissrnd(cfg.sim.delta * nCells);
    if k > 0
        idx = find(occ);
        k = min(k, numel(idx));
        die = idx(randperm(numel(idx), k));
        occ(die) = false;
        S.totalDeaths = S.totalDeaths + numel(die);
    end
end

S.occ = occ;
S.age(~occ)   = 0;
S.clone(~occ) = 0;
end
