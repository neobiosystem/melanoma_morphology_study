function S = init_state(cfg)
%INIT_STATE Part 2 initial state plus the growth-rate field (arm B).
%
%   S = INIT_STATE(cfg) reuses Part 2's CA.INIT_STATE verbatim (never copy the
%   seed/occupancy logic) and adds S.growthField, the smooth static growth-rate
%   field used by the spatial-heterogeneity arm. With heteroOn off the field
%   is all-ones so CA3.STEP_PROLIFERATE reduces exactly to Part 2.

S = ca.init_state(cfg);

N = cfg.sim.N;
if isfield(cfg.sim, 'heteroOn') && cfg.sim.heteroOn
    lambda = cfg.sim.heteroLambda;
    gAmp   = cfg.sim.heteroAmp;
    S.growthField = mech.growth_field(N, lambda, gAmp, cfg.sim.seed + 70000);
else
    S.growthField = ones(N, N);
end
end
