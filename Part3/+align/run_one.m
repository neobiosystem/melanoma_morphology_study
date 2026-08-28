function row = run_one(r, cfg0, descs)
%RUN_ONE Grow, measure, and collapse one Part 3 alignment run to a row.
%
%   row = RUN_ONE(r, cfg0, descs) applies the arm base and the level override
%   (all flat cfg.sim fields -- the Part 3 knobs need no nesting), grows the
%   tumour with CA3.RUN_SIMULATION, scores it clean and degraded via the Part
%   1 bridge, and returns a one-row table.
%
%   Lives in a package so it is visible to parfor workers.

cfg = cfg0;
cfg.sim.seed = r.seed;

b = align.arm_base(r.arm);
for f = fieldnames(b)'
    cfg.sim.(f{1}) = b.(f{1});
end

for f = fieldnames(r.overrides)'
    cfg.sim.(f{1}) = r.overrides.(f{1});
end

S = ca3.run_simulation(cfg);
mask = ca.to_mask(S);
clean = bridge.measure(mask, cfg);
degraded = bridge.measure(bridge.degrade(mask, cfg), cfg);

row = table();
row.arm   = string(r.arm);
row.level = r.level;
row.seed  = r.seed;
row.nCells = S.nCells;
row.dEq    = 2 * sqrt(S.nCells / pi);
row.steps  = S.step;
row.accepted = clean.accepted;
row.DfFlag   = string(clean.DfFlag);
row.strippedAreaFraction = clean.preStrip.strippedAreaFraction;
row.holeAreaFraction     = clean.preStrip.holeAreaFraction;
row.secondCcFraction     = clean.preStrip.secondCcFraction;
for d = descs
    nm = d{1};
    row.(['clean_'    nm]) = clean.(nm);
    row.(['degraded_' nm]) = degraded.(nm);
end
end
