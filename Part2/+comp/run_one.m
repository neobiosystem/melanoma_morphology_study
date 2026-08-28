function row = run_one(r, cfg0, descs)
%RUN_ONE Grow, measure, and collapse one mechanism-comparison run to a row.
%
%   row = RUN_ONE(r, cfg0, descs) applies the arm base and the level override,
%   grows the tumour, scores it clean and degraded via the Part 1 bridge, and
%   returns a one-row table.
%
%   This lives in a package (not a script-local function) so it is visible to
%   parfor workers.

cfg = cfg0;
cfg.sim.seed = r.seed;

b = comp.arm_base(r.arm);

% Arm base: flat fields directly, subclone_* into the nested struct.
for f = fieldnames(b)'
    if startsWith(f{1}, 'subclone_'), continue; end
    cfg.sim.(f{1}) = b.(f{1});
end
cfg.sim.subclone.frac      = b.subclone_frac;
cfg.sim.subclone.alphaMult = b.subclone_alphaMult;

% Level override, mapping the subclone knobs onto the nested struct.
for f = fieldnames(r.overrides)'
    nm = f{1};
    if strcmp(nm, 'subclone_frac')
        cfg.sim.subclone.frac = r.overrides.(nm);
    elseif strcmp(nm, 'subclone_alphaMult')
        cfg.sim.subclone.alphaMult = r.overrides.(nm);
    else
        cfg.sim.(nm) = r.overrides.(nm);
    end
end

S = ca.run_simulation(cfg);
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
