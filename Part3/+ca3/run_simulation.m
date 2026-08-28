function S = run_simulation(cfg)
%RUN_SIMULATION Orchestrate one hybrid-CA growth run (Part 3 extension).
%
%   S = RUN_SIMULATION(cfg) is Part 2's CA.RUN_SIMULATION with exactly two
%   substitutions: the initial state comes from CA3.INIT_STATE (which adds the
%   growth-rate field) and proliferation from CA3.STEP_PROLIFERATE (which adds
%   push-displacement and the field modulation). Every other step
%   (nutrient, migration, apoptosis) is reused verbatim from Part 2.
%
%   Determinism: same cfg.sim.seed gives a byte-identical lattice.

rng(cfg.sim.seed, 'twister');

S = ca3.init_state(cfg);

% With nutrient coupling off the field is uniformly C0, so proliferation is
% not nutrient-limited (arms isolate their single mechanism).
if cfg.sim.nutrientCoupling
    S = ca.step_nutrient(S);
else
    S.Cfine = cfg.sim.C0 * ones(cfg.sim.N, cfg.sim.N);
end

% Metric history buffers.
nSnap = ceil(cfg.sim.maxSteps / cfg.sim.snapshotInterval) + 1;
histStep  = zeros(nSnap, 1);
histN     = zeros(nSnap, 1);
histDEq   = zeros(nSnap, 1);
histMinC  = nan(nSnap, 1);
histMaxC  = nan(nSnap, 1);
histDiv   = zeros(nSnap, 1);
histDeath = zeros(nSnap, 1);
nRec = 0;

maskSnapshots = struct('step', {}, 'mask', {});

for t = 1:cfg.sim.maxSteps
    S.step = t;

    S.age(S.occ) = S.age(S.occ) + 1;                 % advance maturity

    if cfg.sim.nutrientCoupling && mod(t, cfg.sim.nutrientInterval) == 0
        S = ca.step_nutrient(S);
    end

    S = ca3.step_proliferate(S);
    S = ca.step_migrate(S);
    S = ca.step_apoptosis(S);

    S.nCells = nnz(S.occ);

    if mod(t, cfg.sim.snapshotInterval) == 0
        dEq = 2 * sqrt(S.nCells / pi);
        nRec = nRec + 1;
        histStep(nRec)  = t;
        histN(nRec)     = S.nCells;
        histDEq(nRec)   = dEq;
        histDiv(nRec)   = S.totalDivisions;
        histDeath(nRec) = S.totalDeaths;
        if cfg.sim.nutrientCoupling
            histMinC(nRec) = min(S.Cfine(:));
            histMaxC(nRec) = max(S.Cfine(:));
        end
        if mod(t, cfg.sim.maskSaveInterval) == 0
            maskSnapshots(end+1) = struct('step', t, 'mask', S.occ); %#ok<AGROW>
        end
        if dEq >= cfg.sim.stopEqDiameter
            break
        end
    end
end

% Final snapshot is always retained for the measurement step.
maskSnapshots(end+1) = struct('step', S.step, 'mask', S.occ);

S.history = table(histStep(1:nRec), histN(1:nRec), histDEq(1:nRec), ...
    histMinC(1:nRec), histMaxC(1:nRec), histDiv(1:nRec), histDeath(1:nRec), ...
    'VariableNames', {'step','nCells','dEq','minC','maxC', ...
                      'totalDivisions','totalDeaths'});
S.maskSnapshots = maskSnapshots;
S.nCells = nnz(S.occ);
end
