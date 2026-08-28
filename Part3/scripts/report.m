%% Part 3 — Mechanism-to-morphology alignment
% *Two candidate mechanisms -- local mechanical displacement and long-range
% spatial heterogeneity -- are added to the Part 2 cellular automaton and
% scored by the same morphological pipeline, to find which one reproduces the
% measured malignant signature (few large lobes).*
%
% Publish with
%
%   publish('scripts/report.m', struct('format','html','outputDir','results/report'))
%
% Part 2 ended in a clean negative: an empty-neighbour lattice CA, however its
% four knobs are swept, yields near-circular Eden clusters (distance to the
% melanoma median ~1.15), because a surrounded cell has no empty site to divide
% into. The handoff named two concrete gaps -- *displace neighbours* (off-lattice
% / mechanical) or *impose long-range spatial heterogeneity* -- and this part
% implements both, compares them, and aligns the winner.

%% Setup
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'scripts'));
addpath(fullfile(root, 'tests'));
p2 = fullfile(fileparts(root), 'Part2');
addpath(p2); addpath(fullfile(p2, 'config'));
p1 = fullfile(fileparts(root), 'Part1');
addpath(p1); addpath(fullfile(p1, 'config'));

cfg = part3_config();
res = cfg.paths.results;

%% 1. The measured target (from Part 1, not re-typed)
% The alignment objective is the z-like distance across the three primary
% descriptors -- compactness, solidity, convexity -- to the real invasive
% medians. D_f is secondary: Part 1 showed it barely moves between grades.

targets = jsondecode(fileread(cfg.p3.targetsJson));
prim = cfg.p3.primary;
fprintf('Alignment objective (lower dist = closer to melanoma):\n');
for d = prim
    nm = d{1};
    ben = targets.descriptors.(nm).benign.median;
    inv = targets.descriptors.(nm).invasive.median;
    fprintf('  %-12s  benign %.3f  ->  invasive %.3f  (%s)\n', nm, ben, inv, ...
        direction_word(nm, ben, inv));
end

%% 2. What the two new mechanisms do
% Both extend Part 2's CA at the proliferation rule (ca3.step_proliferate):
%
% * *Mechanical displacement* -- a surrounded (layer-2) cell that passes the
%   division test no longer gives up: it pushes one occupied neighbour into an
%   adjacent empty site and divides into the vacated slot (O(1) local push). This
%   lifts the "no empty neighbour" bottleneck so the interior can grow too.
% * *Spatial heterogeneity* -- a smooth static growth-rate field g(x,y) (sum of
%   low-order Fourier modes, mean 1, correlation length lambda) multiplies each
%   cell's division rate. A few broad growth-rate patches leave a few large lobes.
%
% Both knobs are pre-registered; the data, not the author, picks the winner.

%% 3. Validation suite
% Mirrors Part 2: conservation, determinism, growth-field properties, and -- the
% load-bearing check -- that with both mechanisms off the lattice is byte-identical
% to Part 2's CA, so any anisotropy from the new knobs is signal, not a grid
% artefact.

run_tests();   % errors internally if any test fails

%% 4. Phase A — mechanism discovery
% Two arms x 5 levels x 3 seeds at 512^2, scored clean and degraded. The
% success criterion is whether either arm drives (compactness, solidity,
% convexity) toward the melanoma medians and drops dist below Part 2's ~1.15.

csvPath = fullfile(res, 'discovery_trajectories.csv');
T = readtable(csvPath, 'TextType', 'string');
armNames = categories(categorical(T.arm));
fprintf('\nPer-arm best (closest run to melanoma median):\n');
for a = 1:numel(armNames)
    nm = armNames{a};
    sub = T(T.arm == nm & T.accepted == 1, :);
    d = align.target_distance(sub, cfg);
    [dmin, bi] = min(d);
    r = sub(bi, :);
    fprintf('  %-24s dist=%.3f  (compact %.3f, solidity %.3f, convexity %.3f)\n', ...
        nm, dmin, r.clean_compactness, r.clean_solidity, r.clean_convexity);
end

h = jsondecode(fileread(fullfile(res, 'part3_handoff.json')));
fprintf('\nWinner: %s (dist %.3f, vs Part 2 baseline ~1.15)\n', h.winner, h.winnerDistance);

%% 5. Phase B — OAT screening
% One-at-a-time screening of the winning arm's two knobs (growth-field
% correlation length lambda and amplitude gAmp), 5 levels x 3 seeds each.

oatPath = fullfile(res, 'oat_spatial_heterogeneity.csv');
if isfile(oatPath)
    O = readtable(oatPath, 'TextType', 'string');
    base = mean(O.dist(strcmp(O.knob, 'baseline')), 'omitnan');
    fprintf('OAT baseline dist = %.3f\n', base);
    knobs = {'heteroAmp', 'heteroLambda'};
    for ki = 1:numel(knobs)
        idx = strcmp(O.knob, knobs{ki});
        lvls = unique(O.level(idx));
        mx = 0; bestLv = 0;
        for lv = lvls'
            d = mean(O.dist(idx & O.level == lv), 'omitnan');
            if abs(d - base) > mx, mx = abs(d - base); bestLv = lv; end
        end
        fprintf('  %-14s max |Delta dist| = %.3f (level %d)\n', knobs{ki}, mx, bestLv);
    end
end

%% 6. Phase B — reduced Sobol' sensitivity
% Hand-written Saltelli estimator over the two dominant knobs, response = dist.

sobolPath = fullfile(res, 'sobol_spatial_heterogeneity.json');
if isfile(sobolPath)
    SB = jsondecode(fileread(sobolPath));
    fprintf('Sobol'' indices (N=%d, %d evals):\n', SB.N, SB.nEvals);
    for i = 1:numel(SB.params)
        fprintf('  %-14s  S_i=%.3f  S_Ti=%.3f\n', SB.params{i}, SB.Si(i), SB.STi(i));
    end
    fprintf('  total response variance = %.4f\n', SB.varTotal);
end

%% 7. Phase B — phase diagram
% The dist surface over (gAmp, lambda), one seed per cell. The annotated
% minimum is the parameter region that best reproduces the malignant signature.

phasePath = fullfile(res, 'phase_spatial_heterogeneity.csv');
if isfile(phasePath)
    P = readtable(phasePath, 'TextType', 'string');
    [dmin, bi] = min(P.dist);
    r = P(bi, :);
    fprintf('Phase-diagram minimum: dist=%.3f at gAmp=%.3f, lambda=%.1f px\n', ...
        dmin, r.heteroAmp, r.heteroLambda);
    fprintf('  -> compact %.3f, solidity %.3f, convexity %.3f (targets %.3f/%.3f/%.3f)\n', ...
        r.clean_compactness, r.clean_solidity, r.clean_convexity, ...
        targets.descriptors.compactness.invasive.median, ...
        targets.descriptors.solidity.invasive.median, ...
        targets.descriptors.convexity.invasive.median);
end

%% 8. Handoff to Part 4
% The winning mechanism and its aligned parameter region, plus the residual
% alignment error, written to part3_handoff.json.

if isfile(fullfile(res, 'part3_handoff.json'))
    disp(jsondecode(fileread(fullfile(res, 'part3_handoff.json'))));
end

%% 9. Conclusions and limits
% The growth-field mechanism reproduces the malignant signature that Part 2
% could not: a few large lobes (low solidity, low convexity, high compactness)
% at near-constant D_f, landing within dist ~0.3-0.4 of the melanoma median --
% roughly a third of Part 2's ~1.15 baseline. Local mechanical displacement
% does not: it still converges to a near-circle (dist ~1.13), confirming the
% pre-registered expectation that an incompressible pushed droplet relaxes to
% its minimal-energy round shape.
%
% Limits carried into the handoff: (i) the growth field is an imposed smooth
% heterogeneity, not an emergent instability -- Part 4 must ask whether the
% correlated patches can arise from cell-autonomous dynamics; (ii) clean CA
% masks carry no segmentation noise, so absolute descriptor values differ from
% the real Otsu medians; (iii) the alignment is judged on the three primary
% descriptors only, with D_f reported secondarily.

% -------------------------------------------------------------------------
function w = direction_word(nm, ben, inv)
if inv > ben, w = 'increases'; else, w = 'decreases'; end
end
