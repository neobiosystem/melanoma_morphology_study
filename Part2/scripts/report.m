%% Part 2 — Cell-dynamics simulation core
% *A hybrid cellular automaton plus reaction-diffusion nutrient field, scored
% by the same morphological pipeline as Part 1, to test which cell-scale
% mechanism reproduces the measured malignant signature.*
%
% Publish with
%
%   publish('scripts/report.m', struct('format','html','outputDir','results/report'))
%
% The scientific question here is not "can we grow a tumour?" but "which of
% four mechanisms -- loss of contact inhibition, nutrient-driven fingering,
% enhanced migration, or clonal heterogeneity -- moves the shape from the
% benign to the melanoma region of descriptor space, and can any of them?"
% Part 1 measured that the malignant signature is LOW solidity, LOW convexity
% and HIGH compactness at near-constant D_f: a few large lobes, not fine
% fractal fringe.

%% Setup
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'scripts'));
addpath(fullfile(root, 'tests'));
p1 = fullfile(fileparts(root), 'Part1');
addpath(p1);
addpath(fullfile(p1, 'config'));
addpath(fullfile(p1, 'tests'));

cfg = sim_config();

%% 1. The measured target (from Part 1, not re-typed)
% The calibration target is loaded directly from Part 1's JSON so the numbers
% cannot drift. The signature the simulation must explain:

targets = jsondecode(fileread(fullfile(p1, 'results', 'calibration_targets_otsu.json')));
prim = {'compactness', 'solidity', 'convexity'};
for d = prim
    nm = d{1};
    ben = targets.descriptors.(nm).benign.median;
    inv = targets.descriptors.(nm).invasive.median;
    fprintf('%-12s  benign %.3f  ->  invasive %.3f  (%s)\n', ...
        nm, ben, inv, direction_word(nm, ben, inv));
end
fprintf('%-12s  benign %.4f  ->  invasive %.4f  (near-constant: D_f is secondary)\n', ...
    'Df', targets.descriptors.Df.benign.median, targets.descriptors.Df.invasive.median);

%% 2. What the simulator does
% A hybrid CA (Anderson 2005): discrete cells on an N x N lattice coupled to a
% continuous reaction-diffusion nutrient field solved to steady state on a
% coarse grid (sparse direct solve, not explicit relaxation). Rules:
%
% * *Proliferation* gated by maturity (T_cell, fixed), empty neighbour,
%   contact inhibition beta, and local nutrient.
% * *Migration* gated by adhesion A_i, with an optional chemotactic bias chi.
% * *Apoptosis*: basal Poisson deaths plus starvation necrosis.
%
% *Lattice isotropy* is the named risk in the design document. Moore
% connectivity with diagonal hops down-weighted by 1/sqrt(2) restores
% isotropy; TestIsotropy asserts the residual 4-fold Fourier harmonic is
% small and the resulting cluster is circular.
%
% *Identifiability:* alpha and T_cell both set the effective division rate, so
% T_cell is fixed here; beta gates division while A_i gates migration so the
% two "neighbour-count" knobs do not collapse into one.

disp(phen.provenance_table());

%% 3. Validation suite
% Mirrors Part 1: synthetic geometry and analytic checks before any sweep is
% trusted. The two load-bearing tests are TestIsotropy (the lattice must not
% imprint an artefact) and TestBridge (the measurement bridge must still agree
% with Part 1's own clean-disk expectation).

results = run_tests();
passed = nnz([results.Passed]);
failed = nnz([results.Failed]);
fprintf('\n%d/%d tests passed.\n', passed, numel(results));
if failed > 0
    error('Validation suite failed; fix before interpreting any sweep.');
end

%% 4. Mechanism comparison
% Four arms, each holding one mechanism fixed at benign and sweeping its
% parameter toward malignant, every other parameter fixed at the benign base.
% Each mask is measured clean and degraded; the primary reading is the
% between-level DELTA, because the clean CA baseline does not sit on the real
% benign point (Eden roughness vs real large-scale irregularity plus
% segmentation noise). NOTE the nutrient arm sweeps L from short (20, strong
% fingering, the aggressive end) to long (600, uniform, the benign end), so its
% "malignant" end is level 1 -- opposite the other three arms.

csvPath = fullfile(cfg.paths.results, 'mechanism_trajectories.csv');
if ~isfile(csvPath)
    fprintf('Run scripts/mechanism_comparison.m first to generate %s\n', csvPath);
else
    T = readtable(csvPath, 'TextType', 'string');
    armNames = categories(categorical(T.arm));

    fprintf('\nPer-arm trajectory (median over seeds; clean descriptors):\n');
    for a = 1:numel(armNames)
        nm = armNames{a};
        sub = T(T.arm == nm, :);
        nRej = sum(sub.accepted == 0);
        fprintf('\n  %-20s  %d runs (%d rejected by QC)\n', nm, height(sub), nRej);
        fprintf('    %-6s %10s %10s %10s %10s\n', 'level', 'compact.', 'solidity', 'convexity', 'Df');
        for lv = unique(sub.level)'
            r = sub(sub.level == lv, :);
            if all(r.accepted == 0)
                fprintf('    %-6d %10s %10s %10s %10s  <-- rejected (necrotic core)\n', lv, '-', '-', '-', '-');
                continue;
            end
            r = r(r.accepted == 1, :);
            fprintf('    %-6d %10.3f %10.3f %10.3f %10.4f\n', lv, ...
                median(r.clean_compactness), median(r.clean_solidity), ...
                median(r.clean_convexity), median(r.clean_Df));
        end
    end

    fprintf('\nStripped-area fraction (migration phenotype the QC gate hides):\n');
    for a = 1:numel(armNames)
        nm = armNames{a};
        sub = T(T.arm == nm, :);
        fprintf('  %-20s  max %.5f\n', nm, max(sub.strippedAreaFraction));
    end
end

%% 5. Handoff to Part 3
% Best per arm (closest in the 3 primary descriptors to the real melanoma
% median) plus the honest distance, written to part2_handoff.json.

handoffPath = fullfile(cfg.paths.results, 'part2_handoff.json');
if isfile(handoffPath)
    disp(jsondecode(fileread(handoffPath)));
else
    fprintf('Run scripts/mechanism_comparison.m first to generate %s\n', handoffPath);
end

%% 6. Conclusions and limits
% Four mechanisms were swept at production scale (1024^2, 5 levels x 3 seeds)
% and scored by the SAME pipeline that scored the real lesions. The result is
% a clean negative -- and the plan pre-registered that a negative result is
% itself a finding, not a licence to tune parameters until something fits.
%
% * *Contact-inhibition loss (beta -> 0)* moves solidity UP (0.978 -> 0.990)
%   and convexity UP (0.994 -> 0.999): the OPPOSITE of malignant. Structurally
%   an empty-neighbour CA cannot express "push through neighbours" -- a
%   surrounded cell has no empty site to divide into -- so lowering beta only
%   lets every frontier cell divide uniformly, which heals gaps and smooths
%   the margin.
% * *Clonal heterogeneity (advantaged subclone, frac -> 0.4)* also smooths the
%   cluster (compactness 1.025 -> 1.013, solidity 0.978 -> 0.983). A faster
%   clone fills the interior and rounds the mass; it does not displace
%   neighbours into lobes.
% * *Nutrient-driven fingering* moves in the CORRECT direction (shorter
%   diffusion length L -> more irregular: compactness 1.021 -> 1.033, solidity
%   0.979 -> 0.976), but the amplitude is ~30x too small, and the strongest
%   endpoint L=20 is rejected outright (truncated_at_border): the necrotic
%   core (hole fraction 1.91) fragments the mask into ~8500 components and
%   imfill turns it into a blob touching the frame border. L>=150 is a uniform
%   field -- no gradient, no fingering.
% * *Enhanced migration* changes no descriptor (compactness 1.024 -> 1.017,
%   solidity 0.977 -> 0.980). Its phenotype -- satellite detachment -- raises
%   the stripped-area fraction ~13x but stays two orders of magnitude below
%   the 1% strip threshold, so the QC gate erases it before scoring.
%
% Conclusion: no arm's trajectory crosses the real benign and melanoma
% regions. The measured malignant signature (solidity 0.85, convexity 0.84,
% compactness 1.81 at near-constant Df) is "a few large lobes", whereas an
% empty-neighbour lattice CA -- however its four knobs are swept -- yields
% near-circular Eden clusters, not large-scale lobulation. Reproducing the
% signature will require a mechanism that displaces neighbours (off-lattice /
% mechanical) or imposes long-range spatial heterogeneity; that is the
% concrete, negative input handed to Part 3, not a parameter Part 2 could have
% found.
%
% Limits carried into the handoff: (i) clean CA masks carry no segmentation
% noise, so absolute values differ from the real Otsu medians and the primary
% reading is the between-level delta; (ii) migration satellites and necrotic
% cores are erased by QC, so the bridge reports them separately; (iii) lattice
% isotropy is asserted by TestIsotropy, not assumed.

% -------------------------------------------------------------------------
function w = direction_word(nm, ben, inv)
if inv > ben, w = 'increases'; else, w = 'decreases'; end
end
