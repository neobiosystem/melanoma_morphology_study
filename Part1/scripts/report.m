%% Part 1 — Morphological descriptors of melanocytic lesion borders
% *Fractal dimension and compactness from ISIC 2018 expert and automatic
% segmentations, stratified by nevus / melanoma in situ / invasive melanoma.*
%
% This report is generated from |results/descriptors.csv|. Publish it with
%
%   publish('scripts/report.m', struct('format','html','outputDir','results/report'))
%
% MATLAB R2023a cannot open plain-text Live Code (that needs R2025a+), so this
% is a |%%|-sectioned script published to HTML rather than an |.mlx|.
%
% The numbers here are not an end in themselves. They are the calibration
% targets for Part 3's cell-dynamics simulation, so the report leads with what
% would make them *wrong* — rasterisation bias, annotation style, resolution
% heterogeneity — before it reports what they are.

%% Setup
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'scripts'));   % run_tests, analyze_results live here
addpath(fullfile(root, 'tests'));
cfg = default_config();

R = readtable(fullfile(cfg.paths.results, 'descriptors.csv'), 'TextType', 'string');
R.grade = categorical(R.grade, cfg.data.gradeOrder, 'Ordinal', true);
labels  = cfg.data.gradeLabels;

% readtable reads the accepted flags back as numeric 0/1, which cannot be
% used as a logical index. Normalise once rather than sprinkling logical()
% at every use site, and blank the <missing> strings for the same reason.
for v = ["expert_accepted", "otsu_accepted"]
    if ismember(v, R.Properties.VariableNames)
        R.(v) = logical(R.(v));
    end
end
for v = ["expert_flags", "expert_DfFlag", "expert_rejectReason", ...
         "otsu_flags", "otsu_DfFlag", "otsu_rejectReason"]
    if ismember(v, R.Properties.VariableNames)
        R.(v)(ismissing(R.(v))) = "";
    end
end

fprintf('Cohort as measured: %d lesions\n', height(R));

%% 1. Cohort composition
% ISIC 2018 Task1-2 Training, collection 63. A single source was chosen so
% that every image and its mask come from the same acquisition and annotation
% pipeline; pooling cohorts would confound grade with annotation protocol.
%
% Grade is *not* in the challenge archive — the zips carry masks only — so it
% is joined back from the ISIC Archive API by |isic_id| via the current
% hierarchical taxonomy field |diagnosis_3| (the legacy |benign_malignant| and
% |mel_class| fields were retired in 2024).

disp(groupsummary(R, 'grade'));

fprintf('\nHistopathology-confirmed: %d / %d\n', nnz(R.histoConfirmed), height(R));

%% 2. What the measurement pipeline does
% Both segmentation routes funnel through the *same* |morph.process_lesion|,
% differing only in where the mask comes from. Any route disagreement is
% therefore attributable to segmentation, not to two subtly different
% measurement implementations.
%
% Three design choices carry most of the validity burden:
%
% * *Box dimension is an infimum over covers*, so N(eps) is *minimised* over a
%   3x3 lattice of grid origins, never averaged. Averaging estimates a
%   different quantity.
% * *D_f is computed on the raw boundary only.* Smoothing before box-counting
%   would erase precisely the signal being measured. Only perimeter-derived
%   descriptors see the smoothed contour.
% * *eps is made dimensionless* (eps/D_eq) instead of resampling masks to a
%   common size. Resampling a 4500 px lesion down to 512 px would destroy fine
%   boundary detail preferentially on the highest-resolution images — that is,
%   it would homogenise the very thing under study.

disp(cfg.boxcount);
disp(cfg.scaling);

%% 3. Validation against shapes with known dimension
% The decisive test is the *filled disk*: fed a solid mask, a correct
% implementation returns D_f ~ 2.0 because it is covering a 2-D set. An
% implementation that silently measures the boundary of the filled region
% would return ~1.0 and every downstream number would be meaningless.

% Via run_tests rather than a bare runtests: the suite needs root, config and
% tests all on the path, and calling runtests directly silently errors every
% test instead of running them.
results = run_tests();
disp(table({results.Name}', [results.Passed]', [results.Duration]', ...
    'VariableNames', {'test','passed','seconds'}));

%% 4. QC outcomes
% Gates are reported per-reason rather than as a single pass rate, because
% *which* lesions are lost matters more than how many. Border contact is
% judged as a fraction of boundary length rather than any-contact-at-all:
% an any-contact rule discards large lesions preferentially, and lesion size
% tracks severity, so that loss would land differentially across grades and
% manufacture a between-grade difference out of nothing.

fprintf('expert route accepted: %d / %d\n', nnz(R.expert_accepted), height(R));
rr = R.expert_rejectReason(R.expert_rejectReason ~= "");
if ~isempty(rr), disp(groupcounts(table(rr), 'rr')); end

fprintf('\nbox-count scaling-region flags (expert route):\n');
disp(groupcounts(table(R.expert_DfFlag(R.expert_accepted), ...
    'VariableNames', {'flag'}), 'flag'));

%% 5. Descriptor intervals by grade — the headline deliverable
% Median [Q1, Q3] is the primary form because the distributions are
% non-normal and the in-situ group is small. Where n >= 20 a 2.5-97.5
% percentile reference interval is also given; that is the form Part 3
% consumes.

% Primary analysis is the automatic route (see section 8 for why), with the
% expert route reported alongside as the sensitivity arm.
fprintf('########## PRIMARY: automatic route ##########\n');
A = analyze_results('otsu');

fprintf('\n\n########## SENSITIVITY: expert-mask route ##########\n');
Aexpert = analyze_results('expert'); %#ok<NASGU>

%% 6. Perimeter rasterisation bias — the largest false-positive risk
% Raster perimeter is systematically inflated, and inflated *more* for
% irregular boundaries. Left uncorrected it would make melanomas appear less
% compact for purely digital reasons, which is exactly the study's headline
% claim arriving as an artefact.
%
% The check: if a between-grade compactness gap survives fixed-physical-scale
% smoothing it is geometric; if it collapses once sigma exceeds a few pixels,
% it was rasterisation.

% Run on the PRIMARY (automatic) route: that is where the compactness
% difference actually lives, so that is where the rasterisation artefact
% would have to be ruled out. Aacc keeps the expert subset for section 8.
Aacc = R(R.expert_accepted, :);
Aprim = R(R.otsu_accepted, :);
g = categories(Aprim.grade);
rawMed = arrayfun(@(k) median(Aprim.otsu_compactnessRaw(Aprim.grade == g{k}), 'omitnan'), 1:numel(g));
smMed  = arrayfun(@(k) median(Aprim.otsu_compactness(Aprim.grade == g{k}),    'omitnan'), 1:numel(g));

figure; bar([rawMed(:), smMed(:)]);
set(gca, 'XTickLabel', labels);
legend({'raw raster perimeter','fixed-scale smoothed'}, 'Location','best');
ylabel('median compactness'); grid on;
title('Compactness before and after perimeter-bias mitigation');

fprintf('raw    gap benign->invasive: %+.4f\n', rawMed(end) - rawMed(1));
fprintf('smooth gap benign->invasive: %+.4f\n', smMed(end)  - smMed(1));

%% 7. Dual-route agreement
% Agreement is stratified by grade because the stratification is itself a
% result. If automatic-versus-expert agreement falls monotonically with
% severity, that is a biological statement — malignant borders are objectively
% harder to delimit. If agreement is flat across grades, the automatic route
% is a valid substitute on datasets without expert masks. Both outcomes are
% findings, not appendix material.

both = R(logical(R.expert_accepted) & logical(R.otsu_accepted), :);
if ~isempty(both) && any(isfinite(both.dice))
    disp(groupsummary(both, 'grade', {'median','mean'}, {'dice','jaccard','areaRatio'}));

    figure;
    boxplot(both.dice, both.grade, 'Labels', labels, 'Symbol', '.');
    ylabel('Dice (expert vs automatic)'); grid on;
    title('Segmentation agreement by grade');

    fprintf('\n%-18s %6s %8s %9s %9s %12s\n', ...
        'descriptor','n','CCC','bias','sdDiff','reliability');
    fprintf('%s\n', repmat('-', 1, 68));
    for d = 1:numel(cfg.stats.descriptors)
        name = cfg.stats.descriptors{d};
        ec = ['expert_' name];  oc = ['otsu_' name];
        if ~all(ismember({ec, oc}, both.Properties.VariableNames)), continue; end
        ag = agree.descriptor_agreement(both.(ec), both.(oc), both.grade);
        fprintf('%-18s %6d %8.3f %9.4f %9.4f %12.2f\n', ...
            name, ag.n, ag.ccc, ag.bias, ag.sdDiff, ag.reliabilityRatio);
    end

    % Df stratified, since that is the descriptor whose route choice is
    % actually under debate.
    agDf = agree.descriptor_agreement(both.expert_Df, both.otsu_Df, both.grade);
    fprintf('\nD_f agreement by grade:\n');
    disp(agDf.byGroup);
end

%% 8. Why the automatic route is primary
% ISIC 2018 ground-truth masks are substantially polygon-traced rather than
% pixel-accurate: median collinear-boundary-vertex fraction 0.587, with 24.7%
% of masks above 0.85. A polygon boundary is locally straight by construction,
% so its box dimension is pinned near 1 regardless of the lesion underneath.
%
% The effect turns out to be worse than attenuation. For benign versus
% invasive on expert masks, 5 of 7 descriptors run *backwards* relative to
% biological expectation -- invasive melanoma measures more compact, more
% solid and more convex than benign nevi -- and all 5 reversals survive Holm
% correction (p_adj 8.6e-5 to 1.4e-2). This is stated pairwise on purpose: a
% direction claim is a pairwise claim, and the omnibus Kruskal-Wallis p in
% route_comparison.csv is a three-group statistic that cannot support it.
%
% The reversals are significant but SMALL: rank-biserial 0.13 to 0.21,
% eta^2 0.003 to 0.009. At n = 2088 a small effect reaches significance
% easily, so the defensible claim is that the bias is systematic and
% reproducible, not that it is large. Symmetrically, of the 2 descriptors
% pointing the right way only asymmetry is significant; circularity is not
% (p_adj = 0.22). On the automatic route all 7 point the expected way, all 7
% are significant, and the effect sizes are an order of magnitude larger.
%
% A sign reversal cannot be explained by "less detail": missing detail pushes
% an effect toward zero, it does not push it past zero. The most plausible
% mechanism is that annotation coarseness scales with lesion size -- larger,
% more irregular lesions receive proportionally coarser polygons, making their
% masks *smoother* than small benign ones.
%
% That mechanism makes a falsifiable prediction, and the data agree with it:
% the 2 descriptors that do NOT reverse are exactly the 2 that do not depend
% on boundary detail. Asymmetry is a principal-axis fold overlap of the mask
% and circularity is the dispersion of centroid-to-boundary radii; both
% measure where the mass sits and are insensitive to local boundary wiggle.
% All 5 reversed descriptors depend on boundary detail -- perimeter length
% (compactness, convexity, perimeter decay), concavity (solidity), or the
% fine boundary itself (D_f).
%
% The automatic route is therefore primary for ALL shape descriptors, not just
% D_f, with the expert route retained as the sensitivity arm. Two competing
% explanations for the automatic route's advantage were tested rather than
% assumed (section 8b); both survive.
%
% D_f remains the WEAKEST of the seven descriptors (eta^2 0.024 vs 0.087 for
% compactness). The headline morphological signal in this cohort is carried by
% the perimeter- and hull-derived measures, not by fractal dimension.

% Show the pairwise evidence rather than asserting it. export_posthoc writes
% one row per (route, descriptor, group pair) straight from the saved
% analyses, so the directions and p values below cannot drift from them.
export_posthoc();
PH = readtable(fullfile(cfg.paths.results, 'posthoc_all.csv'), 'TextType', 'string');
BI = PH(PH.groupA == "benign" & PH.groupB == "invasive", :);

expectedUp = struct('Df', 1, 'compactness', 1, 'circularity', -1, ...
    'asymmetry', 1, 'solidity', -1, 'convexity', -1, 'perimeterDecaySlope', -1);

fprintf('\nBenign vs invasive, by route (direction vs biological expectation):\n');
fprintf('%-8s %-20s %10s %10s %8s %9s  %s\n', ...
    'route', 'descriptor', 'medBenign', 'medInvasive', 'p_adj', 'rankBis', 'verdict');
fprintf('%s\n', repmat('-', 1, 88));
for r = 1:height(BI)
    name = char(BI.descriptor(r));
    obs  = sign(BI.medianB(r) - BI.medianA(r));
    good = obs == expectedUp.(name);
    sig  = BI.pAdjusted(r) < cfg.stats.alpha;
    if good
        verdict = 'expected';
    else
        verdict = 'REVERSED';
    end
    if ~sig
        verdict = [verdict ' (n.s.)']; %#ok<AGROW>
    end
    fprintf('%-8s %-20s %10.4f %10.4f %8.2g %9.3f  %s\n', ...
        BI.route(r), name, BI.medianA(r), BI.medianB(r), ...
        BI.pAdjusted(r), BI.rankBiserial(r), verdict);
end

fprintf('\nD_f by grade, expert route:\n');
disp(groupsummary(Aacc, 'grade', {'median','mean','std'}, 'expert_Df'));

Botsu = R(R.otsu_accepted, :);
if ~isempty(Botsu)
    fprintf('\nD_f by grade, automatic route:\n');
    disp(groupsummary(Botsu, 'grade', {'median','mean','std'}, 'otsu_Df'));
end

%% 8b. The two confounds, tested
% *Segmentation quality.* If the automatic route simply segments malignant
% lesions worse, its grade separation would be a quality gradient wearing the
% costume of biology. Agreement does decline with severity (Jonckheere-
% Terpstra z = -3.4, p = 0.0007), so the concern is real. But restricted to
% well-segmented lesions only (Dice >= 0.75, n = 1232) all seven descriptors
% remain significant with only modest attenuation. Under the artefact
% explanation the effect sizes would collapse; they do not.
%
% *Residual hair.* Morphological opening is disabled because it eroded the
% lesion periphery, which leaves hair suppression to DullRazor alone. Hair is
% thin, dark and boundary-attached -- exactly what inflates a box dimension --
% and D_f does rise with hair fraction (Spearman rho = 0.168, p = 1.6e-11),
% with hair burden itself differing by grade. So this was tested too:
% restricted to the least-hairy 40%, D_f's effect size nearly DOUBLES
% (eta^2 0.024 -> 0.046). Hair adds noise that masks the signal rather than
% manufacturing it.
%
% Run |route_confound_check| to regenerate both tests.

route_confound_check();

%% 9. Calibration targets for Part 3
% Written to |results/calibration_targets_otsu.json| (primary) and
% |..._expert.json| (sensitivity) so Part 3 loads the intervals
% directly. Hand-transcribing numbers between stages of a multi-part study is
% where quiet errors take root.

targetPath = fullfile(cfg.paths.results, 'calibration_targets_otsu.json');
if isfile(targetPath)
    disp(fileread(targetPath));
else
    fprintf('Run analyze_results(''otsu'') to generate %s\n', targetPath);
end

%% 10. Licensing
% ISIC images are licensed per-image (a CC-0 / CC-BY / CC-BY-NC mix). The
% per-image license and attribution string are retained alongside the results
% so any figure can be traced to its permissions at submission time.

disp(groupcounts(table(R.license, 'VariableNames', {'license'}), 'license'));
