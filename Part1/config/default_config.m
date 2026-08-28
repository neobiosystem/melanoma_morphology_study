function cfg = default_config()
%DEFAULT_CONFIG Single source of truth for every path, threshold and hyperparameter.
%
%   cfg = DEFAULT_CONFIG() returns a struct consumed by every function in
%   Part 1. Nothing in +morph, +seg, +stats or +agree may hard-code a
%   numeric threshold: if a reviewer could question a number, it lives here.
%
%   A copy of this struct is written into results/ on every run so that the
%   outputs are self-describing. Part 3 cites these descriptors as
%   calibration targets, so the provenance of each constant matters.
%
%   See also RUN_ANALYSIS, MORPH.PROCESS_LESION.

% ---- Paths -------------------------------------------------------------
% Resolved relative to this file so the pipeline runs from any working dir.
thisDir     = fileparts(mfilename('fullpath'));
cfg.root    = fileparts(thisDir);

cfg.paths.raw          = fullfile(cfg.root, 'data', 'raw');
cfg.paths.images       = fullfile(cfg.paths.raw, 'ISIC2018_Task1-2_Training_Input');
cfg.paths.masks        = fullfile(cfg.paths.raw, 'ISIC2018_Task1_Training_GroundTruth');
cfg.paths.manifest     = fullfile(cfg.root, 'data', 'manifest.csv');
cfg.paths.metadata     = fullfile(cfg.root, 'data', 'isic_metadata.csv');
cfg.paths.results      = fullfile(cfg.root, 'results');
cfg.paths.figures      = fullfile(cfg.root, 'results', 'figures');
cfg.paths.qcFigures    = fullfile(cfg.root, 'results', 'qc');
% Curated exemplars live at the top level rather than under results/: they are
% a hand-checkable deliverable for figures and talks, not a regenerated output.
cfg.paths.classicExamples = fullfile(cfg.root, 'classical_examples');

% ---- Data source -------------------------------------------------------
% ISIC 2018 Task 1-2 Training == collection 63 in the ISIC Archive API.
% The challenge archives ship masks but NO diagnosis, so the three-tier
% grade must be joined back from Archive metadata by isic_id.
cfg.data.apiBase          = 'https://api.isic-archive.com/api/v2/';
cfg.data.collectionId     = 63;
cfg.data.pageLimit        = 100;      % API page size for cursor pagination
cfg.data.imageZipUrl      = 'https://isic-challenge-data.s3.amazonaws.com/2018/ISIC2018_Task1-2_Training_Input.zip';
cfg.data.maskZipUrl       = 'https://isic-challenge-data.s3.amazonaws.com/2018/ISIC2018_Task1_Training_GroundTruth.zip';
cfg.data.maskSuffix       = '_segmentation.png';

% Three-tier grading. ISIC moved to a hierarchical taxonomy in 2024; the
% legacy benign_malignant / mel_class fields no longer exist and
% diagnosis_3 now carries this split.
cfg.data.gradeField  = 'diagnosis_3';
cfg.data.gradeMap    = { ...
    'Nevus',              'benign'; ...
    'Melanoma in situ',   'in_situ'; ...
    'Melanoma Invasive',  'invasive'};
cfg.data.gradeOrder  = {'benign', 'in_situ', 'invasive'};   % ordinal severity
cfg.data.gradeLabels = {'Nevus', 'Melanoma in situ', 'Invasive melanoma'};

% ---- Box counting ------------------------------------------------------
% Scale bounds are expressed as ratios of the equivalent diameter D_eq so
% that the fitting window is identical across a 10x resolution spread.
% Masks are NEVER resampled to a common size: downsampling a 4500 px lesion
% would erase exactly the fine boundary detail under study, and would do so
% preferentially on the highest-resolution images.
cfg.boxcount.numScales     = 25;    % K, log-spaced
cfg.boxcount.numOffsets    = 3;     % q; q^2 = 9 grid origins per scale
cfg.boxcount.minBoxPx      = 4;     % raster floor: below this a digital
                                    % boundary looks like a straight line
                                    % and D_f collapses toward 1
% 0.005 rather than a coarser ratio: with dEq typically 800-2000 px, a 0.02
% ratio would put the lower bound at 16-40 px and discard the fine-scale
% border detail that distinguishes infiltrative from smooth margins, while
% also squeezing the fit to ~1.1 decades. The 4 px floor above is what
% actually protects against raster noise.
cfg.boxcount.minBoxRatio   = 0.005; % lower bound as a fraction of D_eq
cfg.boxcount.maxBoxRatio   = 0.25;  % above this N(eps) saturates

% ---- Scaling-region detection -----------------------------------------
% Local-slope plateau, NOT bare R^2 maximisation: maximising R^2 alone
% latches onto tiny windows where few points make R^2 trivially high.
cfg.scaling.slopeWindow    = 5;     % points per local-slope OLS window
% 0.10 is the calibrated middle. At 0.05 no plateau is ever found and every
% lesion silently falls back to a full-range fit, reintroducing the
% end-truncation bias the scaling region exists to remove. At 0.15 the
% window is loose enough to swallow genuine curvature: a digitised solid
% disk, whose local slope really does drift 1.98 -> 1.79, gets accepted as
% one plateau and returns 1.92 instead of ~1.97.
cfg.scaling.slopeTol       = 0.10;  % plateau tolerance about run median
cfg.scaling.minDecades     = 0.6;   % minimum span of the fitted region
cfg.scaling.minPoints      = 7;     % minimum points in the fitted region
cfg.scaling.crossCheckTol  = 0.03;  % plateau vs R^2 method disagreement

% ---- Contour smoothing -------------------------------------------------
% Raster perimeter is inflated by staircase effects, and inflated MORE for
% irregular boundaries, which would manufacture a fake between-grade
% difference in compactness. Smoothing windows are a fixed fraction of
% D_eq (not a fixed pixel count) so the correction is resolution-independent.
cfg.smooth.sigmaRatios     = [0, 1/80, 1/40, 1/20, 1/10, 1/5];
cfg.smooth.primaryRatio    = 1/20;  % window used for reported compactness

% ---- Quality-control gates --------------------------------------------
cfg.qc.minEqDiameterPx     = 50;    % reject: too small to resolve boundary
% Backstop only, deliberately set where it cannot fire for a real image.
% ISIC dermoscopy runs up to 6748 px wide, so a lesion with dEq > 4096 px is
% perfectly ordinary there; maxFrameFraction below is the gate that actually
% detects an implausible mask, expressed relative to the frame rather than
% in absolute pixels.
cfg.qc.maxEqDiameterPx     = 8192;
cfg.qc.computeCapPx        = 2048;  % above this, downsample for compute only
cfg.qc.maxFrameFraction    = 0.80;  % reject: lesion fills frame (likely crop)

% Border contact. Rejecting every lesion that grazes the frame discards
% ~12% of the cohort and biases it toward small lesions, since big lesions
% run out of frame more often -- and lesion size correlates with severity,
% so that bias would land differentially across grades. Only substantial
% truncation actually invalidates perimeter and area, so a lesion is
% rejected on the FRACTION of its boundary lying on the frame edge, and
% merely flagged below that.
cfg.qc.maxBorderContactFrac = 0.02;
cfg.qc.rejectSecondCcFrac  = 0.10;  % reject: foreign object present
cfg.qc.flagSecondCcFrac    = 0.01;  % flag: minor speckle, keep largest
cfg.qc.minBoxcountR2       = 0.99;  % flag below this
cfg.qc.maxSlopeStdErr      = 0.03;  % flag above this
cfg.qc.minRouteDice        = 0.50;  % flag Otsu-vs-expert for manual review
cfg.qc.rejectBorderTouch   = true;  % truncated lesion -> P and A invalid

% ---- Segmentation (Otsu reproducibility arm) --------------------------
cfg.seg.hairSeLengthRatio  = 1/40;  % DullRazor line SE length, fraction of image min side
cfg.seg.hairNumOrient      = 6;     % orientations for the closing bank
cfg.seg.hairThreshold      = 10;    % intensity delta marking a hair pixel
cfg.seg.channel            = 'Lab_L';   % 'Lab_L' | 'HSV_V' | 'HSV_H'
% Morphological opening is DISABLED (0). It was intended to strip hairs, but
% measured against expert masks it mostly eroded the lesion periphery: with
% it on, automatic masks recovered ~64% of expert area; with it off, ~82%.
% Hair suppression is left to SEG.REMOVE_HAIR alone.
cfg.seg.openRadiusRatio    = 0;
cfg.seg.closeRadiusRatio   = 1/100;

% Active-contour refinement is DISABLED. It was meant to recover detail lost
% to the opening, but on this cohort it shrank masks slightly and did not
% improve agreement once the two fixes below were in place.
cfg.seg.useActiveContour   = false;
cfg.seg.activeContourIter  = 60;

% Compute the Otsu threshold only over the dermatoscope aperture, so dark
% vignette corners cannot bias it toward the lesion core. This was the single
% largest improvement to MEAN Dice (0.554 -> 0.672), i.e. it rescues
% catastrophic failures rather than nudging typical cases.
cfg.seg.thresholdInsideAperture = true;

% These three settings were chosen on a 90-lesion stratified sample by
% minimising the SPREAD of area ratio across grades, not by maximising Dice.
% Uniform under-segmentation shifts all grades together and largely cancels
% in a between-grade comparison; grade-dependent under-segmentation deforms
% the comparison itself, and would do so in the direction of the hypothesis.
% Spread fell from 0.171 (baseline) to 0.066 with these settings, while Dice
% improved as well, so no tradeoff was actually required.
% See scripts/tune_otsu.m and results/otsu_tuning.csv.

% Vignette removal guard. imclearborder is only allowed to win if what
% survives is still a plausible lesion; below this fraction of the largest
% pre-clearing component we assume the lesion itself touched the frame and
% keep it. Without this, a border-spanning lesion plus one interior speck
% yields the speck.
cfg.seg.minInteriorRetainRatio = 0.30;

% Dice above which a pair of routes is treated as "both segmented this
% lesion sensibly". Used to test whether the automatic route's between-grade
% separation survives when poorly segmented lesions are excluded -- i.e.
% whether that separation is geometry or segmentation failure.
cfg.qc.wellSegmentedDice = 0.75;

% Quantile of hair fraction defining the "least hairy" subset. Used to test
% whether a between-grade D_f difference survives where there is almost no
% hair to confound it, which separates "hair inflates D_f" from "hair causes
% the grade difference in D_f".
cfg.qc.lowHairQuantile = 0.40;

% ---- Statistics --------------------------------------------------------
cfg.stats.alpha            = 0.05;
cfg.stats.numBootstrap     = 2000;
cfg.stats.seed             = 42;
cfg.stats.minNForStats     = 10;    % below this, report n only
cfg.stats.minNForRefInt    = 20;    % below this, no 2.5-97.5 interval
cfg.stats.refIntervalPct   = [2.5, 97.5];
cfg.stats.postHocCorrect   = 'holm';    % within-descriptor pairwise
cfg.stats.omnibusCorrect   = 'bh';      % across-descriptor FDR
cfg.stats.numBalancedReps  = 100;   % balanced-subsample sensitivity reps

% Descriptors carried through the statistical pipeline.
cfg.stats.descriptors = { ...
    'Df', 'compactness', 'circularity', 'asymmetry', ...
    'solidity', 'convexity', 'perimeterDecaySlope'};

% ---- Execution ---------------------------------------------------------
cfg.run.useParfor          = true;
cfg.run.numWorkers         = 8;
cfg.run.smokeTestCount     = 20;    % lesions in the pre-flight subset
cfg.run.saveQcFigures      = true;
cfg.run.qcFigureCount      = 30;    % per-lesion montages to render

end
