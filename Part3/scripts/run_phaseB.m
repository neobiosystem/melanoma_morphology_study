function run_phaseB(armName, N, nSobol)
%RUN_PHASEB Phase B (phase diagram + reduced Sobol') on the winning arm.
%
%   run_phaseB('spatial_heterogeneity') builds the 2D phase diagram over the
%   top-2 OAT knobs (heteroAmp x heteroLambda) and the reduced Saltelli Sobol'
%   indices, both at lattice N. OAT must already be written to
%   results/oat_<arm>.csv (run align.oat first) -- the top-2 knobs are read
%   from its ranking rather than re-swept.
%
%   Outputs: results/phase_<arm>.csv, results/sobol_<arm>.json.
if nargin < 1 || isempty(armName), armName = 'spatial_heterogeneity'; end
if nargin < 2 || isempty(N), N = 512; end

p3 = fileparts(fileparts(mfilename('fullpath')));
p2 = fullfile(fileparts(p3), 'Part2');
p1 = fullfile(fileparts(p3), 'Part1');
addpath(p3); addpath(fullfile(p3, 'config'));
addpath(p2); addpath(fullfile(p2, 'config'));
addpath(p1); addpath(fullfile(p1, 'config'));

cfg0 = part3_config();
cfg0.sim.N = N; cfg0.sim.nc = N/8; cfg0.sim.stopEqDiameter = round(0.74*N);

% Top-2 knobs from the OAT ranking (deterministic, same as run_alignment).
R = cfg0.sim.stopEqDiameter / 2;
K = align.knob_levels(armName, cfg0);
knobs = fieldnames(K);

% Re-derive the OAT ranking cheaply from the saved table so this script is
% self-contained and the knob order never drifts.
T = readtable(fullfile(cfg0.paths.results, ['oat_' armName '.csv']));
base = mean(T.dist(strcmp(T.knob, 'baseline')), 'omitnan');
delta = zeros(numel(knobs), 1);
for ki = 1:numel(knobs)
    nm = knobs{ki};
    idx = strcmp(T.knob, nm);
    if ~any(idx), continue; end
    lvls = unique(T.level(idx));
    mx = 0;
    for lv = lvls'
        d = mean(T.dist(idx & T.level == lv), 'omitnan');
        mx = max(mx, abs(d - base));
    end
    delta(ki) = mx;
end
[~, ord] = sort(delta, 'descend');
params = knobs(ord(1:min(2, numel(knobs))));

fprintf('Phase B: top-2 knobs = %s, %s\n', params{1}, params{2});

% 1. Phase diagram (one seed per cell -- a phase map, not a Monte Carlo mean).
P = align.phase_diagram(armName, params{1}, params{2}, cfg0, cfg0.p3.phase.nPts);
writetable(P, fullfile(cfg0.paths.results, ['phase_' armName '.csv']));
fprintf('Wrote phase_%s.csv (%d rows)\n', armName, height(P));

% 2. Reduced Saltelli Sobol' on the same top-2 knobs.
if nargin < 3 || isempty(nSobol), nSobol = cfg0.p3.sobol.N; end
[Si, STi, info] = align.sobol(armName, params, cfg0, nSobol);
sobolOut = struct('arm', armName, 'params', {params}, 'Si', Si, 'STi', STi, ...
    'lo', info.lo, 'hi', info.hi, 'nEvals', info.nEvals, 'varTotal', info.varTotal);
fid = fopen(fullfile(cfg0.paths.results, ['sobol_' armName '.json']), 'w');
fwrite(fid, jsonencode(sobolOut)); fclose(fid);
fprintf('Wrote sobol_%s.json (%d evals)\n', armName, info.nEvals);
end
