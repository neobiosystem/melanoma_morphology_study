function run_alignment(armName, N, quick)
%RUN_ALIGNMENT Phase B: OAT -> reduced Sobol' -> phase diagram, on the winner.
%
%   run_alignment('spatial_heterogeneity') runs the full Phase B pipeline at
%   the default lattice; quick=true downsizes the Sobol' sample count.
%
%   Outputs: results/oat_<arm>.csv, results/sobol_<arm>.json,
%   results/phase_<arm>.csv.
if nargin < 1 || isempty(armName), armName = 'spatial_heterogeneity'; end
if nargin < 2 || isempty(N), N = 512; end
if nargin < 3, quick = false; end

p3 = fileparts(fileparts(mfilename('fullpath')));
p2 = fullfile(fileparts(p3), 'Part2');
p1 = fullfile(fileparts(p3), 'Part1');
addpath(p3); addpath(fullfile(p3, 'config'));
addpath(p2); addpath(fullfile(p2, 'config'));
addpath(p1); addpath(fullfile(p1, 'config'));

cfg0 = part3_config();
cfg0.sim.N = N; cfg0.sim.nc = N/8; cfg0.sim.stopEqDiameter = round(0.74*N);

% 1. OAT screening.
[T, ranking] = align.oat(armName, cfg0, N);
writetable(T, fullfile(cfg0.paths.results, ['oat_' armName '.csv']));

top = ranking(1:min(2, height(ranking)), :);
params = top.knob;

% 2. Reduced Sobol' on the top-2 knobs.
nSobol = 32;
if ~quick, nSobol = cfg0.p3.sobol.N; end
[Si, STi, info] = align.sobol(armName, params, cfg0, nSobol);
sobolOut = struct('arm', armName, 'params', {params}, 'Si', Si, 'STi', STi, ...
    'lo', info.lo, 'hi', info.hi, 'nEvals', info.nEvals, 'varTotal', info.varTotal);
fid = fopen(fullfile(cfg0.paths.results, ['sobol_' armName '.json']), 'w');
fwrite(fid, jsonencode(sobolOut)); fclose(fid);

% 3. Phase diagram over the top-2 knobs.
P = align.phase_diagram(armName, params{1}, params{2}, cfg0, cfg0.p3.phase.nPts);
writetable(P, fullfile(cfg0.paths.results, ['phase_' armName '.csv']));

fprintf('\nPhase B complete for arm %s:\n', armName);
fprintf('  OAT        -> oat_%s.csv\n', armName);
fprintf('  Sobol''     -> sobol_%s.json (top-2 knobs: %s, %s)\n', armName, params{1}, params{2});
fprintf('  Phase map  -> phase_%s.csv\n', armName);
end
