function [Si, STi, info] = sobol(armName, params, cfg0, N)
%SOBOL Phase B2: hand-written Saltelli estimator of Sobol' indices.
%
%   [Si, STi, info] = SOBOL(armName, params, cfg0, N) estimates the
%   first-order (Si) and total-effect (STi) Sobol' indices of the response
%   dist-to-melanoma w.r.t. the knobs `params` (cell array of names), for the
%   winning arm. It uses the Saltelli (2010) two-matrix scheme: (2 + 2*k) * N
%   evaluations. R2023a has no built-in global sensitivity analysis, so the
%   estimator is written by hand.
%
%   info carries the evaluation count, the parameter ranges, and the total
%   response variance, so the report can annotate convergence.

if nargin < 4 || isempty(N), N = cfg0.p3.sobol.N; end

cfg = cfg0;
cfg.sim.N  = cfg0.p3.sobol.Nrun;
cfg.sim.nc = cfg.sim.N / 8;
cfg.sim.stopEqDiameter = round(0.74 * cfg.sim.N);

K = align.knob_levels(armName, cfg);
k = numel(params);
lo = zeros(1, k); hi = zeros(1, k);
for i = 1:k
    v = K.(params{i});
    lo(i) = min(v); hi(i) = max(v);
end

% Two N x k sample matrices (uniform over each knob's range).
s = RandStream('mt19937ar', 'Seed', cfg.measure.stats.seed);
A = rand(s, N, k) .* (hi - lo) + lo;
B = rand(s, N, k) .* (hi - lo) + lo;

% Design: [fA (N); fB (N); fAB_1..k, fBA_1..k (2kN)].
M = 2 * N + 2 * k * N;
X = zeros(M, k);
X(1:N, :) = A;
X(N+1:2*N, :) = B;
off = 2 * N;
for i = 1:k
    ABi = B; ABi(:, i) = A(:, i);
    BAi = A; BAi(:, i) = B(:, i);
    X(off+1:off+N, :) = ABi; off = off + N;
    X(off+1:off+N, :) = BAi; off = off + N;
end

descs = cfg.measure.stats.descriptors;
y = zeros(M, 1);

useParfor = cfg.comparison.useParfor && M > 1;
if useParfor && isempty(gcp('nocreate'))
    try, parpool(cfg.p3.numWorkers); catch, useParfor = false; end
end
if useParfor
    parfor j = 1:M
        y(j) = align.eval_dist(armName, X(j, :), params, cfg, descs);
    end
else
    for j = 1:M
        y(j) = align.eval_dist(armName, X(j, :), params, cfg, descs);
    end
end

fA = y(1:N);
fB = y(N+1:2*N);
varTot = var(fA);

Si = zeros(1, k); STi = zeros(1, k);
off = 2 * N;
for i = 1:k
    fAB = y(off+1:off+N); off = off + N;
    fBA = y(off+1:off+N); off = off + N;
    Si(i)  = mean(fB .* (fAB - fA)) / varTot;   % first-order
    STi(i) = mean((fA - fBA).^2) / (2 * varTot); % total-effect
end

info = struct('params', {params}, 'lo', lo, 'hi', hi, 'nEvals', M, ...
    'N', N, 'varTotal', varTot);

fprintf('Sobol'' indices (dist-to-melanoma), N=%d, %d evals:\n', N, M);
for i = 1:k
    fprintf('  %-14s  Si=%.3f  STi=%.3f\n', params{i}, Si(i), STi(i));
end
end
