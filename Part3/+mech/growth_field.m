function g = growth_field(N, lambda, gAmp, seed)
%GROWTH_FIELD Smooth static growth-rate field with correlation length lambda.
%
%   g = GROWTH_FIELD(N, lambda, gAmp, seed) returns an N x N field whose mean
%   is 1 and whose amplitude is gAmp. It is built by summing a small number of
%   cosine waves with wavenumber ~ pi/lambda and random phase/direction, which
%   yields a smooth "quenched heterogeneity": a handful of contiguous
%   fast/slow patches of size ~lambda, i.e. a few large lobes rather than fine
%   noise. A LOCAL RandStream is used so the global CA stream (rng(seed) in
%   run_simulation) is left untouched and run determinism is preserved.

s = RandStream('mt19937ar', 'Seed', seed);
[XX, YY] = meshgrid(1:N, 1:N);
k0 = pi / max(lambda, 8);          % characteristic wavenumber (patch size ~lambda)
nModes = 12;

g = zeros(N, N);
for m = 1:nModes
    th = 2 * pi * rand(s);                      % wave direction
    kx = k0 * cos(th) * (0.5 + rand(s));        % wavenumber jitter 0.5..1.5 k0
    ky = k0 * sin(th) * (0.5 + rand(s));
    ph = 2 * pi * rand(s);
    a  = rand(s);
    g = g + a * cos(kx * (XX - 1) + ky * (YY - 1) + ph);
end

g = g - mean(g(:));
g = 1 + gAmp * g / std(g(:));
g = min(max(g, 0.1), 3.0);
end
