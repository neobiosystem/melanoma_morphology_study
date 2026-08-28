function d = eval_dist(armName, x, params, cfg, descs)
%EVAL_DIST Distance-to-melanoma for one parameter vector.
%
%   d = EVAL_DIST(armName, x, params, cfg, descs) maps a single parameter
%   vector x (1 x numel(params)) onto its alignment distance by running one
%   simulation. It is the scalar response evaluated by the Sobol' sampler.

ov = struct();
for i = 1:numel(params)
    ov.(params{i}) = x(i);
end
r = struct('arm', armName, 'level', 0, 'seed', 1, 'overrides', ov);
row = align.run_one(r, cfg, descs);
d = align.target_distance(row, cfg);
% A QC-rejected mask (e.g. the tumour overgrew the frame at large lambda) has
% no measurable descriptors, so dist is NaN. Map it to a large finite value:
% "not melanoma-like" is the honest reading, and it keeps the Sobol' variance
% estimate well-defined instead of silently poisoning it.
if ~isfinite(d), d = 3.0; end
end
