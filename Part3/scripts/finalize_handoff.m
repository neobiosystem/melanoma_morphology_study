function finalize_handoff()
%FINALIZE_HANDOFF Augment part3_handoff.json with the Phase B results.
%
%   After OAT, reduced Sobol' and the phase diagram have run, this script
%   re-reads the Phase A handoff and adds a phaseB block: the winning
%   mechanism's aligned parameter region (from the phase-diagram minimum), the
%   residual alignment error, and the Sobol' indices. The result is the
%   complete Part 3 -> Part 4 handoff.

p3 = fileparts(fileparts(mfilename('fullpath')));
addpath(p3); addpath(fullfile(p3, 'config'));
res = fullfile(p3, 'results');
handoffPath = fullfile(res, 'part3_handoff.json');

H = jsondecode(fileread(handoffPath));

phasePath = fullfile(res, 'phase_spatial_heterogeneity.csv');
if isfile(phasePath)
    P = readtable(phasePath, 'TextType', 'string');
    [dmin, bi] = min(P.dist);
    r = P(bi, :);
    H.phaseB = struct();
    H.phaseB.winner = 'spatial_heterogeneity';
    H.phaseB.bestParams = struct('heteroAmp', r.heteroAmp, 'heteroLambda', r.heteroLambda);
    H.phaseB.minDist = dmin;
    H.phaseB.bestClean = struct('compactness', r.clean_compactness, ...
        'solidity', r.clean_solidity, 'convexity', r.clean_convexity);
end

sobolPath = fullfile(res, 'sobol_spatial_heterogeneity.json');
if isfile(sobolPath)
    SB = jsondecode(fileread(sobolPath));
    H.phaseB.sobol = struct('params', {SB.params}, 'Si', SB.Si, 'STi', SB.STi, ...
        'nEvals', SB.nEvals, 'varTotal', SB.varTotal);
end

fid = fopen(handoffPath, 'w');
fwrite(fid, jsonencode(H));
fclose(fid);
fprintf('Updated %s\n', handoffPath);
end
