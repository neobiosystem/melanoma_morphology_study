function S = step_nutrient(S)
%STEP_NUTRIENT Update the continuous nutrient field from the current cells.
%
%   S = STEP_NUTRIENT(S) builds the coarse occupancy fraction rho from S.occ,
%   solves the steady-state field with CA.SOLVE_NUTRIENT, and interpolates it
%   back to the fine lattice as S.Cfine for the cell-scale rules.

cfg = S.cfg;

% Occupied-cell fraction per coarse cell. Column-major reshape puts each
% (N/nc) x (N/nc) block contiguous, so mean over dims 1 and 3 is the block mean.
bs = cfg.sim.N / cfg.sim.nc;
occ3 = reshape(S.occ, bs, cfg.sim.nc, bs, cfg.sim.nc);
S.rho = squeeze(mean(mean(occ3, 1), 3));   % nc x nc

S.C     = ca.solve_nutrient(S.rho, cfg);
S.Cfine = imresize(S.C, [cfg.sim.N cfg.sim.N], 'bilinear');
end
