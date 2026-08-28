function arms = arm_defs(cfg0)
%ARM_DEFS The two Part 3 mechanism arms and their swept levels.
%
%   arms = ARM_DEFS(cfg0) returns a struct array {name, levels} mirroring Part
%   2's mechanism_comparison arms. Arm A (mechanical_displacement) sweeps beta
%   with push always on; arm B (spatial_heterogeneity) sweeps the correlation
%   length lambda from R/6 to 2R (R = stopEqDiameter/2), so the sweep is
%   scale-invariant across lattice sizes.

R = cfg0.sim.stopEqDiameter / 2;

arms = struct();

arms(1).name   = 'mechanical_displacement';
arms(1).levels = { struct('beta', 0.95), struct('beta', 0.70), ...
                   struct('beta', 0.40), struct('beta', 0.20), ...
                   struct('beta', 0.05) };

arms(2).name   = 'spatial_heterogeneity';
arms(2).levels = { struct('heteroLambda', R/6), struct('heteroLambda', R/4), ...
                   struct('heteroLambda', R/2), struct('heteroLambda', R), ...
                   struct('heteroLambda', 2*R) };
end
