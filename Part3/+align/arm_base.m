function b = arm_base(armName)
%ARM_BASE The benign base parameters for a Part 3 mechanism arm.
%
%   b = ARM_BASE(armName) returns the benign parameter set common to all runs
%   of an arm, including the two new Part 3 knobs (displace / hetero*). Each
%   arm overrides ONE mechanism's parameter(s) from benign toward malignant,
%   so the resulting shape change is attributable to that mechanism alone.
%
%   Lives in a package so it is visible to parfor workers.

b = struct('alpha', 0.4, 'beta', 0.95, 'P_move', 0.02, 'A_i', 0.90, ...
           'T_cell', 4, 'delta', 1e-4, 'nutrientCoupling', false, ...
           'pStarve', 0, 'L', 60, ...
           'displace', false, 'pushRate', 1.0, ...
           'heteroOn', false, 'heteroLambda', 120, 'heteroAmp', 0.5);

switch armName
    case 'mechanical_displacement'
        b.displace = true;              % beta swept by the override
    case 'spatial_heterogeneity'
        b.heteroOn  = true;
        b.heteroAmp = 0.5;              % heteroLambda swept by the override
    otherwise
        error('align:armBase:unknownArm', 'Unknown Part 3 arm: %s', armName);
end
end
