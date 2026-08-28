function b = arm_base(armName)
%ARM_BASE The benign base parameters for a mechanism-comparison arm.
%
%   b = ARM_BASE(armName) returns the benign parameter set common to all runs
%   of an arm. Each arm then overrides ONE mechanism's parameter(s) from
%   benign toward malignant, so the resulting shape change is attributable to
%   that mechanism alone.
%
%   This lives in a package (not a script-local function) so it is visible to
%   parfor workers.

b = struct('alpha', 0.4, 'beta', 0.95, 'P_move', 0.02, 'A_i', 0.90, ...
           'T_cell', 4, 'delta', 1e-4, 'nutrientCoupling', false, ...
           'pStarve', 0, 'L', 60, 'subclone_frac', 0, 'subclone_alphaMult', 3);

switch armName
    case 'nutrient_fingering'
        b.alpha = 0.6; b.nutrientCoupling = true; b.pStarve = 0.01;
    case 'enhanced_migration'
        % benign base; P_move / A_i swept by the override
    case 'clonal_heterogeneity'
        % benign base; subclone_frac swept by the override
    case 'contact_inhibition'
        % benign base; beta swept by the override
end
end
