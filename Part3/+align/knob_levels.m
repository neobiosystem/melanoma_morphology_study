function K = knob_levels(armName, cfg0)
%KNOB_LEVELS Per-knob sweep levels for a winning arm (Phase B).
%
%   K = KNOB_LEVELS(armName, cfg0) returns a struct whose fields are the arm's
%   free knobs, each a 1 x 5 vector from benign-ish to malignant-ish. Used by
%   OAT (one-at-a-time) and by Sobol'/phase-diagram (ranges from these levels).

R = cfg0.sim.stopEqDiameter / 2;
switch armName
    case 'mechanical_displacement'
        K.beta     = [0.95 0.70 0.40 0.20 0.05];
        K.pushRate = [0.2 0.4 0.6 0.8 1.0];
    case 'spatial_heterogeneity'
        K.heteroLambda = [R/6 R/4 R/2 R 2*R];
        K.heteroAmp    = [0.1 0.3 0.5 0.7 1.0];
    otherwise
        error('align:knobLevels:unknownArm', 'Unknown arm: %s', armName);
end
end
