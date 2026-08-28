classdef TestAlignment < matlab.unittest.TestCase
%TESTALIGNMENT The alignment layer's distance metric and arm plumbing.
%
%   target_distance must reproduce the z-distance formula by hand; arm_base /
%   arm_defs / knob_levels must wire the two arms and their knobs correctly.

    methods (Test)
        function testTargetDistanceZeroAtMelanoma(testCase)
            cfg0 = part3_config();
            targets = jsondecode(fileread(cfg0.p3.targetsJson));
            T = table();
            T.clean_compactness = targets.descriptors.compactness.invasive.median;
            T.clean_solidity   = targets.descriptors.solidity.invasive.median;
            T.clean_convexity  = targets.descriptors.convexity.invasive.median;
            d = align.target_distance(T, cfg0);
            testCase.verifyLessThan(d, 1e-9);
        end

        function testTargetDistanceOneIqrOff(testCase)
            cfg0 = part3_config();
            targets = jsondecode(fileread(cfg0.p3.targetsJson));
            T = table();
            T.clean_compactness = targets.descriptors.compactness.invasive.median;
            T.clean_solidity   = targets.descriptors.solidity.invasive.median;
            T.clean_convexity  = targets.descriptors.convexity.invasive.median;
            iqr = diff(targets.descriptors.compactness.invasive.iqr);
            T.clean_compactness = T.clean_compactness + iqr;  % one IQR off
            d = align.target_distance(T, cfg0);
            testCase.verifyEqual(d, 1, 'AbsTol', 1e-6);
        end

        function testArmBaseFields(testCase)
            bA = align.arm_base('mechanical_displacement');
            testCase.verifyTrue(bA.displace);
            bB = align.arm_base('spatial_heterogeneity');
            testCase.verifyTrue(bB.heteroOn);
            testCase.verifyFalse(bB.displace);
        end

        function testArmDefsAndKnobs(testCase)
            cfg0 = part3_config();
            cfg0.sim.stopEqDiameter = 400;
            arms = align.arm_defs(cfg0);
            testCase.verifyEqual(numel(arms), 2);
            testCase.verifyEqual(numel(arms(2).levels), 5);
            K = align.knob_levels('spatial_heterogeneity', cfg0);
            testCase.verifyTrue(all(isfield(K, {'heteroLambda', 'heteroAmp'})));
        end
    end
end
