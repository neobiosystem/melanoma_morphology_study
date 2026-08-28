classdef TestGrowthField < matlab.unittest.TestCase
%TESTGROWTHFIELD The growth-rate field is smooth, centred, and deterministic.
%
%   These tests pin the properties the spatial-heterogeneity arm relies on:
%   a mean-1 field, smooth on the cell scale (no single-cell white noise), a
%   patch size that grows with lambda, and seed determinism.

    methods (Test)
        function testMeanNearOne(testCase)
            g = mech.growth_field(256, 32, 0.5, 1);
            testCase.verifyEqual(mean(g(:)), 1, 'AbsTol', 0.05);
        end

        function testSmoothness(testCase)
            % A smooth field (lambda >> 1 cell) has small adjacent-cell
            % gradients; white noise would give O(gAmp) single-cell jumps.
            g = mech.growth_field(256, 32, 0.5, 1);
            dx = abs(diff(g, 1, 2));
            dy = abs(diff(g, 1, 1));
            testCase.verifyLessThan(max([dx(:); dy(:)]), 0.5);
        end

        function testDeterministic(testCase)
            g1 = mech.growth_field(256, 32, 0.5, 7);
            g2 = mech.growth_field(256, 32, 0.5, 7);
            testCase.verifyEqual(g1, g2);
        end

        function testSeedChangesField(testCase)
            g1 = mech.growth_field(256, 32, 0.5, 7);
            g2 = mech.growth_field(256, 32, 0.5, 8);
            testCase.verifyNotEqual(g1, g2);
        end

        function testLargerLambdaFewerPatches(testCase)
            % Larger correlation length -> fewer contiguous above-mean patches.
            gSmall = mech.growth_field(256, 16, 0.5, 3);
            gLarge = mech.growth_field(256, 64, 0.5, 3);
            nSmall = bwconncomp(gSmall > 1, 8).NumObjects;
            nLarge = bwconncomp(gLarge > 1, 8).NumObjects;
            testCase.verifyLessThan(nLarge, nSmall);
        end
    end
end
