classdef TestBoxCounting < matlab.unittest.TestCase
%TESTBOXCOUNTING Validation of the box-counting fractal dimension.
%
%   These tests gate the entire study. D_f is the calibration target that
%   Part 3 fits its simulation against, so an undetected bias here would
%   propagate into the whole closed loop.
%
%   The single most valuable case is testFilledDiskIsDimensionTwo: it feeds
%   the box counter every pixel of a solid disk rather than its outline. A
%   correct implementation must return ~2, and an implementation that has
%   silently been measuring area instead of boundary is caught immediately.

properties
    cfg
end

methods (TestClassSetup)
    function loadConfig(testCase)
        testCase.cfg = default_config();
    end
end

methods (Test)

    function testDiagonalLineIsDimensionOne(testCase)
        % Long enough that the ceil() quantisation of N(eps) = ceil(L/eps)
        % averages out; a short line biases the slope low for purely
        % arithmetic reasons.
        P  = shapes.diagonal_line(2000);
        out = morph.fractal_boxcount(P, 2000, testCase.cfg);
        testCase.verifyEqual(out.Df, 1.0, 'AbsTol', 0.03, ...
            'A straight digital line must have box dimension 1.');
    end

    function testCircleBoundaryIsDimensionOne(testCase)
        mask = shapes.disk(400);
        P    = morph.boundary_pixels(mask);
        dEq  = morph.equivalent_diameter(mask);
        out  = morph.fractal_boxcount(P, dEq, testCase.cfg);
        testCase.verifyEqual(out.Df, 1.0, 'AbsTol', 0.02, ...
            'A smooth closed curve must have box dimension 1.');
    end

    function testSquareBoundaryIsDimensionOne(testCase)
        mask = shapes.square(600);
        P    = morph.boundary_pixels(mask);
        dEq  = morph.equivalent_diameter(mask);
        out  = morph.fractal_boxcount(P, dEq, testCase.cfg);
        % 0.03 rather than 0.02: a square outline covered by eps-boxes needs
        % 4*S/eps - 4 of them, and that corner deficit grows in relative
        % terms as eps rises, biasing the finite-range slope a little above
        % the asymptotic 1.
        testCase.verifyEqual(out.Df, 1.0, 'AbsTol', 0.03, ...
            'A rectilinear closed curve must have box dimension 1.');
    end

    function testFilledDiskIsDimensionTwo(testCase)
        % THE critical test: proves the counter measures the 1-D point set
        % it is given, not the 2-D region that set encloses. If box counting
        % were accidentally applied to filled areas, every lesion D_f in the
        % study would be inflated toward 2 and the whole result set void.
        % Tolerance is 0.08 rather than 0.02 because a digitised disk is not
        % a true 2-D fractal at these scales: N(eps) = A/eps^2 + c*P/eps,
        % and the boundary term biases the measured slope slightly below 2
        % for any finite eps. The point of the test is the decisive gap
        % against the ~1.0 that the disk's OUTLINE returns.
        P   = shapes.filled_disk_points(400);
        out = morph.fractal_boxcount(P, 800, testCase.cfg);
        testCase.verifyEqual(out.Df, 2.0, 'AbsTol', 0.08, ...
            'A filled disk point set must have box dimension 2.');
        testCase.verifyGreaterThan(out.Df, 1.7, ...
            'Filled area must be unambiguously separated from a boundary.');
    end

    function testKochSnowflakeDimension(testCase)
        % The genuine fractal case. Tolerance acknowledges finite recursion
        % depth and rasterisation.
        [P, dEq] = shapes.koch_points(5, 2048);
        out      = morph.fractal_boxcount(P, dEq, testCase.cfg);
        expected = log(4) / log(3);
        testCase.verifyEqual(out.Df, expected, 'AbsTol', 0.05, ...
            sprintf('Koch snowflake must have D_f = %.4f.', expected));
    end

    function testScalingRegionIsReported(testCase)
        mask = shapes.disk(400);
        P    = morph.boundary_pixels(mask);
        dEq  = morph.equivalent_diameter(mask);
        out  = morph.fractal_boxcount(P, dEq, testCase.cfg);

        testCase.verifyGreaterThanOrEqual(out.numPoints, ...
            testCase.cfg.scaling.minPoints);
        testCase.verifyGreaterThanOrEqual(out.spanDecades, ...
            testCase.cfg.scaling.minDecades);
        testCase.verifyGreaterThan(out.R2, 0.99);
        testCase.verifyTrue(out.epsRange(2) > out.epsRange(1));
    end

    function testCountsDecreaseWithBoxSize(testCase)
        % N(eps) must be monotonically non-increasing in eps. A violation
        % would mean the grid-offset minimisation is inconsistent.
        mask = shapes.disk(400);
        P    = morph.boundary_pixels(mask);
        out  = morph.fractal_boxcount(P, morph.equivalent_diameter(mask), testCase.cfg);
        testCase.verifyTrue(all(diff(out.countsAll) <= 0), ...
            'Box counts must not increase as box size grows.');
    end

    function testInsufficientRangeIsFlagged(testCase)
        % A lesion too small to admit a scaling region must be flagged, not
        % silently given a confident slope.
        P   = shapes.diagonal_line(10);
        out = morph.fractal_boxcount(P, 12, testCase.cfg);
        testCase.verifyEqual(out.flag, 'insufficient_range');
        testCase.verifyTrue(isnan(out.Df));
    end

    function testGridOffsetMinimisationLowersCounts(testCase)
        % More origin offsets can only lower the minimal cover. This guards
        % the infimum semantics that make the estimator unbiased.
        mask = shapes.disk(300);
        P    = morph.boundary_pixels(mask);
        dEq  = morph.equivalent_diameter(mask);

        cfgOne = testCase.cfg;  cfgOne.boxcount.numOffsets = 1;
        cfgMany = testCase.cfg; cfgMany.boxcount.numOffsets = 4;

        one  = morph.fractal_boxcount(P, dEq, cfgOne);
        many = morph.fractal_boxcount(P, dEq, cfgMany);

        testCase.verifyTrue(all(many.countsAll <= one.countsAll + 1e-9), ...
            'Minimising over more grid origins must not raise box counts.');
    end

end
end
