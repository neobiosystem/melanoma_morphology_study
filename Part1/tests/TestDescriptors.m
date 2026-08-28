classdef TestDescriptors < matlab.unittest.TestCase
%TESTDESCRIPTORS Validation of the shape descriptors and QC gates.
%
%   Two groups of tests matter most here:
%
%   * testSmoothingRemovesPerimeterBias proves the staircase mitigation
%     works. Raster perimeter is inflated, and inflated more for irregular
%     boundaries, so without this correction melanomas would appear less
%     compact for purely digital reasons.
%
%   * testAsymmetryAndElongationAreIndependent proves the descriptor set is
%     not secretly redundant, the same failure the reciprocal-circularity
%     definition would have introduced.

properties
    cfg
end

methods (TestClassSetup)
    function loadConfig(testCase)
        testCase.cfg = default_config();
    end
end

methods (Test)

    % ---- Compactness ---------------------------------------------------

    function testCompactnessOfCircleIsOne(testCase)
        C = shapes.circle_contour(400);
        testCase.verifyEqual(morph.compactness(C), 1.0, 'RelTol', 0.001, ...
            'A circle is the minimum-compactness shape, P^2/(4 pi A) = 1.');
    end

    function testCompactnessOfSquareIsFourOverPi(testCase)
        C = shapes.square_contour(600);
        testCase.verifyEqual(morph.compactness(C), 4/pi, 'RelTol', 0.01, ...
            'A square has compactness 16 L^2 / (4 pi L^2) = 4/pi.');
    end

    function testCompactnessOfRasterCircleAfterSmoothing(testCase)
        mask = shapes.disk(400);
        dEq  = morph.equivalent_diameter(mask);
        C    = morph.boundary_trace(mask);
        Cs   = morph.smooth_contour(C, round(testCase.cfg.smooth.primaryRatio * dEq));
        testCase.verifyEqual(morph.compactness(Cs), 1.0, 'AbsTol', 0.02, ...
            'A rasterised circle must recover compactness ~1 once smoothed.');
    end

    % ---- Perimeter bias ------------------------------------------------

    function testSmoothingRemovesPerimeterBias(testCase)
        radius   = 400;
        truePerim = 2 * pi * radius;

        mask = shapes.disk(radius);
        dEq  = morph.equivalent_diameter(mask);
        C    = morph.boundary_trace(mask);
        Cs   = morph.smooth_contour(C, round(testCase.cfg.smooth.primaryRatio * dEq));

        rawErr    = abs(morph.contour_perimeter(C)  - truePerim) / truePerim;
        smoothErr = abs(morph.contour_perimeter(Cs) - truePerim) / truePerim;

        testCase.verifyGreaterThan(rawErr, 0.02, ...
            'Raw raster perimeter should be measurably inflated.');
        testCase.verifyLessThan(smoothErr, 0.02, ...
            'Fixed-scale smoothing must bring perimeter within 2% of truth.');
        testCase.verifyLessThan(smoothErr, rawErr, ...
            'Smoothing must reduce, not increase, perimeter error.');
    end

    function testPerimeterDecreasesWithSmoothing(testCase)
        mask = shapes.disk(400);
        dEq  = morph.equivalent_diameter(mask);
        ms   = morph.multi_scale(mask, dEq, testCase.cfg);
        testCase.verifyTrue(all(diff(ms.perimeters) <= 1e-6), ...
            'Perimeter must be non-increasing as the smoothing scale grows.');
    end

    % ---- Circularity ---------------------------------------------------

    function testCircularityOfCircleIsOne(testCase)
        P = shapes.circle_points(400);
        testCase.verifyEqual(morph.circularity(P, [0 0]), 1.0, 'AbsTol', 0.001, ...
            'All radii equal means zero radial deviation, circularity 1.');
    end

    function testCircularityOfEllipseIsBelowCircle(testCase)
        mask   = shapes.ellipse(400, 200);
        P      = morph.boundary_pixels(mask);
        props  = regionprops(mask, 'Centroid');
        cRC    = [props(1).Centroid(2), props(1).Centroid(1)];
        testCase.verifyLessThan(morph.circularity(P, cRC), 0.95, ...
            'A 2:1 ellipse is radially anisotropic, so circularity drops.');
    end

    function testCircularityIsNotReciprocalCompactness(testCase)
        % Guards the descriptor-redundancy fix: an ellipse must separate the
        % two measures. If circularity were 4 pi A / P^2 this would fail.
        mask = shapes.ellipse(400, 200);
        dEq  = morph.equivalent_diameter(mask);
        P    = morph.boundary_pixels(mask);
        props = regionprops(mask, 'Centroid');
        cRC  = [props(1).Centroid(2), props(1).Centroid(1)];

        C    = morph.boundary_trace(mask);
        Cs   = morph.smooth_contour(C, round(testCase.cfg.smooth.primaryRatio * dEq));

        circ = morph.circularity(P, cRC);
        comp = morph.compactness(Cs);

        testCase.verifyGreaterThan(abs(circ - 1/comp), 0.05, ...
            'Circularity must carry information independent of compactness.');
    end

    % ---- Asymmetry -----------------------------------------------------

    function testAsymmetryOfDiskIsZero(testCase)
        asym = morph.asymmetry_index(shapes.disk(300));
        testCase.verifyLessThan(asym.mean, 0.02);
        testCase.verifyLessThan(asym.central, 0.02);
    end

    function testAsymmetryAndElongationAreIndependent(testCase)
        % A symmetric ellipse is elongated but not asymmetric. If these two
        % descriptors were entangled this separation would collapse.
        asym = morph.asymmetry_index(shapes.ellipse(300, 150));
        testCase.verifyLessThan(asym.mean, 0.03, ...
            'A symmetric ellipse must have near-zero asymmetry.');
        testCase.verifyGreaterThan(asym.elongation, 0.25, ...
            'A 2:1 ellipse must have clearly non-zero elongation.');
    end

    function testAsymmetryOfHalfDiskIsPositive(testCase)
        % A half disk is strongly asymmetric about its MAJOR axis (the flat
        % cut) yet perfectly symmetric about its minor axis. The two-fold
        % mean therefore understates it, so the discriminating assertion is
        % on max. Asserting the minor fold is ~0 pins down that the axis
        % alignment is working rather than accidentally passing.
        mask = shapes.disk(300);
        mask(1:round(size(mask,1)/2), :) = false;    % lop off half
        asym = morph.asymmetry_index(mask);

        testCase.verifyGreaterThan(asym.max, 0.15, ...
            'A half disk must be strongly asymmetric about one axis.');
        testCase.verifyLessThan(asym.minor, 0.05, ...
            'A half disk remains symmetric about its minor axis.');
    end

    % ---- Multi-scale ---------------------------------------------------

    function testSolidityAndConvexityOfDiskAreOne(testCase)
        mask = shapes.disk(400);
        ms   = morph.multi_scale(mask, morph.equivalent_diameter(mask), testCase.cfg);
        testCase.verifyEqual(ms.solidity,  1.0, 'AbsTol', 0.02);
        testCase.verifyEqual(ms.convexity, 1.0, 'AbsTol', 0.03);
    end

    function testSolidityAndConvexityNeverExceedOne(testCase)
        % Both are ratios against a convex hull, so both are bounded above by
        % 1 by definition: no set has more area than its hull, and no closed
        % curve is shorter than its hull's boundary. Values above 1 are not
        % noise, they mean numerator and denominator were measured in
        % different representations -- which is exactly the bug this pins
        % down (pixel-count area vs contour-polygon area, and raw-contour
        % hull perimeter vs smoothed contour perimeter).
        shapesToCheck = {shapes.disk(400), shapes.ellipse(400, 200), ...
                         shapes.square(500)};
        for k = 1:numel(shapesToCheck)
            mask = shapesToCheck{k};
            ms = morph.multi_scale(mask, morph.equivalent_diameter(mask), testCase.cfg);
            testCase.verifyLessThanOrEqual(ms.solidity, 1 + 1e-9, ...
                sprintf('shape %d: solidity must not exceed 1.', k));
            testCase.verifyLessThanOrEqual(ms.convexity, 1 + 1e-9, ...
                sprintf('shape %d: convexity must not exceed 1.', k));
        end
    end

    function testNotchedShapeHasLowSolidity(testCase)
        mask = shapes.disk(300);
        [h, w] = size(mask);
        mask(round(h/2):end, round(w/2)-30:round(w/2)+30) = false;   % deep notch
        ms = morph.multi_scale(mask, morph.equivalent_diameter(mask), testCase.cfg);
        testCase.verifyLessThan(ms.solidity, 0.95, ...
            'A notched shape must lose solidity against its convex hull.');
    end

    function testCircleHasFlatPerimeterDecay(testCase)
        mask = shapes.disk(400);
        ms   = morph.multi_scale(mask, morph.equivalent_diameter(mask), testCase.cfg);
        testCase.verifyLessThan(abs(ms.perimeterDecaySlope), 0.05, ...
            'A smooth circle sheds almost no perimeter under smoothing.');
    end

    % ---- QC gates ------------------------------------------------------

    function testEmptyMaskIsRejected(testCase)
        out = morph.process_lesion(false(100), testCase.cfg);
        testCase.verifyFalse(out.accepted);
        testCase.verifyEqual(out.rejectReason, 'empty_mask');
    end

    function testBorderTouchingMaskIsRejected(testCase)
        mask = false(400);
        mask(1:200, 1:200) = true;     % runs into the top-left corner
        out = morph.process_lesion(mask, testCase.cfg);
        testCase.verifyFalse(out.accepted);
        testCase.verifyEqual(out.rejectReason, 'truncated_at_border');
    end

    function testLesionGrazingBorderIsKeptAndFlagged(testCase)
        % A lesion that barely reaches the frame edge keeps valid geometry.
        % Rejecting on any contact at all would preferentially discard large
        % lesions, and lesion size tracks severity, so that loss would land
        % differentially across grades.
        % A square with a narrow stalk reaching the bottom edge: contact is
        % ~21 px against a boundary of ~2300, well under the 2% gate.
        mask = false(700, 700);
        mask(100:600, 100:600) = true;
        mask(550:700, 300:320) = true;
        out  = morph.process_lesion(mask, testCase.cfg);

        testCase.verifyTrue(out.accepted, out.rejectReason);
        testCase.verifyTrue(ismember('grazes_border', out.flags), ...
            'Minor border contact must still be recorded as a flag.');
        testCase.verifyGreaterThan(out.borderContactFraction, 0);
        testCase.verifyLessThanOrEqual(out.borderContactFraction, ...
            testCase.cfg.qc.maxBorderContactFrac);
    end

    function testTinyLesionIsRejected(testCase)
        out = morph.process_lesion(shapes.disk(10), testCase.cfg);
        testCase.verifyFalse(out.accepted);
        testCase.verifyEqual(out.rejectReason, 'lesion_too_small');
    end

    function testTwoComponentsAreRejected(testCase)
        mask = false(900, 900);
        mask(100:400, 100:400) = true;
        mask(500:800, 500:800) = true;    % equal-sized foreign object
        out = morph.process_lesion(mask, testCase.cfg);
        testCase.verifyFalse(out.accepted);
        testCase.verifyEqual(out.rejectReason, 'multiple_components');
    end

    function testCleanDiskIsAcceptedWithSaneDescriptors(testCase)
        out = morph.process_lesion(shapes.disk(400), testCase.cfg);
        testCase.verifyTrue(out.accepted, out.rejectReason);
        testCase.verifyEqual(out.Df,          1.0, 'AbsTol', 0.03);
        testCase.verifyEqual(out.compactness, 1.0, 'AbsTol', 0.03);
        testCase.verifyEqual(out.circularity, 1.0, 'AbsTol', 0.02);
        testCase.verifyLessThan(out.asymmetry, 0.02);
    end

end
end
