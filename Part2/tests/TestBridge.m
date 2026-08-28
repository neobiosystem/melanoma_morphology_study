classdef TestBridge < matlab.unittest.TestCase
%TESTBRIDGE The measurement bridge reuses Part 1's pipeline and reports QC side-effects.
%
%   testDiskRegression is the Part 3 alignment guarantee: a clean synthetic
%   disk passed through the bridge must score exactly as Part 1's own
%   validation suite expects. If it drifts, the bridge (or a copied +morph)
%   has diverged and the alignment is void.

properties
    cfg
end

methods (TestClassSetup)
    function loadConfig(testCase)
        addpaths();
        testCase.cfg = sim_config();
    end
end

methods (Test)

    function testDiskRegression(testCase)
        out = bridge.measure(shapes.disk(400), testCase.cfg);
        testCase.verifyTrue(out.accepted, out.rejectReason);
        testCase.verifyEqual(out.Df,          1.0, 'AbsTol', 0.03, ...
            'Bridge D_f must match Part 1''s clean-disk expectation.');
        testCase.verifyEqual(out.compactness, 1.0, 'AbsTol', 0.03, ...
            'Bridge compactness must match Part 1''s clean-disk expectation.');
        testCase.verifyEqual(out.circularity, 1.0, 'AbsTol', 0.02);
    end

    function testSatelliteIsStrippedAndReported(testCase)
        % A main mass plus a ~2.3% satellite: below the 10% reject gate, above
        % the 1% flag gate, so it is accepted but stripped to the largest
        % component. The bridge must report the stripping rather than hide it.
        mask = false(900, 900);
        mask(100:500, 100:500) = true;      % main mass 401^2
        mask(600:660, 600:660) = true;      % satellite 61^2 (~2.3%)

        out = bridge.measure(mask, testCase.cfg);
        testCase.verifyTrue(out.accepted, out.rejectReason);
        testCase.verifyGreaterThan(out.preStrip.secondCcFraction, 0, ...
            'A second component must be detected before stripping.');
        testCase.verifyGreaterThan(out.preStrip.strippedAreaFraction, 0, ...
            'Stripped area must be reported.');
    end

    function testNecroticHoleIsRecorded(testCase)
        % imfill hides an interior cavity from every boundary descriptor; the
        % bridge must record it so the nutrient mechanism's necrotic core is
        % not silently invisible.
        mask = shapes.disk(300);
        mask(250:260, 250:260) = false;     % interior hole
        out = bridge.measure(mask, testCase.cfg);
        testCase.verifyTrue(out.accepted, out.rejectReason);
        testCase.verifyGreaterThan(out.preStrip.holeAreaFraction, 0, ...
            'A filled hole must be recorded before imfill erases it.');
    end

end
end
