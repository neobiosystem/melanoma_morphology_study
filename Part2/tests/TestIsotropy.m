classdef TestIsotropy < matlab.unittest.TestCase
%TESTISOTROPY The headline lattice-anisotropy guard.
%
%   A pure-proliferation cluster with no nutrient field and no migration must
%   grow into a circle, not a diamond or square. A square lattice with a
%   careless neighbour rule grows a spurious 4-fold shape, which would show up
%   as a fake compactness/solidity difference between arms -- i.e. an artefact
%   wearing the costume of a mechanism. This test is a hard gate: if it fails,
%   every morphological conclusion is suspect.

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

    function testProliferationIsIsotropic(testCase)
        cfg = testCase.cfg;
        cfg.sim.N = 256; cfg.sim.nc = 64;
        cfg.sim.nutrientCoupling = false;   % no field -> no drift source
        cfg.sim.T_cell = 1; cfg.sim.alpha = 1.0; cfg.sim.beta = 0.0;
        cfg.sim.P_move = 0; cfg.sim.A_i = 0.9;
        cfg.sim.delta = 0; cfg.sim.pStarve = 0;
        cfg.sim.seed = 7;
        cfg.sim.snapshotInterval = 10; cfg.sim.maskSaveInterval = 500;
        cfg.sim.maxSteps = 400; cfg.sim.stopEqDiameter = 120;
        cfg.sim.initRadius = 5;

        S = ca.run_simulation(cfg);
        mask = S.occ;

        % Centroid in image coordinates [x y].
        props = regionprops(mask, 'Centroid');
        cx = props(1).Centroid(1);
        cy = props(1).Centroid(2);

        P = morph.boundary_pixels(mask);
        ang = atan2(P(:, 2) - cy, P(:, 1) - cx);   % image-coordinate angle
        r   = hypot(P(:, 2) - cy, P(:, 1) - cx);

        % Radial profile in 36 angular bins, then the 4-fold Fourier harmonic.
        nb = 36;
        edges = linspace(-pi, pi, nb + 1);
        [~, ~, bin] = histcounts(ang, edges);
        rMean = accumarray(bin, r, [nb 1], @mean, NaN);
        rMean(isnan(rMean)) = mean(r, 'omitnan');

        f = abs(fft(rMean));
        ratio4 = f(5) / f(1);        % harmonic 4 = index 5 (index 1 is DC)

        testCase.verifyLessThan(ratio4, 0.08, ...
            sprintf(['4-fold Fourier harmonic %.3f too large: the lattice is ' ...
            'imprinting a square/diamond artefact.'], ratio4));

        % Smoothed compactness must stay near 1 (circle), not a polygon.
        out = bridge.measure(mask, cfg);
        testCase.verifyTrue(out.accepted, out.rejectReason);
        testCase.verifyLessThan(out.compactness, 1.15, ...
            sprintf('compactness %.3f: cluster is not close to circular.', out.compactness));
    end

end
end
