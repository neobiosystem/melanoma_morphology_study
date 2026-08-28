classdef TestCA < matlab.unittest.TestCase
%TESTCA Correctness of the CA update rules, independent of morphology.
%
%   These tests do not care what the tumour looks like; they pin down the
%   invariants any cellular automaton must satisfy: cells are conserved, the
%   cluster grows, and a fixed seed reproduces bit-for-bit.

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

    function testCellCountConservation(testCase)
        cfg = testCase.cfg;
        cfg.sim.N = 128; cfg.sim.nc = 64; cfg.sim.maxSteps = 60;
        cfg.sim.nutrientCoupling = false;
        cfg.sim.T_cell = 2; cfg.sim.alpha = 0.5; cfg.sim.beta = 0.9;
        cfg.sim.P_move = 0.1; cfg.sim.A_i = 0.5;
        cfg.sim.delta = 1e-3; cfg.sim.pStarve = 0;
        cfg.sim.seed = 11; cfg.sim.snapshotInterval = 20;
        cfg.sim.maskSaveInterval = 100; cfg.sim.stopEqDiameter = 60;

        S = ca.run_simulation(cfg);
        expected = S.initCells + S.totalDivisions - S.totalDeaths;
        testCase.verifyEqual(S.nCells, expected, ...
            'Cells must satisfy n = n0 + divisions - deaths.');
    end

    function testGrowthIsMonotonic(testCase)
        cfg = testCase.cfg;
        cfg.sim.N = 128; cfg.sim.nc = 64; cfg.sim.maxSteps = 80;
        cfg.sim.nutrientCoupling = false;
        cfg.sim.T_cell = 2; cfg.sim.alpha = 0.7; cfg.sim.beta = 0.9;
        cfg.sim.P_move = 0; cfg.sim.A_i = 0.9;
        cfg.sim.delta = 0; cfg.sim.pStarve = 0;
        cfg.sim.seed = 12; cfg.sim.snapshotInterval = 10;
        cfg.sim.maskSaveInterval = 100; cfg.sim.stopEqDiameter = 70;

        S = ca.run_simulation(cfg);
        testCase.verifyGreaterThan(S.history.nCells(end), S.history.nCells(1), ...
            'The cluster must grow.');
        testCase.verifyGreaterThan(S.history.dEq(end), S.history.dEq(1), ...
            'The equivalent diameter must grow.');
    end

    function testDeterminism(testCase)
        cfg = testCase.cfg;
        cfg.sim.N = 128; cfg.sim.nc = 64; cfg.sim.maxSteps = 40;
        cfg.sim.nutrientCoupling = false;
        cfg.sim.T_cell = 2; cfg.sim.alpha = 0.6; cfg.sim.beta = 0.8;
        cfg.sim.P_move = 0.05; cfg.sim.A_i = 0.6;
        cfg.sim.delta = 1e-3; cfg.sim.pStarve = 0;
        cfg.sim.seed = 13; cfg.sim.snapshotInterval = 20;
        cfg.sim.maskSaveInterval = 100; cfg.sim.stopEqDiameter = 60;

        S1 = ca.run_simulation(cfg);
        S2 = ca.run_simulation(cfg);
        testCase.verifyTrue(isequal(S1.occ, S2.occ), ...
            'Same seed must produce an identical lattice.');
    end

end
end
