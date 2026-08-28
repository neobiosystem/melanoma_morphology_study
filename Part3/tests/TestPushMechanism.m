classdef TestPushMechanism < matlab.unittest.TestCase
%TESTPUSHMECHANISM The push-displacement arm conserves cells and stays exact.
%
%   Four properties: (1) push conserves nCells = init + div - death; (2) push
%   lets a surrounded shell divide, so totalDivisions is >= the no-push run;
%   (3) the simulator is deterministic; (4) with both new mechanisms off, the
%   lattice is byte-identical to Part 2's CA -- the "homogeneous limit" that
%   guarantees anisotropy from the new knobs is signal, not a grid artefact.

    methods (Test)
        function testConservationWithPush(testCase)
            cfg = small_cfg(true, 0.2, 1);
            S = ca3.run_simulation(cfg);
            expected = S.initCells + S.totalDivisions - S.totalDeaths;
            testCase.verifyEqual(S.nCells, expected, ...
                'Conservation violated: nCells = initCells + divisions - deaths');
        end

        function testPushEnablesBulkGrowth(testCase)
            % At a FIXED number of steps (no early stop), push lets a second
            % shell divide, so the tumour gains strictly more cells than the
            % surface-only (Eden) baseline.
            S_off = ca3.run_simulation(small_cfg(false, 0.2, 1));
            S_on  = ca3.run_simulation(small_cfg(true, 0.2, 1));
            testCase.verifyGreaterThan(S_on.nCells, S_off.nCells);
        end

        function testDeterministic(testCase)
            cfg = small_cfg(true, 0.2, 1);
            cfg.sim.seed = 5;
            S1 = ca3.run_simulation(cfg);
            S2 = ca3.run_simulation(cfg);
            testCase.verifyEqual(S1.occ, S2.occ);
        end

        function testHomogeneousLimitMatchesPart2(testCase)
            % displace=false AND heteroOn=false must reduce to Part 2 exactly.
            cfg = small_cfg(false, 0.95, 1);
            cfg.sim.displace = false;
            cfg.sim.heteroOn = false;
            S3 = ca3.run_simulation(cfg);
            S2 = ca.run_simulation(cfg);
            testCase.verifyEqual(S3.occ, S2.occ);
        end
    end
end

function cfg = small_cfg(displace, beta, pushRate)
% Fixed-step run: stopEqDiameter is set unreachably high so no early-stop
% overshoot confounds the comparison between runs of different growth rate.
cfg = part3_config();
cfg.sim.N  = 96;
cfg.sim.nc = 12;
cfg.sim.stopEqDiameter = 1000;
cfg.sim.maxSteps = 200;
cfg.sim.alpha = 0.5;
cfg.sim.beta  = beta;
cfg.sim.displace  = displace;
cfg.sim.pushRate  = pushRate;
cfg.sim.heteroOn  = false;
cfg.sim.seed = 2;
end
