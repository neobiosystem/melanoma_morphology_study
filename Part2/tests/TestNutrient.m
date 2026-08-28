classdef TestNutrient < matlab.unittest.TestCase
%TESTNUTRIENT The steady-state field solver against known behaviour.
%
%   Two properties are checked: with no consumption the Dirichlet solution
%   is the constant C0 field, and with uniform consumption the field is
%   symmetric and depressed toward the centre. These pin down the sparse
%   solve's Dirichlet handling and its consumption coefficient without
%   needing a full closed-form 2-D solution.

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

    function testNoConsumptionIsUniform(testCase)
        cfg = testCase.cfg;
        rho = zeros(cfg.sim.nc);
        C = ca.solve_nutrient(rho, cfg);
        testCase.verifyEqual(C, cfg.sim.C0 * ones(size(rho)), 'AbsTol', 1e-6, ...
            'Zero consumption must yield the constant Dirichlet field C0.');
    end

    function testUniformConsumptionDepressesCenter(testCase)
        cfg = testCase.cfg;
        rho = 0.5 * ones(cfg.sim.nc);
        C = ca.solve_nutrient(rho, cfg);

        % Centre must be depleted relative to the far-field boundary.
        c = round(cfg.sim.nc / 2);
        testCase.verifyLessThan(C(c, c), cfg.sim.C0 - 1e-9, ...
            'Consumption must depress the field below C0 at the centre.');

        % The solution is symmetric about both midlines (square domain).
        testCase.verifyEqual(C, fliplr(C), 'AbsTol', 1e-9, ...
            'Left-right symmetry expected on a square domain.');
        testCase.verifyEqual(C, flipud(C), 'AbsTol', 1e-9, ...
            'Top-bottom symmetry expected on a square domain.');
    end

end
end
