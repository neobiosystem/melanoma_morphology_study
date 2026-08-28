function run_tests()
%RUN_TESTS Run the Part 3 validation suite.
%   run_tests() runs the class-based tests under Part3/tests and errors if any
%   fail. The "homogeneous limit" and "distance hand-computation" tests are the
%   load-bearing ones: they prove the new mechanisms vanish cleanly and the
%   alignment metric is exact.
import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.TestRunProgressPlugin;

p3 = fileparts(fileparts(mfilename('fullpath')));
p2 = fullfile(fileparts(p3), 'Part2');
p1 = fullfile(fileparts(p3), 'Part1');
addpath(p3); addpath(fullfile(p3, 'config')); addpath(fullfile(p3, 'tests'));
addpath(p2); addpath(fullfile(p2, 'config')); addpath(fullfile(p2, 'tests'));
addpath(p1); addpath(fullfile(p1, 'config')); addpath(fullfile(p1, 'tests'));

suite = TestSuite.fromFolder(fullfile(p3, 'tests'), 'IncludingSubfolders', false);
runner = TestRunner.withTextOutput();
results = runner.run(suite);

nPass = sum([results.Passed]);
nTot  = numel(results);
fprintf('\nPart 3 tests: %d/%d passed.\n', nPass, nTot);
if nPass < nTot
    error('run_tests:failures', 'Part 3 tests failed (%d/%d).', nPass, nTot);
end
end
