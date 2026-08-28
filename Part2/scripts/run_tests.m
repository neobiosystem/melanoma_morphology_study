function results = run_tests()
%RUN_TESTS Execute the Part 2 validation suite.
%
%   results = RUN_TESTS() runs every test under Part2/tests and prints a
%   summary. The suite must be green before any mechanism comparison is
%   trusted -- in particular TestIsotropy proves the lattice is not imprinting
%   a spurious anisotropy, and TestBridge proves the measurement bridge still
%   agrees with Part 1's own validation.

p2 = fileparts(fileparts(mfilename('fullpath')));
addpath(p2);
addpath(fullfile(p2, 'config'));
addpath(fullfile(p2, 'tests'));
p1 = fullfile(fileparts(p2), 'Part1');
addpath(p1);
addpath(fullfile(p1, 'config'));
addpath(fullfile(p1, 'tests'));

results = runtests(fullfile(p2, 'tests'));

disp(' ');
disp(table(results));

fprintf('\n%d passed, %d failed, %d incomplete (%.1f s)\n', ...
    nnz([results.Passed]), nnz([results.Failed]), ...
    nnz([results.Incomplete]), sum([results.Duration]));
end
