function results = run_tests()
%RUN_TESTS Execute the Part 1 validation suite.
%
%   results = RUN_TESTS() runs every test under Part1/tests and prints a
%   summary. The suite must be green before any cohort result is trusted:
%   in particular TestBoxCounting/testFilledDiskIsDimensionTwo proves the
%   box counter measures boundaries rather than areas, and without it every
%   reported D_f would be meaningless.

root = fileparts(fileparts(mfilename('fullpath')));

addpath(root);                              % +morph, +seg, +stats, +agree
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'tests'));

results = runtests(fullfile(root, 'tests'));

disp(' ');
disp(table(results));

fprintf('\n%d passed, %d failed, %d incomplete (%.1f s)\n', ...
    nnz([results.Passed]), nnz([results.Failed]), ...
    nnz([results.Incomplete]), sum([results.Duration]));
end
