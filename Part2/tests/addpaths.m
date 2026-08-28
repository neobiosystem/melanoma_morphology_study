function addpaths()
%ADDPATHS Put Part 2 and Part 1 on the MATLAB path for tests and scripts.
%
%   Part 2's +ca and +bridge call Part 1's +morph, and the measurement bridge
%   calls Part 1's default_config. Those live in the sibling Part1/ tree and
%   are added here (never copied), so the alignment with Part 3 cannot drift.

p2 = fileparts(fileparts(mfilename('fullpath')));          % Part2/
p1 = fullfile(fileparts(p2), 'Part1');

addpath(p2);
addpath(fullfile(p2, 'config'));
addpath(fullfile(p2, 'tests'));
addpath(p1);
addpath(fullfile(p1, 'config'));
addpath(fullfile(p1, 'tests'));
end
