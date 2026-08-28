function addpaths()
%ADDPATHS Add Part 1/2/3 roots and config/tests dirs to the MATLAB path.
p3 = fileparts(fileparts(mfilename('fullpath')));   % Part3/
p2 = fullfile(fileparts(p3), 'Part2');
p1 = fullfile(fileparts(p3), 'Part1');
paths = {p3, fullfile(p3, 'config'), fullfile(p3, 'tests'), ...
         p2, fullfile(p2, 'config'), fullfile(p2, 'tests'), ...
         p1, fullfile(p1, 'config'), fullfile(p1, 'tests')};
addpath(paths{:});
end
