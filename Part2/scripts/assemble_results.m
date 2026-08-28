function assemble_results()
%ASSEMBLE_RESULTS Merge per-arm chunk CSVs into the combined result set.
%
%   ASSEMBLE_RESULTS() reads every results/arm_*.csv written by the chunked
%   mechanism_comparison runs, concatenates and sorts them, writes
%   results/mechanism_trajectories.csv, and rebuilds results/part2_handoff.json
%   via +comp/write_handoff. Run it once after all four arm chunks finish.

p2 = fileparts(fileparts(mfilename('fullpath')));
addpath(p2); addpath(fullfile(p2, 'config'));
p1 = fullfile(fileparts(p2), 'Part1');
addpath(p1); addpath(fullfile(p1, 'config'));

res = fullfile(p2, 'results');
files = dir(fullfile(res, 'arm_*.csv'));
if isempty(files)
    error('No arm_*.csv chunks found in %s', res);
end

T = [];
for k = 1:numel(files)
    sub = readtable(fullfile(res, files(k).name), 'TextType', 'string');
    T = [T; sub]; %#ok<AGROW>
end
T = sortrows(T, {'arm', 'level', 'seed'});

cfg0 = sim_config();
writetable(T, fullfile(res, 'mechanism_trajectories.csv'));
comp.write_handoff(T, cfg0.sim.N, cfg0.sim.stopEqDiameter, p1, ...
    fullfile(res, 'part2_handoff.json'));

fprintf('Assembled %d rows from %d arm chunk(s) -> mechanism_trajectories.csv + part2_handoff.json\n', ...
    height(T), numel(files));
end
