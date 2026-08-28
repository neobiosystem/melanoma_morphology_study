function robustness_dual_route()
%ROBUSTNESS_DUAL_ROUTE Sign-reversal evidence: does the signature depend on route?
%
%   ROBUSTNESS_DUAL_ROUTE() regenerates the pairwise benign-vs-invasive effect
%   (rank-biserial correlation, sign = direction) for all seven descriptors on
%   BOTH segmentation routes, and renders them as a diverging bar so the
%   sign-reversal is visible rather than asserted. Bars coloured green point
%   the biologically expected direction; bars coloured red point the wrong way.
%
%   Outputs:
%     results/figures/route_sign_reversal.png  - the divergence figure
%     results/robustness_dual_route.csv        - the consolidated table
%
%   The numbers come from results/posthoc_all.csv (written by export_posthoc),
%   so the directions and p-values cannot drift from the saved analyses.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root, 'config'));
cfg = default_config();

PH = readtable(fullfile(cfg.paths.results, 'posthoc_all.csv'), 'TextType', 'string');
BI = PH(PH.groupA == "benign" & PH.groupB == "invasive", :);

% Biological expectation: does invasive move UP (+1) or DOWN (-1) vs benign?
expectedUp = struct('Df', 1, 'compactness', 1, 'circularity', -1, ...
    'asymmetry', 1, 'solidity', -1, 'convexity', -1, 'perimeterDecaySlope', -1);

names = cfg.stats.descriptors;
n = numel(names);

% Pull the otsu effect size to order the bars by strength (weakest last).
otsuRb = zeros(n, 1);
for i = 1:n
    row = BI(strcmp(BI.descriptor, names{i}) & strcmp(BI.route, 'otsu'), :);
    if ~isempty(row), otsuRb(i) = row.rankBiserial(1); end
end
[~, ord] = sort(abs(otsuRb), 'descend');
names = names(ord);

% Build the table rows.
rows = struct('descriptor', {}, 'expected', {}, 'expert_rb', {}, 'expert_pAdj', {}, ...
    'expert_dir', {}, 'otsu_rb', {}, 'otsu_pAdj', {}, 'otsu_dir', {});
expertRb = zeros(n, 1); expertP = zeros(n, 1); expertDir = zeros(n, 1);
otsuRb   = zeros(n, 1); otsuP   = zeros(n, 1); otsuDir   = zeros(n, 1);
for i = 1:n
    nm = names{i};
    exp = expectedUp.(nm);
    re = BI(strcmp(BI.descriptor, nm) & strcmp(BI.route, 'expert'), :);
    ro = BI(strcmp(BI.descriptor, nm) & strcmp(BI.route, 'otsu'), :);
    if isempty(re), continue; end

    expertRb(i) = re.rankBiserial(1); expertP(i) = re.pAdjusted(1);
    otsuRb(i)   = ro.rankBiserial(1); otsuP(i)   = ro.pAdjusted(1);
    expertDir(i) = double(sign(expertRb(i)) == exp);
    otsuDir(i)   = double(sign(otsuRb(i))   == exp);

    rows{end+1} = struct('descriptor', string(nm), 'expected', exp, ...
        'expert_rb', expertRb(i), 'expert_pAdj', expertP(i), 'expert_dir', expertDir(i), ...
        'otsu_rb', otsuRb(i), 'otsu_pAdj', otsuP(i), 'otsu_dir', otsuDir(i)); %#ok<AGROW>
end
T = struct2table([rows{:}]);

fprintf('Sign reversal by route (rank-biserial, benign -> invasive):\n');
fprintf('%-20s %4s | %9s %9s %8s | %9s %9s %8s\n', ...
    'descriptor', 'exp', 'expert_rb', 'p_adj', 'dir', 'otsu_rb', 'p_adj', 'dir');
fprintf('%s\n', repmat('-', 1, 72));
for i = 1:height(T)
    fprintf('%-20s %+4d | %9.3f %9.2g %8s | %9.3f %9.2g %8s\n', ...
        T.descriptor(i), T.expected(i), T.expert_rb(i), T.expert_pAdj(i), ...
        dirword(T.expert_dir(i)), T.otsu_rb(i), T.otsu_pAdj(i), dirword(T.otsu_dir(i)));
end

nExpertRev = nnz(T.expert_dir == 0);
nOtsuRev   = nnz(T.otsu_dir == 0);
fprintf('\n%d/%d reversed on the expert route; %d/%d on the automatic route.\n', ...
    nExpertRev, height(T), nOtsuRev, height(T));

% ---- Diverging bar figure --------------------------------------------
figDir = fullfile(cfg.paths.results, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

f = figure('Visible', 'off', 'Position', [100 100 820 520]);
hold on;
ymax = max(abs([T.expert_rb; T.otsu_rb])) * 1.15;
plot([0 0], [0.4 n+0.6], 'k-', 'LineWidth', 0.8);           % zero line
dy = 0.34;
for i = 1:height(T)
    y = n - i + 1;
    % expert bar (left half) and otsu bar (right half)
    cExp = [0.85 0.33 0.30];   % red
    cOts = [0.23 0.55 0.34];   % green
    if T.expert_dir(i), cExp = cOts; end          % expected -> green
    if ~T.otsu_dir(i),   cOts = [0.85 0.33 0.30]; end

    barh(y, T.expert_rb(i), 0.40, 'FaceColor', cExp, 'EdgeColor', 'none');
    barh(y - dy, T.otsu_rb(i), 0.40, 'FaceColor', cOts, 'EdgeColor', 'none');

    text(-sign(T.expert_rb(i))*ymax*0.5, y, sprintf('%.2f', T.expert_rb(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8);
    text(-sign(T.otsu_rb(i))*ymax*0.5, y - dy, sprintf('%.2f', T.otsu_rb(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8);
end
set(gca, 'YTick', (1:n) - dy/2, 'YTickLabel', T.descriptor, 'YTickLabelRotation', 0);
xlim([-ymax ymax]);
xlabel('rank-biserial, benign -> invasive (sign = direction)');
title('Between-grade effect by segmentation route: expert masks reverse the sign');
legend({'expert (polygon) mask', 'automatic (Otsu) mask'}, 'Location', 'southwest');
grid on; box on;
hold off;

saveas(f, fullfile(figDir, 'route_sign_reversal.png'));
saveas(f, fullfile(figDir, 'route_sign_reversal.pdf'));
close(f);

writetable(T, fullfile(cfg.paths.results, 'robustness_dual_route.csv'));
fprintf('Figure written to %s\n', figDir);
end

% -------------------------------------------------------------------------
function w = dirword(d)
if d, w = 'ok'; else, w = 'REV'; end
end
