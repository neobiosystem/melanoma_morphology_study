function make_figures()
%MAKE_FIGURES Render the Part 3 report figures from analysis outputs.
%
%   Four figures:
%     1. Discovery trajectories: the two mechanism arms in (compactness,
%        solidity) space with the real benign/invasive medians overlaid.
%     2. Morphology tiling: the homogeneous limit, the push endpoint, and two
%        growth-field endpoints so the reader sees the shape each mechanism
%        actually produces.
%     3. Phase diagram: dist-to-melanoma over the two dominant knobs
%        (heteroAmp x heteroLambda) with the minimum annotated.
%     4. Sensitivity: OAT tornado and Sobol' indices side by side.
%
%   Figures are written to results/figures as PNG.

p3 = fileparts(fileparts(mfilename('fullpath')));
addpath(p3); addpath(fullfile(p3, 'config')); addpath(fullfile(p3, 'tests'));
p2 = fullfile(fileparts(p3), 'Part2');
addpath(p2); addpath(fullfile(p2, 'config'));
p1 = fullfile(fileparts(p3), 'Part1');
addpath(p1); addpath(fullfile(p1, 'config'));

outDir = fullfile(p3, 'results', 'figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end
res = fullfile(p3, 'results');
targets = jsondecode(fileread(fullfile(p1, 'results', 'calibration_targets_otsu.json')));

benC = targets.descriptors.compactness.benign.median;
benS = targets.descriptors.solidity.benign.median;
invC = targets.descriptors.compactness.invasive.median;
invS = targets.descriptors.solidity.invasive.median;

% ---- Figure 1: discovery trajectories --------------------------------
csvPath = fullfile(res, 'discovery_trajectories.csv');
if isfile(csvPath)
    T = readtable(csvPath, 'TextType', 'string');
    armNames = categories(categorical(T.arm));
    cmap = lines(numel(armNames));
    f1 = figure('Visible', 'off', 'Position', [100 100 760 560]);
    hold on;
    for a = 1:numel(armNames)
        nm = armNames{a};
        sub = T(T.arm == nm & T.accepted == 1, :);
        scatter(sub.clean_compactness, sub.clean_solidity, 40, cmap(a, :), ...
            'filled', 'DisplayName', char(nm));
        for k = 1:height(sub)
            text(sub.clean_compactness(k), sub.clean_solidity(k), ...
                sprintf('%d', sub.level(k)), 'FontSize', 7, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
    end
    scatter(benC, benS, 160, '^', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
        'DisplayName', 'real benign');
    scatter(invC, invS, 160, 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
        'DisplayName', 'real invasive');
    xlabel('compactness'); ylabel('solidity');
    title('Part 3 arms: growth field reaches the melanoma region (level numbers)');
    legend('Location', 'best'); grid on;
    saveas(f1, fullfile(outDir, 'discovery_trajectories.png'));
    saveas(f1, fullfile(outDir, 'discovery_trajectories.pdf'));
    close(f1);
end

% ---- Figure 2: morphology tiling --------------------------------------
cfgBase = part3_config();
cfgBase.sim.N = 384; cfgBase.sim.nc = 48;
cfgBase.sim.stopEqDiameter = round(0.74 * cfgBase.sim.N);
cfgBase.sim.maxSteps = 4000; cfgBase.sim.seed = 1;
R = cfgBase.sim.stopEqDiameter / 2;

tiles = { ...
    'homogeneous limit', struct('displace', false, 'heteroOn', false); ...
    'push (beta=0.05)', struct('displace', true, 'heteroOn', false, 'beta', 0.05); ...
    'growth field  l=R/2, a=0.5', struct('displace', false, 'heteroOn', true, ...
        'heteroLambda', R/2, 'heteroAmp', 0.5); ...
    'growth field  l=R/4, a=1.0', struct('displace', false, 'heteroOn', true, ...
        'heteroLambda', R/4, 'heteroAmp', 1.0) };

f2 = figure('Visible', 'off', 'Position', [100 100 1000 280]);
tl = tiledlayout(1, 4, 'Padding', 'tight', 'TileSpacing', 'compact');
for a = 1:size(tiles, 1)
    cfg = cfgBase;
    ov = tiles{a, 2};
    for f = fieldnames(ov)'
        cfg.sim.(f{1}) = ov.(f{1});
    end
    S = ca3.run_simulation(cfg);
    nexttile;
    imagesc(ca.to_mask(S)); axis image off; colormap(gray);
    title(strrep(tiles{a, 1}, '_', '\_'), 'Interpreter', 'tex');
end
title(tl, 'Part 3 mechanisms: growth field produces large lobes', 'FontWeight', 'bold');
saveas(f2, fullfile(outDir, 'mechanism_morphology.png'));
close(f2);

% ---- Figure 3: phase diagram ------------------------------------------
phasePath = fullfile(res, 'phase_spatial_heterogeneity.csv');
if isfile(phasePath)
    P = readtable(phasePath, 'TextType', 'string');
    p1 = 'heteroAmp'; p2 = 'heteroLambda';
    amp = unique(P.(p1)); lam = unique(P.(p2));
    [A, L] = ndgrid(amp, lam);
    D = nan(numel(amp), numel(lam));
    for i = 1:height(P)
        ia = find(amp == P.(p1)(i));
        il = find(lam == P.(p2)(i));
        D(ia, il) = P.dist(i);
    end
    f3 = figure('Visible', 'off', 'Position', [100 100 680 560]);
    imagesc(lam, amp, D); axis xy; colorbar;
    colormap(flipud(hot));
    set(gca, 'YDir', 'normal');
    hold on;
    [~, idx] = min(D(:));
    [ia, il] = ind2sub(size(D), idx);
    plot(lam(il), amp(ia), 'wo', 'MarkerSize', 14, 'LineWidth', 2);
    text(lam(il), amp(ia), sprintf('  min dist %.2f', D(idx)), ...
        'Color', 'w', 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    xlabel(sprintf('%s (px)', p2)); ylabel(p1);
    title('Distance to melanoma median vs. growth-field knobs');
    saveas(f3, fullfile(outDir, 'phase_diagram.png'));
    saveas(f3, fullfile(outDir, 'phase_diagram.pdf'));
    close(f3);
end

% ---- Figure 4: sensitivity --------------------------------------------
oatPath = fullfile(res, 'oat_spatial_heterogeneity.csv');
sobolPath = fullfile(res, 'sobol_spatial_heterogeneity.json');
f4 = figure('Visible', 'off', 'Position', [100 100 980 420]);
tl2 = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;   % OAT tornado
if isfile(oatPath)
    T = readtable(oatPath, 'TextType', 'string');
    base = mean(T.dist(strcmp(T.knob, 'baseline')), 'omitnan');
    knobs = {'heteroAmp', 'heteroLambda'};
    delta = zeros(numel(knobs), 1);
    for ki = 1:numel(knobs)
        idx = strcmp(T.knob, knobs{ki});
        lvls = unique(T.level(idx));
        mx = 0;
        for lv = lvls'
            d = mean(T.dist(idx & T.level == lv), 'omitnan');
            mx = max(mx, abs(d - base));
        end
        delta(ki) = mx;
    end
    barh(flip(delta)); set(gca, 'YTickLabel', flip(knobs));
    xlabel('max |\Deltadist| from baseline'); title('OAT screening');
    grid on;
end

nexttile;   % Sobol' bars
if isfile(sobolPath)
    SB = jsondecode(fileread(sobolPath));
    np = numel(SB.params);
    xx = 1:np;
    b1 = bar(xx - 0.15, SB.Si, 0.3); hold on;
    b2 = bar(xx + 0.15, SB.STi, 0.3);
    set(gca, 'XTick', xx, 'XTickLabel', SB.params);
    ylabel('index'); title('Sobol'' (dist response)');
    legend([b1 b2], {'first-order S_i', 'total S_Ti'}, 'Location', 'best');
    grid on;
end

title(tl2, 'Part 3 sensitivity: both growth-field knobs matter');
saveas(f4, fullfile(outDir, 'sensitivity.png'));
saveas(f4, fullfile(outDir, 'sensitivity.pdf'));
close(f4);

fprintf('Figures written to %s\n', outDir);
end
