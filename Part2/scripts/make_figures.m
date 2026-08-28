function make_figures()
%MAKE_FIGURES Render the Part 2 report figures from analysis outputs.
%
%   Two figures:
%     1. Arm trajectories in (compactness, solidity) with the real benign and
%        invasive medians overlaid, read from mechanism_trajectories.csv.
%     2. A morphology tiling: one short simulation per arm at its malignant
%        endpoint, so the reader sees the shape each mechanism produces.
%
%   Figures are written to results/figures as PNG (and PDF where useful).

p2 = fileparts(fileparts(mfilename('fullpath')));
addpath(p2); addpath(fullfile(p2, 'config')); addpath(fullfile(p2, 'tests'));
p1 = fullfile(fileparts(p2), 'Part1');
addpath(p1); addpath(fullfile(p1, 'config')); addpath(fullfile(p1, 'tests'));

outDir = fullfile(p2, 'results', 'figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end

csvPath = fullfile(p2, 'results', 'mechanism_trajectories.csv');
targets = jsondecode(fileread(fullfile(p1, 'results', 'calibration_targets_otsu.json')));

% ---- Figure 1: trajectories -------------------------------------------
if isfile(csvPath)
    T = readtable(csvPath, 'TextType', 'string');
    armNames = categories(categorical(T.arm));

    benC  = targets.descriptors.compactness.benign.median;
    benS  = targets.descriptors.solidity.benign.median;
    invC  = targets.descriptors.compactness.invasive.median;
    invS  = targets.descriptors.solidity.invasive.median;

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
    scatter(benC, benS, 140, '^', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
        'DisplayName', 'real benign');
    scatter(invC, invS, 140, 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
        'DisplayName', 'real invasive');
    xlabel('compactness'); ylabel('solidity');
    title('Mechanism arms: benign -> malignant (level numbers)');
    legend('Location', 'best'); grid on;
    saveas(f1, fullfile(outDir, 'mechanism_trajectories.png'));
    saveas(f1, fullfile(outDir, 'mechanism_trajectories.pdf'));
    close(f1);
end

% ---- Figure 2: morphology tiling ---------------------------------------
endpoints = { ...
    'contact_inhibition', struct('alpha',0.4,'beta',0.05,'P_move',0.02,'A_i',0.90, ...
        'nutrientCoupling',false,'pStarve',0,'subclone_frac',0,'subclone_alphaMult',3); ...
    'nutrient_fingering', struct('alpha',0.6,'beta',0.95,'P_move',0.02,'A_i',0.90, ...
        'nutrientCoupling',true,'pStarve',0.01,'L',20,'subclone_frac',0,'subclone_alphaMult',3); ...
    'enhanced_migration', struct('alpha',0.4,'beta',0.95,'P_move',0.35,'A_i',0.05, ...
        'nutrientCoupling',false,'pStarve',0,'subclone_frac',0,'subclone_alphaMult',3); ...
    'clonal_heterogeneity', struct('alpha',0.4,'beta',0.95,'P_move',0.02,'A_i',0.90, ...
        'nutrientCoupling',false,'pStarve',0,'subclone_frac',0.40,'subclone_alphaMult',3) };

f2 = figure('Visible', 'off', 'Position', [100 100 1000 280]);
tl = tiledlayout(1, 4, 'Padding', 'tight', 'TileSpacing', 'compact');
for a = 1:size(endpoints, 1)
    cfg = sim_config();
    cfg.sim.N = 384; cfg.sim.nc = 48;
    cfg.sim.stopEqDiameter = round(0.74 * cfg.sim.N);
    cfg.sim.maxSteps = 4000; cfg.sim.seed = 42;
    ov = endpoints{a, 2};
    for f = fieldnames(ov)'
        if strcmp(f{1}, 'subclone_frac')
            cfg.sim.subclone.frac = ov.(f{1});
        elseif strcmp(f{1}, 'subclone_alphaMult')
            cfg.sim.subclone.alphaMult = ov.(f{1});
        else
            cfg.sim.(f{1}) = ov.(f{1});
        end
    end
    S = ca.run_simulation(cfg);
    nexttile;
    imagesc(ca.to_mask(S)); axis image off; colormap(gray);
    title(strrep(endpoints{a, 1}, '_', '\_'), 'Interpreter', 'tex');
end
title(tl, 'Final morphology per mechanism (malignant endpoint)', 'FontWeight', 'bold');
saveas(f2, fullfile(outDir, 'mechanism_morphology.png'));
close(f2);

fprintf('Figures written to %s\n', outDir);
end
