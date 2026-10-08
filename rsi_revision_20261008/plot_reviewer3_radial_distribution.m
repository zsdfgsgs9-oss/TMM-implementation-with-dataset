function plot_reviewer3_radial_distribution()
% Reproduce the 12.67 kV image-plane radial distributions from paired particles.
root = fileparts(mfilename('fullpath'));
set(groot, 'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman');
ci = readtable(fullfile(root, 'RSI第四版结论数据报告归档_20260916', ...
    '04_分项数据', 'D1_统计区间', 'baseline_bootstrap.csv'), 'TextType', 'string');
tmm_file = fullfile(root, 'RSI第三版任务1-5成果归档_20260915', ...
    '02_任务1-3_MATLAB算子修正与回归', ...
    'table_iii_corrected_tmm_coordinates_with_pidx.csv');
comsol_root = fullfile(root, 'RSI大修仿真统一归档_20260911', ...
    '03_COMSOL完整结果', 'revision_comsol', ...
    'baseline_three_model_20260904');
t = readtable(tmm_file, 'TextType', 'string');
t = t(t.case_name == "01_intrinsic_error", :);
nr = readmatrix(fullfile(comsol_root, 'NR_full', ...
    'poincare_zrel_1000um_with_pidx.csv'), 'FileType', 'text');
rel = readmatrix(fullfile(comsol_root, 'REL_full', ...
    'poincare_zrel_1000um_with_pidx.csv'), 'FileType', 'text');
nr = nr(all(isfinite(nr(:, [1 2 4])), 2), :);
rel = rel(all(isfinite(rel(:, [1 2 4])), 2), :);
[ids, ia, ib] = intersect(round(t.pidx), round(nr(:, 4)), 'sorted');
[ids2, ~, ic] = intersect(ids, round(rel(:, 4)), 'stable');
assert(isequal(ids, ids2) && numel(ids) == 50000, ...
    'The three models must share the same 50,000 particle IDs.');
xy = {[t.x_nm(ia), t.y_nm(ia)], nr(ib, 1:2) * 1e3, ...
    rel(ic, 1:2) * 1e3};
names = {'NR-TMM', 'NR-full', 'REL-full'};
colors = [0.12 0.32 0.67; 0.85 0.33 0.13; 0.11 0.55 0.42];
styles = {'-', '--', '-.'};
r = cell(1, 3);
q90 = zeros(1, 3); q99 = zeros(1, 3);
for k = 1:3
    center = mean(xy{k}, 1);
    r{k} = hypot(xy{k}(:, 1) - center(1), xy{k}(:, 2) - center(2));
    q90(k) = sample_quantile(r{k}, 0.90);
    q99(k) = sample_quantile(r{k}, 0.99);
end
assert(abs(2*q99(1) - 39.1508460633318) < 1e-5);
assert(abs(2*q99(2) - 39.2218841700284) < 1e-5);
assert(abs(2*q99(3) - 37.9058818088727) < 1e-5);
fprintf('Model       FW90 (nm)   FW99 (nm)\n');
for k = 1:3
    fprintf('%-9s %9.4f %11.4f\n', names{k}, 2*q90(k), 2*q99(k));
end

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'inches', 'Position', [0.5 0.5 7.2 3.35]);
layout = tiledlayout(fig, 5, 3, 'TileSpacing', 'compact', 'Padding', 'loose');

ax = nexttile(layout, 1, [4 1]); hold(ax, 'on');
edges = 0:0.5:ceil(max(cellfun(@max, r)) + 0.5);
for k = 1:3
    c = histcounts(r{k}, edges);
    stairs(ax, edges, [c, c(end)] / (numel(r{k}) * 0.5), ...
        'Color', colors(k, :), 'LineStyle', styles{k}, ...
        'LineWidth', 1.45, 'DisplayName', names{k});
end
xlim(ax, [0 23]); ylim(ax, [0 0.10]);
xlabel(ax, 'Radius (nm)'); ylabel(ax, 'Radial PDF (nm^{-1})');
title(ax, '(a) 0.50-nm bins');
grid(ax, 'on');

ax = nexttile(layout, 2, [4 1]); hold(ax, 'on');
for k = 1:3
    sr = sort(r{k});
    stairs(ax, [sr(1); sr], [0; (1:numel(sr))' / numel(sr)], ...
        'Color', colors(k, :), 'LineStyle', styles{k}, 'LineWidth', 1.4);
    plot(ax, q90(k), 0.90, 'o', 'Color', colors(k, :), ...
        'MarkerFaceColor', colors(k, :), 'MarkerSize', 3.5);
    plot(ax, q99(k), 0.99, 's', 'Color', colors(k, :), ...
        'MarkerFaceColor', colors(k, :), 'MarkerSize', 3.5);
    for p = [90 99]
        row = ci(ci.metric == "FW" + p & ci.model_a == names{k}, :);
        if ~isempty(row)
            limits = [row.model_a_ci95_low_nm(1) row.model_a_ci95_high_nm(1)]/2;
        else
            row = ci(ci.metric == "FW" + p & ci.model_b == names{k}, :);
            assert(~isempty(row));
            limits = [row.model_b_ci95_low_nm(1) row.model_b_ci95_high_nm(1)]/2;
        end
        y = p/100;
        plot(ax, limits, [y y], 'Color', colors(k,:), 'LineWidth', 1.0);
        plot(ax, limits(1)*[1 1], y+[-.002 .002], 'Color', colors(k,:));
        plot(ax, limits(2)*[1 1], y+[-.002 .002], 'Color', colors(k,:));
    end
end
xlim(ax, [15 22]); ylim(ax, [0.85 1.005]);
yticks(ax, [0.85 0.90 0.95 0.99 1.00]);
xlabel(ax, 'Radius (nm)'); ylabel(ax, 'Empirical CDF');
title(ax, '(b) CDF tail and 95% CIs');
grid(ax, 'on');

ax = nexttile(layout, 3, [4 1]); hold(ax, 'on');
widths = [0.25 0.50 1.00];
bin_colors = [0.60 0.16 0.47; 0.00 0.58 0.70; 0.20 0.20 0.20];
bin_styles = {':', '-.', '-'};
for k = 1:numel(widths)
    w = widths(k);
    bin_edges = 0:w:(ceil(max(r{2}) / w) * w + w);
    c = histcounts(r{2}, bin_edges);
    stairs(ax, bin_edges, [c, c(end)] / (numel(r{2}) * w), ...
        'Color', bin_colors(k, :), 'LineStyle', bin_styles{k}, ...
        'LineWidth', 1.3, 'DisplayName', sprintf('%.2f', w));
end
xlim(ax, [15 22]); ylim(ax, [0 0.11]);
xlabel(ax, 'Radius (nm)'); ylabel(ax, 'Radial PDF (nm^{-1})');
title(ax, '(c) NR-full bin sensitivity');
grid(ax, 'on');

% Dedicated key tiles identify the panels without covering the data.
key_ax = nexttile(layout, 13, [1 2]); hold(key_ax, 'on'); axis(key_ax, 'off');
for k = 1:3
    plot(key_ax, nan, nan, 'Color', colors(k,:), ...
        'LineStyle', styles{k}, 'LineWidth', 1.45, 'DisplayName', names{k});
end
model_legend = legend(key_ax, 'Location', 'north', 'Orientation', 'horizontal', ...
    'NumColumns', 3, 'FontSize', 8, 'Box', 'off');
model_legend.Title.String = 'Panels (a), (b): transport model';
model_legend.Title.FontWeight = 'normal';

key_ax = nexttile(layout, 15); hold(key_ax, 'on'); axis(key_ax, 'off');
for k = 1:3
    plot(key_ax, nan, nan, 'Color', bin_colors(k,:), ...
        'LineStyle', bin_styles{k}, 'LineWidth', 1.3, ...
        'DisplayName', sprintf('%.2f', widths(k)));
end
bin_legend = legend(key_ax, 'Location', 'north', 'Orientation', 'horizontal', ...
    'NumColumns', 3, 'FontSize', 8, 'Box', 'off');
bin_legend.ItemTokenSize = [14 10];
bin_legend.Title.String = 'Panel (c): bin width (nm)';
bin_legend.Title.FontWeight = 'normal';
for ax = findall(fig, 'Type', 'axes')'
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 8.5, ...
        'LineWidth', 0.65, 'Box', 'on');
end
set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
outbase = fullfile(root, 'figures', 'fig_radial_distribution_reviewer3');
exportgraphics(fig, [outbase '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [outbase '.png'], 'Resolution', 300);
close(fig);
fprintf('Saved %s.pdf and .png\n', outbase);
end

function q = sample_quantile(r, p)
sr = sort(r(:));
pos = 1 + (numel(sr) - 1) * p;
lo = floor(pos); hi = ceil(pos);
q = sr(lo) + (pos - lo) * (sr(hi) - sr(lo));
end
