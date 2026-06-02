% export_figs_600dpi.m
% 打开 FIG，设置为期刊目标页宽，统一字体 ≥10pt，导出 600 DPI
%   - 双栏/全页宽: 7.0 英寸 (178 mm)
%   - 单栏宽:      3.5 英寸 (89 mm)
clear; close all;

target_dir = 'exported_figs';
if ~exist(target_dir, 'dir')
    mkdir(target_dir);
end

% {文件名, 目标宽度(英寸)}
% 这三栏对比图建议全页宽; 密度分布图可用全页或单栏
fig_config = {
    '聚焦图像示意图.fig',   7.0   % 三栏子图, 建议全页宽
    '5e7密度分布.fig',      7.0   % 三栏子图, 全页宽
    '5e4密度分布.fig',      7.0   % 三栏子图, 全页宽
    '系统级畸变对比.fig',   7.0   % 三栏子图, 全页宽
    };

for i = 1:size(fig_config, 1)
    fname     = fig_config{i, 1};
    tgt_width = fig_config{i, 2};  % 英寸

    fprintf('处理: %s (目标宽度 %.1f in)\n', fname, tgt_width);

    fig = openfig(fname, 'invisible');
    fig.Visible = 'off';

    % === 1. 将 figure 设为目标物理尺寸 ===
    orig_pos = fig.Position;  % [left bottom width height] in pixels
    aspect = orig_pos(4) / orig_pos(3);  % height/width
    tgt_height = tgt_width * aspect;

    set(fig, 'Units', 'inches', ...
             'Position', [1, 1, tgt_width, tgt_height], ...
             'PaperUnits', 'inches', ...
             'PaperSize', [tgt_width, tgt_height], ...
             'PaperPosition', [0, 0, tgt_width, tgt_height]);

    fprintf('  Figure: %.2f x %.2f in (宽高比 %.2f)\n', tgt_width, tgt_height, aspect);

    % === 2. 统一字体 ≥ 10pt ===
    ax_list = findobj(fig, 'Type', 'axes');
    for j = 1:length(ax_list)
        ax = ax_list(j);
        if ax.FontSize < 10
            ax.FontSize = 10;
        end
        if ax.XLabel.FontSize < 10
            ax.XLabel.FontSize = 10;
        end
        if ax.YLabel.FontSize < 10
            ax.YLabel.FontSize = 10;
        end
        if ~isempty(ax.Title.String) && ax.Title.FontSize < 11
            ax.Title.FontSize = 11;
        end
        ax.LineWidth = 0.8;
    end

    % 图例 ≥ 9pt
    lg_list = findobj(fig, 'Type', 'legend');
    for j = 1:length(lg_list)
        if lg_list(j).FontSize < 9
            lg_list(j).FontSize = 9;
        end
    end

    % === 3. 消除白边 ===
    % 对每个 axes 设置紧凑的 OuterPosition
    n_axes = length(ax_list);
    if n_axes == 3
        for j = 1:n_axes
            ax = ax_list(j);
            left = 0.04 + (j-1) * 0.33;
            ax.OuterPosition = [left, 0.08, 0.30, 0.88];
        end
    end

    % === 4. 导出 ===
    [~, base_name] = fileparts(fname);

    % PNG — 600 DPI
    png_path = fullfile(target_dir, [base_name '.png']);
    print(fig, png_path, '-dpng', '-r600', '-opengl');
    fprintf('  -> %s\n', png_path);

    % TIFF — 600 DPI
    tif_path = fullfile(target_dir, [base_name '.tif']);
    print(fig, tif_path, '-dtiff', '-r600', '-opengl');
    fprintf('  -> %s\n', tif_path);

    % PDF — 矢量 (线条/文字无损, 但散点图可能文件大)
    try
        pdf_path = fullfile(target_dir, [base_name '.pdf']);
        set(fig, 'Renderer', 'painters');
        exportgraphics(fig, pdf_path, 'ContentType', 'vector');
        fprintf('  -> %s\n', pdf_path);
    catch
        fprintf('  [跳过] PDF矢量导出失败(散点太多), 请用TIFF\n');
    end
    set(fig, 'Renderer', 'opengl');

    close(fig);
end

fprintf('\n===== 字体验证 =====\n');
fprintf('目标页宽 %.1f in → 导出的图中 10pt 字体在论文中即 10pt\n', tgt_width);
fprintf('若需单栏宽(3.5in): 将 config 中 7.0 改为 3.5 重跑\n');
fprintf('期刊投稿推荐: TIFF (600 DPI, 无损, DPI元数据正确)\n');
fprintf('文件在 %s/\n', target_dir);
