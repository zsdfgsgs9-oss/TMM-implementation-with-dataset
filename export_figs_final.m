% export_figs_final.m
% 重新排版 FIG 并导出 600 DPI
%  1. 密度图: 5e4/5e7 合并 ab 图, 保留 alpha + U0, 删除 spatial
%  2. 聚焦图像: 删除重复图例, 双行标题防重叠
%  3. 畸变对比: 3子图竖排
clear; close all;

% === 字体: JAP 接受 Helvetica 或 Times New Roman, Helvetica 渲染最稳定 ===
% 如需 Times New Roman: 改为 true (但 TeX 符号如 \alpha 可能显示异常)
use_times = false;
font_name = 'Helvetica';
if use_times
    font_name = 'Times New Roman';
end

target_dir = 'exported_figs';
if ~exist(target_dir, 'dir')
    mkdir(target_dir);
end

%% ========================================================================
%%  任务 1: 密度图 — 删alpha子图 + 合并为ab图
%% ========================================================================
fprintf('=== 处理密度分布图 ===\n');

% 打开两张源图, 在内存中删除 alpha 子图, 提取剩余子图数据
[fig_5e4, data_5e4] = load_and_clean_density('5e4密度分布.fig');
[fig_5e7, data_5e7] = load_and_clean_density('5e7密度分布.fig');

% 创建合并图: 2×2 (上排 alpha, 下排 U0; 左 5e4, 右 5e7)
fig_density = figure('Name', '密度分布对比', 'Visible', 'off');
set(fig_density, 'Units', 'inches', 'Position', [1 1 7.0 5.5]);

pos = {
    [0.08 0.55 0.41 0.40]   % (a) 左上: 5e4 alpha
    [0.55 0.55 0.41 0.40]   % (a) 右上: 5e4 U0
    [0.08 0.08 0.41 0.40]   % (b) 左下: 5e7 alpha
    [0.55 0.08 0.41 0.40]   % (b) 右下: 5e7 U0
    };
titles = {
    '\alpha distribution (5\times10^4)', 'U_0 distribution (5\times10^4)', ...
    '\alpha distribution (5\times10^7)', 'U_0 distribution (5\times10^7)'
    };
sources = {data_5e4.alpha, data_5e4.u0, data_5e7.alpha, data_5e7.u0};

for i = 1:4
    ax_new = axes('Parent', fig_density, 'OuterPosition', pos{i});
    copy_axes_content(ax_new, sources{i});
    title(ax_new, titles{i}, 'FontSize', 14);
    ax_new.FontSize = 13;
    ax_new.XLabel.FontSize = 13;
    ax_new.YLabel.FontSize = 13;
end

% (a) 上行 (5e4), (b) 下行 (5e7)
annotation(fig_density, 'textbox', [0.04, 0.92, 0.06, 0.04], ...
    'String', '(a)', 'FontSize', 15, 'FontWeight', 'bold', ...
    'EdgeColor', 'none', 'VerticalAlignment', 'top');
annotation(fig_density, 'textbox', [0.04, 0.45, 0.06, 0.04], ...
    'String', '(b)', 'FontSize', 15, 'FontWeight', 'bold', ...
    'EdgeColor', 'none', 'VerticalAlignment', 'top');

export_one(fig_density, fullfile(target_dir, '密度分布对比'), 7.0, font_name);
close(fig_density);
close(fig_5e4);
close(fig_5e7);

%% ========================================================================
%%  任务 2: 聚焦图像示意图 — 图例移到底部, 标题对齐
%% ========================================================================
fprintf('\n=== 处理聚焦图像示意图 ===\n');

fig_focus = openfig('聚焦图像示意图.fig', 'invisible');
fig_focus.Visible = 'off';
set(fig_focus, 'Units', 'inches', 'Position', [1 1 7.0 3.6]);

ax_focus = flipud(findobj(fig_focus, 'Type', 'axes'));
n_ax_f = min(length(ax_focus), 3);

% 收集全局数据范围, 统一坐标轴
x_all = []; y_all = [];
for j = 1:n_ax_f
    kids = ax_focus(j).Children;
    for c = 1:length(kids)
        if isa(kids(c), 'matlab.graphics.chart.primitive.Scatter')
            x_all = [x_all, kids(c).XData];
            y_all = [y_all, kids(c).YData];
        end
    end
end
margin = 0.08;
x_range = [min(x_all)-margin*range(x_all), max(x_all)+margin*range(x_all)];
y_range = [min(y_all)-margin*range(y_all), max(y_all)+margin*range(y_all)];

% 3子图横排, 右侧留更多空间防标题截断
ax_positions = {
    [0.03 0.16 0.28 0.70]   % 左: Initial
    [0.34 0.16 0.28 0.70]   % 中: TMM
    [0.65 0.16 0.28 0.70]   % 右: FEM (右边到0.93, 留0.07给标题)
    };

for j = 1:n_ax_f
    ax = ax_focus(j);
    ax.OuterPosition = ax_positions{j};
    ax.FontSize = 13;
    ax.XLabel.FontSize = 13;
    ax.YLabel.FontSize = 13;

    % 统一坐标轴范围 (先设范围, 再用 daspect 锁定比例)
    ax.XLim = x_range;
    ax.YLim = y_range;
    daspect(ax, [1 1 1]);

    % 长标题分两行, 保留括号
    t_str = ax.Title.String;
    if contains(t_str, 'Final electron distribution')
        ax.Title.String = strrep(t_str, ' (', [newline '(']);
        ax.Title.FontSize = 12;
    else
        ax.Title.FontSize = 13;
    end

    kids = ax.Children;
    for c = 1:length(kids)
        if isa(kids(c), 'matlab.graphics.chart.primitive.Scatter')
            % 降采样降低密度, 避免遮挡 FWxx 圆
            xd = kids(c).XData; yd = kids(c).YData;
            step = max(1, round(length(xd) / 800));
            idx = 1:step:length(xd);
            kids(c).XData = xd(idx);
            kids(c).YData = yd(idx);
            kids(c).SizeData = 22;
        end
        % FW50/FW90/FW99 圆: 点划线
        if isa(kids(c), 'matlab.graphics.chart.primitive.Line')
            if kids(c).LineWidth < 2
                kids(c).LineWidth = 1.5;
            end
            % FW99 改为点划线
            if contains(get(kids(c), 'DisplayName'), 'FW99') || ...
               (isempty(get(kids(c), 'DisplayName')) && kids(c).LineStyle == '-')
                % 最外圈 (FW99): 点划线, 较粗
                kids(c).LineStyle = '-.';
                kids(c).LineWidth = 1.8;
            end
        end
    end
end

% 强制所有子图的 Position 高度/宽度/bottom 对齐 (横排)
ref_pos = ax_focus(1).Position;
for j = 2:n_ax_f
    cur = ax_focus(j).Position;
    ax_focus(j).Position = [cur(1), ref_pos(2), ref_pos(3), ref_pos(4)];
end

% 添加 (a)(b)(c) 标签 — 按子图标题排序, 放左上角内侧
abc = {'(a)', '(b)', '(c)'};
% 按 OuterPosition left 排序, 确保左→右顺序
lefts = zeros(1, n_ax_f);
for j = 1:n_ax_f; lefts(j) = ax_focus(j).OuterPosition(1); end
[~, sort_idx] = sort(lefts);
label_y_f = ax_focus(sort_idx(1)).OuterPosition(2) + ax_focus(sort_idx(1)).OuterPosition(4) + 0.015;
for j = 1:n_ax_f
    ax_j = ax_focus(sort_idx(j));
    op = ax_j.OuterPosition;
    annotation(fig_focus, 'textbox', ...
        [op(1), label_y_f, 0.06, 0.03], ...
        'String', abc{j}, 'FontSize', 15, 'FontWeight', 'bold', ...
        'EdgeColor', 'none', 'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'left');
end

% 删除重复图例, 横向居中于底部
lg_list = findobj(fig_focus, 'Type', 'legend');
if length(lg_list) >= 2
    delete(lg_list(2));
end
lg = lg_list(1);
lg.FontSize = 10;
lg.Orientation = 'horizontal';
lg.Box = 'on';
lg.Position = [0.30, 0.03, 0.40, 0.05];

export_one(fig_focus, fullfile(target_dir, '聚焦图像示意图'), 7.0, font_name);
close(fig_focus);

%% ========================================================================
%%  任务 3: 系统级畸变对比 — 竖排
%% ========================================================================
fprintf('\n=== 处理系统级畸变对比 ===\n');

fig_dist = openfig('系统级畸变对比.fig', 'invisible');
fig_dist.Visible = 'off';
set(fig_dist, 'Units', 'inches', 'Position', [1 1 7.0 8.0]);

ax_dist = flipud(findobj(fig_dist, 'Type', 'axes'));
n_ad = min(length(ax_dist), 3);

% 收集全局数据范围
x_all = []; y_all = [];
for j = 1:n_ad
    kids = ax_dist(j).Children;
    for c = 1:length(kids)
        if isa(kids(c), 'matlab.graphics.chart.primitive.Scatter')
            x_all = [x_all, kids(c).XData];
            y_all = [y_all, kids(c).YData];
        end
    end
end
margin = 0.10;
x_range_d = [min(x_all)-margin*range(x_all), max(x_all)+margin*range(x_all)];
y_range_d = [min(y_all)-margin*range(y_all), max(y_all)+margin*range(y_all)];

% 3子图竖排, 完全相同的宽度和高度
vert_pos = {
    [0.10 0.68 0.85 0.26]   % 上: TMM
    [0.10 0.37 0.85 0.26]   % 中: FEM
    [0.10 0.06 0.85 0.26]   % 下: comparison
    };

for j = 1:n_ad
    ax = ax_dist(j);
    ax.OuterPosition = vert_pos{j};
    ax.FontSize = 13;
    ax.XLabel.FontSize = 13;
    ax.YLabel.FontSize = 13;
    ax.Title.FontSize = 14;

    % 统一坐标轴范围 (先设范围, 再用 daspect 锁定比例)
    ax.XLim = x_range_d;
    ax.YLim = y_range_d;
    daspect(ax, [1 1 1]);

    kids = ax.Children;
    for c = 1:length(kids)
        if isa(kids(c), 'matlab.graphics.chart.primitive.Scatter')
            if kids(c).SizeData < 20
                kids(c).SizeData = 28;
            end
        end
    end
end

% 所有子图: 降采样 + 统一标记风格 (TMM=蓝o, FEM=红x)
n_sample = 500;
for j = 1:n_ad
    kids = ax_dist(j).Children;
    scatters = [];
    for c = 1:length(kids)
        if isa(kids(c), 'matlab.graphics.chart.primitive.Scatter')
            scatters = [scatters, kids(c)];
        end
    end
    n_s = length(scatters);
    if n_s == 0; continue; end

    % 降采样
    for s = 1:n_s
        xd = scatters(s).XData; yd = scatters(s).YData;
        step = max(1, round(length(xd) / n_sample));
        idx = 1:step:length(xd);
        scatters(s).XData = xd(idx);
        scatters(s).YData = yd(idx);
    end

    if n_s == 1
        % 单数据集子图 (a)或(b)
        t = ax_dist(j).Title.String;
        if contains(t, 'TMM')
            scatters(1).Marker = 'o';
            scatters(1).MarkerEdgeColor = [0.15 0.35 0.75];
            scatters(1).MarkerFaceColor = 'none';
            scatters(1).SizeData = 50;
            scatters(1).LineWidth = 0.6;
        else
            scatters(1).Marker = 'x';
            scatters(1).MarkerEdgeColor = [0.85 0.25 0.25];
            scatters(1).SizeData = 28;
            scatters(1).LineWidth = 0.8;
        end
    elseif n_s >= 2
        % 双数据集子图 (c): TMM=蓝o, FEM=红x
        scatters(1).Marker = 'o';
        scatters(1).MarkerEdgeColor = [0.15 0.35 0.75];
        scatters(1).MarkerFaceColor = 'none';
        scatters(1).SizeData = 50;
        scatters(1).LineWidth = 0.6;
        scatters(2).Marker = 'x';
        scatters(2).MarkerEdgeColor = [0.85 0.25 0.25];
        scatters(2).SizeData = 28;
        scatters(2).LineWidth = 0.8;
    end
end

% 强制所有子图的 Position left/宽度/高度对齐
ref_pos = ax_dist(1).Position;
for j = 2:n_ad
    cur = ax_dist(j).Position;
    ax_dist(j).Position = [ref_pos(1), cur(2), ref_pos(3), ref_pos(4)];
end

% 添加 (a)(b)(c) 标签 (统一 x 位置)
abc = {'(a)', '(b)', '(c)'};
label_x = ax_dist(1).OuterPosition(1) + 0.005;
for j = 1:n_ad
    op = ax_dist(j).OuterPosition;
    annotation(fig_dist, 'textbox', ...
        [label_x, op(2)+op(4)-0.02, 0.06, 0.04], ...
        'String', abc{j}, 'FontSize', 15, 'FontWeight', 'bold', ...
        'EdgeColor', 'none', 'VerticalAlignment', 'top');
end

% 图例移到第一个子图上方 (远离边缘)
lg_list = findobj(fig_dist, 'Type', 'legend');
for k = 1:length(lg_list)
    lg_list(k).FontSize = 11;
    lg_list(k).Box = 'on';
    lg_list(k).Position = [0.72, 0.93, 0.20, 0.05];  % 上移到 y=0.93-0.98
end

export_one(fig_dist, fullfile(target_dir, '系统级畸变对比'), 7.0, font_name);
close(fig_dist);

fprintf('\n===== 全部完成 =====\n');
fprintf('文件在 %s/\n', target_dir);

%% ========================================================================
%%  辅助函数
%% ========================================================================

function [fig, data] = load_and_clean_density(fig_name)
    % 打开密度图, 在内存中删除 spatial 子图, 保留 alpha + U0
    fig = openfig(fig_name, 'invisible');
    fig.Visible = 'off';

    ax_all = findobj(fig, 'Type', 'axes');
    for j = 1:length(ax_all)
        t = ax_all(j).Title.String;
        if contains(t, 'spatial') || contains(t, 'Initial')
            delete(ax_all(j));
        end
    end

    % 提取剩余的 alpha 和 U0 子图
    ax_remain = findobj(fig, 'Type', 'axes');
    data = struct('alpha', [], 'u0', []);
    for j = 1:length(ax_remain)
        t = ax_remain(j).Title.String;
        if contains(t, 'alpha') || contains(t, '\alpha')
            data.alpha = ax_remain(j);
        elseif contains(t, 'U_0') || contains(t, 'U0')
            data.u0 = ax_remain(j);
        end
    end
end

function copy_axes_content(ax_dest, ax_src)
    if isempty(ax_src) || ~isvalid(ax_src)
        return;
    end
    children = ax_src.Children;
    copyobj(children, ax_dest);

    % 复制关键属性
    ax_dest.XLabel.String = ax_src.XLabel.String;
    ax_dest.YLabel.String = ax_src.YLabel.String;
    ax_dest.XScale = ax_src.XScale;
    ax_dest.YScale = ax_src.YScale;
    try; ax_dest.XLim = ax_src.XLim; catch; end
    try; ax_dest.YLim = ax_src.YLim; catch; end
    ax_dest.XGrid = ax_src.XGrid;
    ax_dest.YGrid = ax_src.YGrid;
end

function export_one(fig, path_no_ext, width_inch, font_name)
    % 统一字体
    set(findall(fig, '-property', 'FontName'), 'FontName', font_name);

    pos = fig.Position;
    aspect = pos(4) / pos(3);
    height = width_inch * aspect;

    % 加 10% 边距防截断
    margin_inch = 0.3;
    paper_w = width_inch + 2*margin_inch;
    paper_h = height + 2*margin_inch;

    set(fig, 'Units', 'inches', ...
             'Position', [1, 1, width_inch, height], ...
             'PaperUnits', 'inches', ...
             'PaperSize', [paper_w, paper_h], ...
             'PaperPosition', [margin_inch, margin_inch, width_inch, height]);

    fprintf('  Figure: %.1f x %.1f in\n', width_inch, height);

    png_path = [path_no_ext '.png'];
    print(fig, png_path, '-dpng', '-r600', '-opengl');
    fprintf('  -> %s\n', png_path);

    tif_path = [path_no_ext '.tif'];
    print(fig, tif_path, '-dtiff', '-r600', '-opengl');
    fprintf('  -> %s\n', tif_path);

    eps_path = [path_no_ext '.eps'];
    try
        print(fig, eps_path, '-depsc', '-painters', '-r600');
        fprintf('  -> %s\n', eps_path);
    catch
        fprintf('  [跳过] EPS导出失败, 请用TIFF\n');
    end
end
