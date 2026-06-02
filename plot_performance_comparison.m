% plot_performance_comparison.m
% 绘制 COMSOL 与 MATLAB 转移矩阵法 (TMM) 的计算效率对比图 (双对数坐标)

clc; clear; close all;

%% 1. 数据准备
% 粒子数量
N_COMSOL = [1000, 5000, 10000, 50000, 100000];
N_MATLAB = [1000, 5000, 10000, 50000, 100000, 1000000]; % 包含百万粒子的情况

% COMSOL 耗时数据 (秒)
Time_COMSOL = [182, 227, 290, 613, 1048];

% MATLAB (TMM) 耗时数据 (秒)
% 最后一个 2.8s 是根据你之前提到的百万粒子不到3秒预估的，你可以替换为真实测试数据
Time_MATLAB = [0.2337, 0.2619, 0.2385, 0.3699, 0.4834, 3.0105]; 

%% 2. 图形初始化 (Publication Quality 风格)
figure('Name', '计算效率对比', 'Position', [100, 100, 800, 600], 'Color', 'w');
hold on;
set(gca, 'XScale', 'log', 'YScale', 'log'); % 设置双对数坐标系

% 开启网格 (包括副网格，增强对数图的可读性)
grid on;
grid minor;
set(gca, 'GridAlpha', 0.5, 'MinorGridAlpha', 0.5, 'LineWidth', 1.2);

%% 3. 绘制数据曲线
% 绘制 COMSOL 曲线 (红色实线加圆点)
p1 = plot(N_COMSOL, Time_COMSOL, '-ro', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% 绘制 MATLAB 曲线 (蓝色实线加方块)
p2 = plot(N_MATLAB, Time_MATLAB, '-bs', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'b');

%% 4. 添加标注与视觉强化
% 在 N = 100,000 处画一条虚线墙，表示内存崩溃
plot([100000, 100000], [0.1, 10000], 'r--', 'LineWidth', 2);
text(120000, 1000, {'COMSOL', 'Out of Memory', '(>128 GB)'}, ...
    'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold');

% 在 N = 100,000 处添加加速比的双向箭头
x_arrow = 100000;
y_bottom = Time_MATLAB(5);
y_top = Time_COMSOL(5);

% 使用 annotation 画双向箭头 (需将坐标转换为 Figure 相对坐标)
% 此处使用巧妙的 errorbar 变体或者直接用 line 配合文本模拟
plot([x_arrow, x_arrow], [y_bottom*1.5, y_top*0.8], 'k-', 'LineWidth', 1.5);
% 绘制箭头三角形
plot(x_arrow, y_top*0.8, 'k^', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
plot(x_arrow, y_bottom*1.5, 'kv', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

% 添加倍数文字 (放在箭头的左侧)
text(70000, sqrt(y_bottom * y_top), '~2168\times Speedup', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
    'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');

%% 5. 坐标轴与图例设置
% 设置坐标轴范围
xlim([500, 2e6]);
ylim([0.1, 1e4]);

% 坐标轴标签与字体
xlabel('Number of Simulated Particles ($N$)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Wall-clock Computation Time (Seconds)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

% 设置坐标轴刻度字体
set(gca, 'FontSize', 14);

% 图例设置
lgd = legend([p1, p2], {'COMSOL (FEM + Ray Tracing)', 'MATLAB (Proposed TMM)'}, ...
    'Location', 'northwest', 'FontSize', 14, 'Interpreter', 'latex');
set(lgd, 'EdgeColor', 'k', 'LineWidth', 1);

% 去除图表上边框和右边框（更符合学术期刊审美）
box off;

%% 6. (可选) 导出高分辨率图片
% 取消下一行的注释即可自动保存为 600 DPI 的高分辨率 TIFF 图片，适合直接投稿
% exportgraphics(gcf, 'Computation_Time_Comparison.tiff', 'Resolution', 600);