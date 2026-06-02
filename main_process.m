% main_init_test.m

clear; clc; close all;

%% === 0. 开启全局计时器 ===
total_timer = tic;
fprintf('==================================================\n');
fprintf('开始运行电子束光学转移矩阵仿真...\n');
fprintf('==================================================\n');

%% 粒子初始化
% 设定模拟数量
num_electrons = 50000;

% 配置初始化参数
config.load_from_file = true; % 设为 true 时尝试读取文件
config.file_path = 'existing_particles.mat';
config.pattern_mode = "spot";

% 条纹几何参数 (假设单位为米 m，这里设置微米级别演示)
config.stripe_center_x = 1e-3;
config.stripe_center_y = 0;
config.stripe_width = 36e-6;    % 线宽 2 um
config.stripe_length = 100e-6;  % 线长 50 um
config.stripe_spacing = 72e-6;  % 周期 6 um

% --- 圆斑发射参数设置 ---
config.spot_radius = 18e-9;  % 设定圆斑半径，例如 500 nm
config.spot_center_x = 0;     % 圆心 X 坐标
config.spot_center_y = 0;     % 圆心 Y 坐标

% 物理分布参数
config.U0_max = 0.5;           % 初始能量上限 (eV/V)
config.alpha_max_deg = 10.7;   % 最大发射角 (度)

% 调用初始化函数
[V_init, p] = init_particles(num_electrons, config);

% --- 数据验证与可视化 ---
x_plot = squeeze(V_init(1,1,:)) * 1e9; % 转为微米
y_plot = squeeze(V_init(2,1,:)) * 1e9;

figure('Name', '电子初始分布');

% 1. 空间分布散点图
subplot(1,3,1);
scatter(x_plot, y_plot, 1, 'b', '.');
axis equal;
title('Initial spatial distribution');
xlabel('x (nm)'); ylabel('y (nm)');

% 2. 能量 U0 直方图
subplot(1,3,2);
histogram(p.U0, 50, 'Normalization', 'pdf');
title('U_0 distribution (0-0.5V)');
xlabel('U_0 (V)'); ylabel('probability density');

% 3. 发射角 Alpha 直方图
subplot(1,3,3);
histogram(rad2deg(p.alpha), 50, 'Normalization', 'pdf');
title('\alpha distribution (\leq 10.7^{\circ})');
xlabel('\alpha (deg)'); ylabel('probability density');

%% 矩阵构造

% 1. 定义系统与误差参数 (可根据需要提取到 config 结构体中)
% --- 方案 A: 使用解析函数定义场分布 ---
% 例如：假设磁场中心为 B0，边缘有抛物线型的衰减
sys.B0 = 1.2;
sys.L = 1e-3;    % 极间距
sys.U = 12.67e3;  % 轴向总加速电压 (用于能量计算)
err.G = 6*0;

sys.func_B_vec = @(X, Y, Z) [
    -0.5 * err.G .* X;             % Bx (由 -0.5*G*r 分解得到)
    -0.5 * err.G .* Y ;             % By
    sys.B0 + err.G .* ( Z - 0.5 * sys.L ) .* ones(size(X))            % Bz (主磁场)
    ];

% 假设电场在 x 方向有微小的偏移引起不平行
sys.func_E_vec = @(X, Y) [
    (sys.U / sys.L) * ones(size(X))*1e-3*0.8*0;           % Ex 误差
    (sys.U / sys.L) * ones(size(X))*1e-3*0.6*0;                % Ey
    (sys.U / sys.L) * ones(size(X))% Ez
    ];

sys.N_turn = 1;      % 聚焦圈数 N

err.delta_E = 1.267e4*0;   % 电场0阶误差 (V/m)
err.delta_B = 1.2e-3*0;  % 磁场0阶误差 (T)

% 设定局部电场畸变发生的 Z 坐标区间 (例如在掩模附近 0.01m 到 0.015m 处)
err.z1 = 0;
err.z2 = 1e-3;

% 定义局部电场误差函数 \delta Ez(x, y, z)
% 假设由于掩模孔径效应，中心电场无误差，越靠近边缘误差越大，
% 并且在 z 方向上呈现抛物线或高斯分布的衰减。
err.func_dEz = @(x, y, z) ...
    -6.335e7*z*0;
% 5000 是强度系数，这只是一个模拟边缘场增强的数学假设函数

% 在主脚本参数定义部分
err.I_beam = 3e-15;   % 束流大小 1uA
err.R0 = 18e-9;       % 初始光斑半径 1um
err.ee_method = 'gauss';

% 2. 预计算所有电子共用的物理量 (完全向量化)
% 提前计算好时间 t、相位 beta、半径 R 等，避免在各个矩阵函数里重复计算
phys = calc_common_physics(p, V_init, sys);

fprintf('正在构建转移矩阵...\n');

% 3. 调用矩阵构建接口 (每个函数返回 3x3xN 的数组)
% 注意：这里的子函数我们暂时留空，只定义输入输出接口
M_spread = calc_M_spread(p, phys, num_electrons);           % 能散与发散角矩阵 [cite: 119]
M_EB0    = calc_M_EB0(p, phys, sys, err, num_electrons);    % 电磁场0阶误差矩阵 [cite: 140]
M_Br1    = calc_M_Br1(p, phys, sys, err, num_electrons);    % 磁场r方向1阶误差(旋转畸变) [cite: 175]
M_EBrho  = calc_M_EB_angle(p, phys, sys, err, num_electrons); % 电磁场夹角误差(拉伸畸变) [cite: 183]
M_ee     = calc_M_ee(p, phys, sys, err, num_electrons);     % 空间电荷效应(电子间排斥) [cite: 200]
M_E1     = calc_M_E1(p, phys, sys, err, num_electrons);     % 电场Ez的z方向1阶及高阶误差
M_B1     = calc_M_B1(p, phys, sys, err, num_electrons);     % 磁场Bz的z方向1阶及高阶误差

% 4. 矩阵组合 (使用 pagemtimes 进行 3D 数组的高效页乘)
fprintf('正在进行批量矩阵运算...\n');
% 乘法顺序极度重要：靠近初始状态的在右边，最后发生的在左边
% 假设叠加顺序为: 初始 -> 发散 -> 0阶误差 -> 旋转畸变 -> 角度误差 -> 空间电荷
M_temp1 = pagemtimes(M_EB0, M_spread);
M_temp2 = pagemtimes(M_Br1, M_temp1);
M_temp3 = pagemtimes(M_E1, M_temp2);
M_temp4 = pagemtimes(M_B1, M_temp3);
M_temp5 = pagemtimes(M_EBrho, M_temp4);
M_total = pagemtimes(M_ee, M_temp5);

% 5. 计算最终位置向量
% M_total 是 3x3xN, V_init 是 3x1xN，结果 V_final 是 3x1xN
V_final = pagemtimes(M_total, V_init);

% 6. 数据分离与提取
% squeeze 会把大小为 1 的维度压缩掉，使得 3x1xN 变成 3xN 的标准二维矩阵
Pos_final = squeeze(V_final);
X_final = Pos_final(1, :);
Y_final = Pos_final(2, :);

fprintf('位置计算完成！\n');

%% 7. 初始与最终电子分布对照图可视化
fprintf('正在绘制初始与最终电子分布对照图...\n');

% 将坐标单位转换为纳米 (nm) 以方便对比
scale_factor = 1e9;
X_init_nm = squeeze(V_init(1, 1, :)) * scale_factor;
Y_init_nm = squeeze(V_init(2, 1, :)) * scale_factor;
X_final_nm = X_final * scale_factor;
Y_final_nm = Y_final * scale_factor;

% 创建新的对比窗口
figure('Name', '电子束传输前后分布对照', 'Position', [100, 100, 1000, 450]);

% --- 左图：初始分布 ---
subplot(1, 2, 1);
scatter(X_init_nm, Y_init_nm, 1, 'b', '.');
axis equal; % 保持 X 和 Y 比例一致
grid on;
title('初始电子空间分布 (光电阴极)');
xlabel('X (nm)');
ylabel('Y (nm)');

% --- 右图：最终分布 ---
subplot(1, 2, 2);
scatter(X_final_nm, Y_final_nm, 1, 'r', '.');
axis equal; % 保持 X 和 Y 比例一致
grid on;
title('最终电子空间分布 (硅片靶面)');
xlabel('X (nm)');
ylabel('Y (nm)');

% --- 计算并打印定量的模糊度/扩展量 ---
RMS_blur_X = std(X_final_nm - X_init_nm');
RMS_blur_Y = std(Y_final_nm - Y_init_nm');
fprintf('----------------------------------------\n');
fprintf('传输完成！靶面处 X 方向的 RMS 模糊度: %.2f nm\n', RMS_blur_X);
fprintf('传输完成！靶面处 Y 方向的 RMS 模糊度: %.2f nm\n', RMS_blur_Y);
fprintf('----------------------------------------\n');

%% --- 8. 单图初始与最终分布对比连线图 ---
% 在同一张图上绘制初始点(空心圆)、最终点(实心点)并用直线连接

% 创建新窗口
figure('Name', '电子初始与最终位置对比');
hold on;
grid on;

% 绘制连线 (用稍淡的颜色和细线)
% 向量化绘制连线比较快的方法是用 plot 绘制矩阵，但需要调整矩阵形状
% 这里用一个循环，虽然有循环但只有线，绘图引擎处理很快，或者用线段矩阵
% 为了最高效绘图，我们构建交叉的矩阵
line_x = [X_init_nm'; X_final_nm];
line_y = [Y_init_nm'; Y_final_nm];

% 绘制初始点 (蓝色空心圆)
scatter(X_init_nm, Y_init_nm, 10, 'b', 'o');

% 绘制最终点 (红色实心小点)
scatter(X_final_nm, Y_final_nm, 15, 'r', 'filled');

% 坐标轴设置
axis equal;
title('电子空间分布变化追踪');
xlabel('X (nm)');
ylabel('Y (nm)');
legend('初始位置 (蓝圆)', '最终位置 (红点)', 'Location', 'best');
hold off;

%% === 8. 结束全局计时并输出运行效率 ===
total_time_seconds = toc(total_timer);

fprintf('\n==================================================\n');
fprintf('               仿真计算与可视化全部完成！             \n');
fprintf('--------------------------------------------------\n');
fprintf('总计模拟粒子数 : %d 个\n', num_electrons);
fprintf('程序运行总耗时 : %.4f 秒\n', total_time_seconds);
fprintf('平均单粒子耗时 : %.2e 秒/个\n', total_time_seconds / num_electrons);
fprintf('==================================================\n');

%% 9.圆斑参数分析
if config.pattern_mode == "spot"
    % === 靶面圆斑分辨率分析与量化 ===

    % 1. 计算电子束的质心 (Beam Centroid)
    cx = mean(X_final);
    cy = mean(Y_final);

    % 2. 计算每个电子到质心的径向距离 (Radial Distance)
    r_final = sqrt((X_final - cx).^2 + (Y_final - cy).^2);

    % 3. 计算 RMS 半径
    RMS_radius = sqrt(mean(r_final.^2));

    % 4. 计算 FW50 和 FW99 (对距离进行排序)
    r_sorted = sort(r_final);
    N_total = length(r_sorted);

    % 寻找包含 50% 和 99% 电子的临界半径
    r_50 = r_sorted(round(0.50 * N_total));
    r_90 = r_sorted(round(0.90 * N_total));
    r_99 = r_sorted(round(0.99 * N_total));

    % Full Width (直径) 是半径的两倍
    FW50 = 2 * r_50;
    FW90 = 2 * r_90;
    FW99 = 2 * r_99;

    % --- 打印输出量化结果 ---
    fprintf('\n==================================================\n');
    fprintf('                MATLAB靶面光斑分辨率量化指标               \n');
    fprintf('--------------------------------------------------\n');
    fprintf('光斑质心偏移 (cx, cy): (%.2f nm, %.2f nm)\n', (cx -config.spot_center_x) * 1e9, (cy-config.spot_center_y) * 1e9);
    fprintf('RMS 光斑半径 (r_rms) : %.2f nm\n', RMS_radius * 1e9);
    fprintf('FW50 (包含50%%电子直径): %.2f nm\n', FW50 * 1e9);
    fprintf('FW90 (包含90%%电子直径): %.2f nm\n', FW90 * 1e9);
    fprintf('FW99 (包含99%%电子直径): %.2f nm\n', FW99 * 1e9);
    fprintf('==================================================\n');

    % --- 绘制专业的束斑量化分析图 ---
    figure('Name', 'MATLAB电子束斑分辨率分析', 'Position', [100, 100, 1000, 450]);

    % 子图 1: 最终散点图及包络圆圈
    subplot(1, 2, 1);
    hold on; grid on; axis equal;
    % 画出电子散点
    scatter(X_final * 1e9, Y_final * 1e9, 2, 'b', '.');

    % 绘制质心
    plot(cx * 1e9, cy * 1e9, 'k+', 'MarkerSize', 10, 'LineWidth', 2);

    % 绘制 FW50 和 FW99 的边界圆
    theta = linspace(0, 2*pi, 100);
    plot((cx + r_50 * cos(theta)) * 1e9, (cy + r_50 * sin(theta)) * 1e9, 'g-', 'LineWidth', 2);
    plot((cx + r_90 * cos(theta)) * 1e9, (cy + r_90 * sin(theta)) * 1e9, 'y--', 'LineWidth', 2);
    plot((cx + r_99 * cos(theta)) * 1e9, (cy + r_99 * sin(theta)) * 1e9, 'r--', 'LineWidth', 2);

    title('靶面光斑分布与 FW50/FW99 包络');
    xlabel('X (nm)'); ylabel('Y (nm)');
    legend('电子落点', '束流质心', 'FW50 边界', 'FW90 边界', 'FW99 边界', 'Location', 'best');

    % 子图 2: 径向电子累积分布函数 (CDF)
    subplot(1, 2, 2);
    cumulative_fraction = (1:N_total) / N_total;
    plot(r_sorted * 1e9, cumulative_fraction * 100, 'b-', 'LineWidth', 2);
    hold on; grid on;
    % 标出 FW50 和 FW99 的截断点
    yline(50, 'g--', 'LineWidth', 1.5);
    yline(90, 'b--', 'LineWidth', 1.5);
    yline(99, 'r--', 'LineWidth', 1.5);
    xline(r_50 * 1e9, 'g--', 'LineWidth', 1.5);
    xline(r_90 * 1e9, 'b--', 'LineWidth', 1.5);
    xline(r_99 * 1e9, 'r--', 'LineWidth', 1.5);

    title('径向电子流强度累积分布');
    xlabel('距质心径向距离 r (nm)');
    ylabel('包含的电子比例 (%)');
    ylim([0 100]);

    if config.load_from_file == true

        % 1. 读取 COMSOL 导出的 CSV 文件
        fprintf('正在读取 COMSOL 数据...\n');
        filename = '06_b_field_gradient_edge_image.csv';
        % readmatrix 能够自动忽略以 '%' 开头的 COMSOL 注释行
        data = readmatrix(filename);

        % data 的列对应关系：
        % 第 1 列: qx (um)
        % 第 2 列: qy (um)

        N = size(data, 1); % 电子总数

        % 2. 提取并转换空间坐标 (转换为米 m)
        % COMSOL 导出的单位是微米 (um)，需乘以 1e-6 统一到国际标准单位
        x = -data(:, 1)' * 1e-6; % 转为 1xN 行向量
        y = data(:, 2)' * 1e-6;

        % 构造 3x1xN 的齐次坐标矩阵 V_init
        X_c_result(1, :) = x;
        Y_c_result(1, :) = y;

        % === 靶面圆斑分辨率分析与量化 ===

        % 1. 计算电子束的质心 (Beam Centroid)
        cx_c_result = mean(X_c_result);
        cy_c_result = mean(Y_c_result);

        % 2. 计算每个电子到质心的径向距离 (Radial Distance)
        r_c_result = sqrt((X_c_result - cx_c_result).^2 + (Y_c_result - cy_c_result).^2);

        % 3. 计算 RMS 半径
        RMS_radius_c_result = sqrt(mean(r_c_result.^2));

        % 4. 计算 FW50 和 FW99 (对距离进行排序)
        r_sorted_c_result = sort(r_c_result);
        N_total_c_result = length(r_sorted_c_result);

        % 寻找包含 50% 和 99% 电子的临界半径
        r_50_c_result = r_sorted_c_result(round(0.50 * N_total_c_result));
        r_90_c_result = r_sorted_c_result(round(0.90 * N_total_c_result));
        r_99_c_result = r_sorted_c_result(round(0.99 * N_total_c_result));

        % Full Width (直径) 是半径的两倍
        FW50_c_result = 2 * r_50_c_result;
        FW90_c_result = 2 * r_90_c_result;
        FW99_c_result = 2 * r_99_c_result;

        % --- 打印输出量化结果 ---
        fprintf('\n==================================================\n');
        fprintf('                COMSOL靶面光斑分辨率量化指标               \n');
        fprintf('--------------------------------------------------\n');
        fprintf('光斑质心偏移 (cx, cy): (%.2f nm, %.2f nm)\n', cx_c_result * 1e9, cy_c_result * 1e9);
        fprintf('RMS 光斑半径 (r_rms) : %.2f nm\n', RMS_radius_c_result * 1e9);
        fprintf('FW50 (包含50%%电子直径): %.2f nm\n', FW50_c_result * 1e9);
        fprintf('FW90 (包含90%%电子直径): %.2f nm\n', FW90_c_result * 1e9);
        fprintf('FW99 (包含99%%电子直径): %.2f nm\n', FW99_c_result * 1e9);
        fprintf('==================================================\n');

        % --- 绘制专业的束斑量化分析图 ---
        figure('Name', 'COMSOL电子束斑分辨率分析', 'Position', [100, 100, 1000, 450]);

        % 子图 1: 最终散点图及包络圆圈
        subplot(1, 2, 1);
        hold on; grid on; axis equal;
        % 画出电子散点
        scatter(X_c_result * 1e9, Y_c_result * 1e9, 2, 'b', '.');

        % 绘制质心
        plot(cx_c_result * 1e9, cy_c_result * 1e9, 'k+', 'MarkerSize', 10, 'LineWidth', 2);

        % 绘制 FW50 和 FW90 的边界圆
        theta = linspace(0, 2*pi, 100);
        plot((cx_c_result + r_50_c_result * cos(theta)) * 1e9, (cy_c_result + r_50_c_result * sin(theta)) * 1e9, 'g-', 'LineWidth', 2);
        plot((cx_c_result + r_90_c_result * cos(theta)) * 1e9, (cy_c_result + r_90_c_result * sin(theta)) * 1e9, 'y--', 'LineWidth', 2);
        plot((cx_c_result + r_99_c_result * cos(theta)) * 1e9, (cy_c_result + r_99_c_result * sin(theta)) * 1e9, 'r--', 'LineWidth', 2);

        title('靶面光斑分布与 FW50/FW99 包络');
        xlabel('X (nm)'); ylabel('Y (nm)');
        legend('电子落点', '束流质心', 'FW50 边界', 'FW90 边界', 'FW99 边界', 'Location', 'best');

        % 子图 2: 径向电子累积分布函数 (CDF)
        subplot(1, 2, 2);
        cumulative_fraction_c_result = (1:N_total_c_result) / N_total_c_result;
        plot(r_sorted_c_result * 1e9, cumulative_fraction_c_result * 100, 'b-', 'LineWidth', 2);
        hold on; grid on;
        % 标出 FW50 和 FW90 的截断点
        yline(50, 'g--', 'LineWidth', 1.5);
        yline(90, 'b--', 'LineWidth', 1.5);
        yline(99, 'r--', 'LineWidth', 1.5);
        xline(r_50_c_result * 1e9, 'g--', 'LineWidth', 1.5);
        xline(r_90_c_result * 1e9, 'b--', 'LineWidth', 1.5);
        xline(r_99_c_result * 1e9, 'r--', 'LineWidth', 1.5);

        title('径向电子流强度累积分布');
        xlabel('距质心径向距离 r (nm)');
        ylabel('包含的电子比例 (%)');
        ylim([0 100]);

        % 创建新的对比窗口
        figure('Name', '电子束传输前后分布对照', 'Position', [100, 100, 1000, 450]);

        % --- 左图：初始分布 ---
        subplot(1, 3, 1);
        scatter(X_init_nm, Y_init_nm, 1, 'b', '.');
        axis equal; % 保持 X 和 Y 比例一致
        grid on;
        title('Initial electron distribution');
        xlabel('x (nm)');
        ylabel('y (nm)');

        % --- 中图：MATLAB最终分布 ---
        subplot(1, 3, 2);
        scatter(X_final_nm, Y_final_nm, 1, 'r', '.');
        axis equal; % 保持 X 和 Y 比例一致
        grid on;
        title('Final electron distribution (TMM)');
        xlabel('x (nm)');
        ylabel('y (nm)');
        hold on;
        theta = linspace(0, 2*pi, 100);
        plot((cx + r_50 * cos(theta)) * 1e9, (cy + r_50 * sin(theta)) * 1e9, 'g-', 'LineWidth', 2);
        plot((cx + r_90 * cos(theta)) * 1e9, (cy + r_90 * sin(theta)) * 1e9, 'y--', 'LineWidth', 2);
        plot((cx + r_99 * cos(theta)) * 1e9, (cy + r_99 * sin(theta)) * 1e9, 'b.', 'LineWidth', 2);
        legend('electron landing point', 'FW50', 'FW90', 'FW99', 'Location', 'best');

        % --- 右图：COMSOL最终分布 ---
        subplot(1, 3, 3);
        scatter(X_c_result * 1e9, Y_c_result * 1e9, 1, 'r', '.');
        axis equal; % 保持 X 和 Y 比例一致
        grid on;
        title('Final electron distribution (FEM)');
        xlabel('x (nm)');
        ylabel('y (nm)');
        hold on;
        theta = linspace(0, 2*pi, 100);
        plot((cx_c_result + r_50_c_result * cos(theta)) * 1e9, (cy_c_result + r_50_c_result * sin(theta)) * 1e9, 'g-', 'LineWidth', 2);
        plot((cx_c_result + r_90_c_result * cos(theta)) * 1e9, (cy_c_result + r_90_c_result * sin(theta)) * 1e9, 'y--', 'LineWidth', 2);
        plot((cx_c_result + r_99_c_result * cos(theta)) * 1e9, (cy_c_result + r_99_c_result * sin(theta)) * 1e9, 'b.', 'LineWidth', 2);
        legend('electron landing point', 'FW50', 'FW90', 'FW99', 'Location', 'best');

        % --- 可选：统一坐标轴范围以直观体现缩放和畸变 ---
        % 计算包含所有初始和最终点的最大范围
        max_val = max([abs(X_init_nm'), abs(Y_init_nm'), abs(X_final_nm), abs(Y_final_nm)]);
        axis_limit = max_val * 1.1; % 留出 10% 的空白边缘
    end
else
    if config.pattern_mode == "stripe"
        if config.load_from_file == true
            % 1. 读取 COMSOL 导出的 CSV 文件
            fprintf('正在读取 COMSOL 数据...\n');
            filename = '12_system_stripe_1mm_image.csv';
            % readmatrix 能够自动忽略以 '%' 开头的 COMSOL 注释行
            data = readmatrix(filename);

            % data 的列对应关系：
            % 第 1 列: qx (um)
            % 第 2 列: qy (um)

            N = size(data, 1); % 电子总数

            % 2. 提取并转换空间坐标 (转换为米 m)
            % COMSOL 导出的单位是微米 (um)，需乘以 1e-6 统一到国际标准单位
            x = data(:, 1)' * 1e-6; % 转为 1xN 行向量
            y = data(:, 2)' * 1e-6;

            % 构造 3x1xN 的齐次坐标矩阵 V_init
            X_c_result(1, :) = x * 1e9;
            Y_c_result(1, :) = y * 1e9;

            % === 改进后的条纹畸变对比绘图部分 ===
            
            % 1. 将画板拉长变扁，适应条纹分布 (X宽Y窄) [X位置, Y位置, 宽度, 高度]
            figure('Name', 'MATLAB与COMSOL计算条纹对比', 'Position', [100, 200, 1500, 400]);
            
            % 提取 TMM 数据并进行中心平移修正
            x_tmm_plot = X_final_nm - config.stripe_center_x * 1e9;
            y_tmm_plot = Y_final_nm - config.stripe_center_y * 1e9;
            
            % 2. 计算所有数据的整体边界，并加上 5% 的留白裕量
            x_min = min([x_tmm_plot, X_c_result]);
            x_max = max([x_tmm_plot, X_c_result]);
            y_min = min([y_tmm_plot, Y_c_result]);
            y_max = max([y_tmm_plot, Y_c_result]);
            
            x_pad = (x_max - x_min) * 0.05;
            y_pad = (y_max - y_min) * 0.05;
            
            % 设置统一的坐标轴范围
            plot_xlim = [x_min - x_pad, x_max + x_pad];
            plot_ylim = [y_min - y_pad, y_max + y_pad];

            % --- 左图：TMM 结果 ---
            subplot(1,3,1)
            scatter(x_tmm_plot, y_tmm_plot, 10, 'b', 'o');
            axis equal;
            xlim(plot_xlim); ylim(plot_ylim); % 强制收缩白边
            title('TMM results');
            xlabel('x (nm)');
            ylabel('y (nm)');
            grid on;

            % --- 中图：FEM (COMSOL) 结果 ---
            subplot(1,3,2)
            scatter(X_c_result, Y_c_result, 15, 'r', 'filled');
            axis equal;
            xlim(plot_xlim); ylim(plot_ylim); % 强制收缩白边
            title('FEM results');
            xlabel('x (nm)');
            ylabel('y (nm)');
            grid on;

            % --- 右图：对比图 ---
            subplot(1,3,3)
            hold on;
            % 绘制初始点 (蓝色空心圆)
            scatter(x_tmm_plot, y_tmm_plot, 10, 'b', 'o');
            % 绘制最终点 (红色实心小点)
            scatter(X_c_result, Y_c_result, 15, 'r', 'filled');
            legend('TMM results', 'FEM results', 'Location', 'best');
            axis equal;
            xlim(plot_xlim); ylim(plot_ylim); % 强制收缩白边
            title('system-level distortion comparison');
            xlabel('x (nm)');
            ylabel('y (nm)');
            grid on;
            hold off;
        end
    end
end