% main_optimize_U.m
clear; clc; close all;

%% === 1. 初始化设置与粒子生成 (仅执行一次) ===
total_timer = tic;
fprintf('==================================================\n');
fprintf('开始运行电子束光学系统优化 (寻找最佳加速电压 U)...\n');
fprintf('==================================================\n');

% 设定模拟数量 (建议在优化时适当降低数量以节省时间，比如10万)
num_electrons = 100000; 

% 配置初始化参数 (仅保留圆斑)
config.load_from_file = false; 
config.spot_radius = 18e-9;   % 设定圆斑半径，18 nm
config.spot_center_x = 0;     
config.spot_center_y = 0;     

% 物理分布参数
config.U0_max = 1.2;           % 初始能量上限 (eV/V)
config.alpha_max_deg = 17.57;   % 最大发射角 (度)

fprintf('正在生成 %d 个初始粒子...\n', num_electrons);
% 调用初始化函数 (全局仅调用一次，消除优化过程中的随机噪声)
[V_init, p] = init_particles(num_electrons, config);

% --- 数据验证与可视化 ---
x_plot = squeeze(V_init(1,1,:)) * 1e9; % 转为纳米 (nm)
y_plot = squeeze(V_init(2,1,:)) * 1e9;

% 创建窗口并设置一个合适的比例 (宽1200，高400) 适合1x3排列
figure('Name', '电子初始分布', 'Position', [100, 200, 1200, 400]);

% 1. 空间分布散点图
subplot(1,3,1);
scatter(x_plot, y_plot, 1, 'b', '.');
axis equal;
% 动态提取 config 中的半径参数并转换为 nm
title(sprintf('初始空间分布 (半径 r = %.1f nm)', config.spot_radius * 1e9));
xlabel('X (nm)'); ylabel('Y (nm)');

% 2. 能量 U0 直方图
subplot(1,3,2);
histogram(p.U0, 50, 'Normalization', 'pdf');
% 动态提取 config 中的能量上限
title(sprintf('U_0 能量分布 (0 - %.2f V)', config.U0_max));
xlabel('U_0 (V)'); ylabel('概率密度');

% 3. 发射角 Alpha 直方图
subplot(1,3,3);
histogram(rad2deg(p.alpha), 50, 'Normalization', 'pdf');
% 动态提取 config 中的最大角度
title(sprintf('Alpha 发射角分布 (\\leq %.1f^{\\circ})', config.alpha_max_deg));
xlabel('\alpha (deg)'); ylabel('概率密度');

% 可选：为整个窗口添加一个总标题 (适用于 MATLAB R2018b 及更高版本)
% sgtitle(sprintf('电子束初始状态验证 (总粒子数: %d)', num_electrons));

%% === 2. 系统固定参数与误差定义 ===
% 固定系统参数
sys.B0 = 1.4;         % 磁场中心 B0 (T)
sys.L = 1.2e-3;       % 极间距 (m)
sys.N_turn = 1;       % 聚焦圈数 N

% 误差参数
err.G = 0;
err.delta_E = 0;             % 电场0阶误差 (V/m)
err.delta_B = 0;   % 磁场0阶误差 (T)
err.I_beam = 3e-15;          % 束流大小 1uA (这里你的注释是1uA，但数值是3e-15A，保留原值)
err.R0 = 5e-9;              % 初始光斑半径
err.ee_method = 'gauss';
err.z1 = 0.010;
err.z2 = 0.015;
err.func_dEz = @(x, y, z) (5000 * (x.^2 + y.^2)) .* exp(-((z - 0.0125)/0.002).^2) * 0;

% 定义磁场函数 (不依赖于U，可以固定)
sys.func_B_vec = @(X, Y, Z) [
    -0.5 * err.G .* X;             
    -0.5 * err.G .* Y;             
    sys.B0 + err.G .* ( Z - 0.5 * sys.L ) .* ones(size(X))            
];

%% === 3. 运行优化算法寻找最佳 U ===
% 设定电压 U 的搜索范围 (你需要根据物理实际情况设定上下限)
% 例如，原设为 24.927735e3，我们可以在 20000 到 30000 之间寻找
U_min = 20e3;
U_max = 30e3;

fprintf('\n开始在 [%.0f V, %.0f V] 范围内优化加速电压 U...\n', U_min, U_max);

% 设置优化器选项：显示每次迭代的信息
options = optimset('Display', 'iter', 'TolX', 1e-3);

% 定义目标函数：输入为 U_test，输出为 FW99
objective_func = @(U_test) eval_FW99(U_test, V_init, p, sys, err, num_electrons);

% 使用 fminbnd 进行一维有界优化
[opt_U, min_FW99] = fminbnd(objective_func, U_min, U_max, options);

fprintf('\n==================================================\n');
fprintf('优化完成！\n');
fprintf('最佳加速电压 U = %.6f V\n', opt_U);
fprintf('最小 FW99 = %.2f nm\n', min_FW99 * 1e9);
fprintf('==================================================\n');

%% === 4. 使用最佳 U 进行最终模拟与量化 ===
fprintf('\n正在使用最佳电压 U = %.2f V 生成最终分布图...\n', opt_U);

% 更新 sys 结构体为最佳电压
sys.U = opt_U;
sys.func_E_vec = @(X, Y) [
    zeros(size(X));           
    zeros(size(X));                
    (sys.U / sys.L) * ones(size(X))
];

% 重新计算最终状态
phys = calc_common_physics(p, V_init, sys);
M_spread = calc_M_spread(p, phys, num_electrons);           
M_EB0    = calc_M_EB0(p, phys, sys, err, num_electrons);    
M_Br1    = calc_M_Br1(p, phys, sys, err, num_electrons);    
M_EBrho  = calc_M_EB_angle(p, phys, sys, err, num_electrons); 
M_ee     = calc_M_ee(p, phys, sys, err, num_electrons);     
M_E1     = calc_M_E1(p, phys, sys, err, num_electrons);     
M_B1     = calc_M_B1(p, phys, sys, err, num_electrons);     

M_temp1 = pagemtimes(M_EB0, M_spread);
M_temp2 = pagemtimes(M_Br1, M_temp1);
M_temp3 = pagemtimes(M_EBrho, M_temp2);
M_temp4 = pagemtimes(M_E1, M_temp3);
M_temp5 = pagemtimes(M_B1, M_temp4);
M_total = pagemtimes(M_ee, M_temp5);

V_final = pagemtimes(M_total, V_init);
Pos_final = squeeze(V_final);
X_final = Pos_final(1, :);
Y_final = Pos_final(2, :);

% 提取结果数据用于绘图和报告
cx = mean(X_final);
cy = mean(Y_final);
r_final = sqrt((X_final - cx).^2 + (Y_final - cy).^2);
r_sorted = sort(r_final);
N_total = length(r_sorted);

% 计算各项半径与直径
RMS_radius = sqrt(mean(r_final.^2));
r_50 = r_sorted(round(0.50 * N_total));
r_99 = r_sorted(round(0.99 * N_total));
FW50 = 2 * r_50;
FW99 = 2 * r_99;

%% === 5. 格式化命令行输出报告 ===
total_time_seconds = toc(total_timer);

fprintf('\n==================================================\n');
fprintf('                电子束系统优化结果报告               \n');
fprintf('--------------------------------------------------\n');
fprintf('[ 系统输入参数 ]\n');
fprintf('磁场中心强度 B0 : %.4f T\n', sys.B0);
fprintf('系统极间距 L    : %.3f mm\n', sys.L * 1000);
fprintf('初始光斑半径    : %.2f nm\n', config.spot_radius * 1e9);
fprintf('模拟粒子总数    : %d 个\n', num_electrons);
fprintf('--------------------------------------------------\n');
fprintf('[ 优化输出结果 ]\n');
fprintf('最佳加速电压 U  : %.2f V\n', opt_U);
fprintf('靶面 RMS 半径   : %.2f nm\n', RMS_radius * 1e9);
fprintf('靶面 FW50       : %.2f nm\n', FW50 * 1e9);
fprintf('靶面 FW99       : %.2f nm\n', FW99 * 1e9);
fprintf('--------------------------------------------------\n');
fprintf('总计算耗时      : %.2f 秒\n', total_time_seconds);
fprintf('==================================================\n');

%% === 6. 绘图与结果展示 ===
% --- 靶面光斑分辨率分析图 ---
figure('Name', '优化后靶面电子束斑分辨率分析', 'Position', [100, 100, 1000, 450]);

% 构造包含输入输出参数的通用标题字符串
title_str = sprintf('B_0 = %.4f T, L = %.3f mm  |  最佳 U = %.1f V', sys.B0, sys.L * 1000, opt_U);

% 子图 1: 最终散点图及包络圆圈
subplot(1, 2, 1);
hold on; grid on; axis equal;
scatter(X_final * 1e9, Y_final * 1e9, 2, 'b', '.');
plot(cx * 1e9, cy * 1e9, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
theta = linspace(0, 2*pi, 100);
plot((cx + r_50 * cos(theta)) * 1e9, (cy + r_50 * sin(theta)) * 1e9, 'g-', 'LineWidth', 2);
plot((cx + r_99 * cos(theta)) * 1e9, (cy + r_99 * sin(theta)) * 1e9, 'r--', 'LineWidth', 2);
title({'优化后靶面电子落点分布', title_str}); % 使用 cell 数组实现多行标题
xlabel('X (nm)'); ylabel('Y (nm)');
legend('电子落点', '束流质心', sprintf('FW50 (%.2fnm)', FW50*1e9), sprintf('FW99 (%.2fnm)', FW99*1e9), 'Location', 'best');

% 子图 2: 径向电子累积分布函数 (CDF)
subplot(1, 2, 2);
cumulative_fraction = (1:N_total) / N_total;
plot(r_sorted * 1e9, cumulative_fraction * 100, 'b-', 'LineWidth', 2);
hold on; grid on;
yline(50, 'g--', 'LineWidth', 1.5);
yline(99, 'r--', 'LineWidth', 1.5);
xline(r_50 * 1e9, 'g--', 'LineWidth', 1.5);
xline(r_99 * 1e9, 'r--', 'LineWidth', 1.5);
title({'径向电子流强度累积分布', title_str}); % 同样增加输入参数标识
xlabel('距质心径向距离 r (nm)');
ylabel('包含的电子比例 (%)');
ylim([0 100]);
total_time_seconds = toc(total_timer);
fprintf('程序运行总耗时 : %.4f 秒\n', total_time_seconds);


%% ========================================================================
%  辅助函数：计算特定 U 下的 FW99 
%  (放在脚本末尾或存为单独文件均可，R2016b及以上支持脚本内嵌函数)
%% ========================================================================
function FW99 = eval_FW99(U_test, V_init, p, sys, err, num_electrons)
    % 1. 将当前测试的电压赋值给系统
    sys.U = U_test;
    
    % 2. 动态更新依赖于 U 的电场函数
    sys.func_E_vec = @(X, Y) [
        zeros(size(X));           
        zeros(size(X));                
        (sys.U / sys.L) * ones(size(X))
    ];

    % 3. 预计算物理量
    phys = calc_common_physics(p, V_init, sys);

    % 4. 计算各个转移矩阵
    M_spread = calc_M_spread(p, phys, num_electrons);           
    M_EB0    = calc_M_EB0(p, phys, sys, err, num_electrons);    
    M_Br1    = calc_M_Br1(p, phys, sys, err, num_electrons);    
    M_EBrho  = calc_M_EB_angle(p, phys, sys, err, num_electrons); 
    M_ee     = calc_M_ee(p, phys, sys, err, num_electrons);     
    M_E1     = calc_M_E1(p, phys, sys, err, num_electrons);     
    M_B1     = calc_M_B1(p, phys, sys, err, num_electrons);     

    % 5. 矩阵相乘 (顺序与主脚本保持一致)
    M_temp1 = pagemtimes(M_EB0, M_spread);
    M_temp2 = pagemtimes(M_Br1, M_temp1);
    M_temp3 = pagemtimes(M_EBrho, M_temp2);
    M_temp4 = pagemtimes(M_E1, M_temp3);
    M_temp5 = pagemtimes(M_B1, M_temp4);
    M_total = pagemtimes(M_ee, M_temp5);

    % 6. 计算最终位置
    V_final = pagemtimes(M_total, V_init);
    Pos_final = squeeze(V_final);
    X_final = Pos_final(1, :);
    Y_final = Pos_final(2, :);

    % 7. 计算质心并统计 FW99
    cx = mean(X_final);
    cy = mean(Y_final);
    r_final = sqrt((X_final - cx).^2 + (Y_final - cy).^2);
    
    r_sorted = sort(r_final);
    N_total = length(r_sorted);
    
    % 取 99% 分位数的半径
    r_99 = r_sorted(round(0.99 * N_total));
    
    % 返回直径 FW99 (单位为米，fminbnd直接以此进行优化)
    FW99 = 2 * r_99;
end