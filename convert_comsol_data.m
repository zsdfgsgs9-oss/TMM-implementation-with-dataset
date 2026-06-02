% convert_comsol_data.m
% 用于读取 COMSOL 导出的电子数据，并转化为转移矩阵程序所需的 V_init 和 p 结构

clc; clear;

%% 1. 读取 COMSOL 导出的 CSV 文件
fprintf('正在读取 COMSOL 数据...\n');
filename = '28_second_order_2H_initial.csv';
% readmatrix 能够自动忽略以 '%' 开头的 COMSOL 注释行
data = readmatrix(filename);

% data 的列对应关系：
% 第 1 列: 索引
% 第 2 列: qx (um)
% 第 3 列: qy (um)
% 第 4 列: 动能 Ep (J)
% 第 5 列: vx (m/s)
% 第 6 列: vy (m/s)
% 第 7 列: vz (m/s)

N = size(data, 1); % 电子总数

%% 2. 提取并转换空间坐标 (转换为米 m)
% COMSOL 导出的单位是微米 (um)，需乘以 1e-6 统一到国际标准单位
x = data(:, 2)' * 1e-6; % 转为 1xN 行向量
y = data(:, 3)' * 1e-6; 

dx = 0;
dy = 0;

% 构造 3x1xN 的齐次坐标矩阵 V_init
V_init = zeros(3, 1, N);
V_init(1, 1, :) = x + dx;
V_init(2, 1, :) = y + dy;
V_init(3, 1, :) = 1; % 第 3 维度设为 1（转移矩阵齐次坐标要求）

%% 3. 提取速度分量与能量
e = 1.602e-19;
Ep_Joule = data(:, 4)'; % 1xN 行向量 (焦耳)

% 能量 U0 转为电子伏特 (eV)
p.U0 = Ep_Joule / e;

% 提取速度分量 (m/s)
vx = data(:, 5)';
vy = data(:, 6)';
vz_comsol = data(:, 7)';

%% 4. 坐标系匹配与角度计算
% 【核心操作】: 因为 COMSOL 中电子沿 -z 方向发射，而 MATLAB 模型主轴是 +z
% 必须将 z 方向的速度反转
vz_new = vz_comsol;

% 横向速度大小
v_xy = sqrt(vx.^2 + vy.^2);

% 计算发射发散角 alpha (与新 z 轴的夹角，范围 0 到 pi/2)
% atan2(Y, X) 可确保不会出现相限错误
p.alpha = atan2(v_xy, vz_new);

% 计算方位角 phi (在 x-y 平面内的投影角，范围 -pi 到 pi)
p.phi = atan2(vy, vx);

%% 5. 验证与保存
% 简单剔除由于计算精度可能产生的 NaN（安全保障）
invalid_idx = isnan(p.U0) | isnan(p.alpha) | isnan(p.phi);
if any(invalid_idx)
    fprintf('警告：发现 %d 个无效数据，正在清理...\n', sum(invalid_idx));
    V_init(:,:,invalid_idx) = [];
    p.U0(invalid_idx) = [];
    p.alpha(invalid_idx) = [];
    p.phi(invalid_idx) = [];
    N = length(p.U0);
end

% 保存为 MATLAB 数据文件
save('existing_particles.mat', 'V_init', 'p');
fprintf('成功将 %d 个 COMSOL 电子转化为矩阵格式！\n', N);
fprintf('文件已保存为 existing_particles.mat\n');

%% (可选) 画一张简单的数据验证直方图
figure('Name', 'COMSOL 电子初始参数分布');
subplot(1,3,1); histogram(x * 1e9, 50); title('x 坐标分布 (nm)');
subplot(1,3,2); histogram(p.U0, 50); title('初始动能 U_0 (eV)');
subplot(1,3,3); histogram(rad2deg(p.alpha), 50); title('发射发散角 \alpha (deg)');