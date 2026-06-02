function [V_init, p] = init_particles(num_electrons, config)
% INIT_PARTICLES 初始化电子束的初始状态
% 输入:
%   num_electrons - 需要初始化的电子总数
%   config - 包含各种初始化参数的结构体
% 输出:
%   V_init - 3 x 1 x N 的位置向量矩阵 (齐次坐标 [x; y; 1])
%   p      - 结构体，包含所有电子的运动参数 (U0, alpha, phi)

%% 1. 读取现有文件接口
if isfield(config, 'load_from_file') && config.load_from_file == true
    if exist(config.file_path, 'file')
        fprintf('正在从文件读取粒子数据: %s\n', config.file_path);
        % 假设文件中存有 X, Y, U0, Alpha, Phi 变量
        loaded_data = load(config.file_path);

        % 将加载的数据赋给输出变量
        V_init = loaded_data.V_init;
        p = loaded_data.p;
        return;
    else
        warning('未找到指定的文件，将降级为随机生成模式。');
    end
end

% === 新增：检查并设置默认生成模式 ===
if ~isfield(config, 'pattern_mode')
    config.pattern_mode = "spot";
    fprintf('未指定 pattern_mode，默认使用圆斑模式 (spot)。\n');
end
% ===================================

%% 2. 位置初始化 (5条条纹分布)
fprintf('正在随机生成 %d 个粒子的初始状态...\n', num_electrons);

switch config.pattern_mode
    case "stripe"
        % 解析条纹参数
        n_stripes = 5;
        center_x = config.stripe_center_x; % 条纹整体分布的中心 X
        center_y = config.stripe_center_y; % 条纹整体分布的中心 Y
        width = config.stripe_width;       % 单个条纹的宽度 (X方向)
        length = config.stripe_length;     % 单个条纹的长度 (Y方向)
        spacing = config.stripe_spacing;   % 条纹中心之间的间距

        % 计算每条条纹的中心坐标 (向量化)
        % 例如 5 条条纹，索引为 -2, -1, 0, 1, 2
        stripe_indices = -(n_stripes-1)/2 : (n_stripes-1)/2;
        stripe_centers_x = center_x + stripe_indices * spacing;

        % 步骤 A: 为每个电子随机分配一个条纹 (1 到 5 的随机整数)
        assigned_stripes = randi([1, n_stripes], 1, num_electrons);

        % 步骤 B: 在选定条纹的内部生成局部均匀分布
        local_x = (rand(1, num_electrons) - 0.5) * width;  % [-width/2, width/2]
        local_y = (rand(1, num_electrons) - 0.5) * length; % [-length/2, length/2]

        % 步骤 C: 计算绝对坐标
        x_init = stripe_centers_x(assigned_stripes) + local_x;
        y_init = center_y + local_y;

        % 构建 3x1xN 齐次坐标矩阵
        V_init = zeros(3, 1, num_electrons);
        V_init(1, 1, :) = x_init;
        V_init(2, 1, :) = y_init;
        V_init(3, 1, :) = 1;

    case "spot"
        % 1. 均匀圆斑坐标生成 (核心数学技巧：对 rand 开根号保证面积均匀)
        r_random = config.spot_radius * sqrt(rand(1, num_electrons)); % 径向分布
        theta_random = 2 * pi * rand(1, num_electrons);               % 角向分布

        % 2. 转换为笛卡尔坐标
        X_init = config.spot_center_x + r_random .* cos(theta_random);
        Y_init = config.spot_center_y + r_random .* sin(theta_random);

        % 3. 压入转移矩阵所需的 3x1xN 齐次坐标结构
        V_init = zeros(3, 1, num_electrons);
        V_init(1, 1, :) = X_init;
        V_init(2, 1, :) = Y_init;
        V_init(3, 1, :) = 1;

        otherwise
        % 建议加上 otherwise 以防输入了不支持的模式字符串
        error('不支持的 pattern_mode: %s。请使用 "spot" 或 "stripe"。', config.pattern_mode);
end

%% 3. 能量分布初始化 (0 ~ 0.5V 均匀分布)
p.U0 = config.U0_max * rand(1, num_electrons);

%% 4. 方位角 Phi 初始化 (0 ~ 360° 均匀分布)
p.phi = 2 * pi * rand(1, num_electrons);

%% 5. 发射角 Alpha 初始化 (截断朗伯分布)
% 朗伯分布的 CDF 为 F(alpha) = sin^2(alpha)
% 设 u 为 [0, 1] 的均匀随机数，则 alpha = asin(sqrt(u))
% 为了将 alpha 限制在 alpha_max 内，我们只需让 u 的最大值截断在 sin^2(alpha_max) 即可

alpha_max_rad = deg2rad(config.alpha_max_deg);
u_max = sin(alpha_max_rad)^2;

% 在受限区间内生成均匀随机数
u = u_max * rand(1, num_electrons);

% 反解出 alpha
p.alpha = asin(sqrt(u));

end