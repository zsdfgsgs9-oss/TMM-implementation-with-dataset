function M = calc_M_ee(p, phys, sys, err, N)
% CALC_M_EE 计算电子间相互作用(图案扩大)的转移矩阵
% 支持 'sum' (10000步离散求和) 和 'gauss' (50阶高斯求积) 两种方法对比
% 输入:
%   p    - 粒子参数
%   phys - 共用物理量 (包含 t, v1, U_local, L_local, omega, R)
%   sys  - 系统参数
%   err  - 误差参数 (需包含 I_beam, R0, 以及方法开关 ee_method)
%   N    - 电子总数

    % 1. 初始化 3x3xN 单位矩阵
    M = repmat(eye(3), 1, 1, N);
    
    % 物理常量
    eps0 = 8.854e-12;
    e = 1.602e-19;
    m = 9.109e-31;
    
    % 从 err 提取系统参数
    I = err.I_beam; 
    R0 = err.R0; 
    
    % 提取并计算基础运动学向量 (1 x N 向量)
    t_flight = phys.t; 
    a_z = (e .* phys.U_local) ./ (m .* phys.L_local );
    
    % 默认使用 gauss 方法，如果指定了则读取指定方法
    if isfield(err, 'ee_method')
        method = err.ee_method;
    else
        method = 'gauss';
    end
    
    % 2. 核心算法分支
    tic; % 开始计时，方便你对比速度
    switch lower(method)
        case 'sum'
            % ==== 方法一：10000 步矩阵求和法 ====
            n_steps = 10000;
            tau = linspace(0, 1, n_steps)'; % 10000 x 1 列向量
            dt_vec = t_flight ./ n_steps;   % 1 x N 行向量
            
            % 构建时间网格 (10000 x N 矩阵，利用隐式展开)
            t_matrix = tau * t_flight; 
            
            % 运动学演化
            v_z_matrix = phys.v1 + a_z .* t_matrix;
            R_matrix = R0 + 2 .* phys.R .* sin(phys.omega .* t_matrix ./ 2);
            E_perp_matrix = I ./ (2 * pi * eps0 .* R_matrix .* v_z_matrix );
            
            % 离散求和：sum( (T - t_current) * E_perp * dt )
            T_minus_t = t_flight - t_matrix;
            delta_ee = (e / m) .* sum(E_perp_matrix .* T_minus_t, 1) .* dt_vec;
            
            fprintf('空间电荷 (求和法 10000步) 计算完毕，耗时: %.4f 秒\n', toc);
            
        case 'gauss'
            % ==== 方法二：50 阶高斯-勒让德求积法 ====
            n_pts = 50;
            
            % 获取 [0, 1] 区间的高斯节点 (tau) 和权重 (w)
            [tau, w] = get_gauss_nodes(n_pts); 
            % tau 是 n_pts x 1 列向量, w 是 1 x n_pts 行向量
            
            % 构建时间网格 (50 x N 矩阵)
            t_matrix = tau * t_flight; 
            
            % 运动学演化 (与求和法公式完全一致，但维度极小)
            v_z_matrix = phys.v1 + a_z .* t_matrix;
            R_matrix = R0 + 2 .* phys.R .* sin(phys.omega .* t_matrix ./ 2);
            E_perp_matrix = I ./ (2 * pi * eps0 .* R_matrix .* v_z_matrix );
            
            % 计算被积函数: (1 - tau) * E_perp
            % 隐式展开: tau(50x1) 扩展后与 E_perp_matrix(50xN) 点乘
            integrand_matrix = (1 - tau) .* E_perp_matrix;
            
            % 矩阵相乘实现高斯积分求和: w(1x50) * integrand(50xN) = (1xN)
            integral_result = w * integrand_matrix;
            
            % 乘以外部系数: (e/m) * T^2 * 积分结果
            delta_ee = (e / m) .* (t_flight.^2) .* integral_result;
            
            fprintf('空间电荷 (高斯法 50阶) 计算完毕，耗时: %.4f 秒\n', toc);
            
        otherwise
            error('未知的电子间相互作用计算方法，请设定 err.ee_method 为 ''sum'' 或 ''gauss''');
    end
    
    % 3. 赋予随机的扩散方向并填入转移矩阵
    epsilon = 2 * pi * rand(1, N);
    
    M(1, 3, :) = delta_ee .* cos(epsilon);
    M(2, 3, :) = delta_ee .* sin(epsilon);
end

% ==========================================================
% 内部子函数：生成 [0, 1] 区间的高斯-勒让德节点和权重
% ==========================================================
function [x_mapped, w_mapped] = get_gauss_nodes(N_pts)
    % 使用 Golub-Welsch 算法计算标准 [-1, 1] 区间的节点和权重
    i = 1:N_pts-1;
    a = i ./ sqrt(4*i.^2 - 1);
    CM = diag(a, 1) + diag(a, -1); % 构建伴随矩阵
    
    [V, L] = eig(CM); % 求特征值和特征向量
    [x, ind] = sort(diag(L));
    V = V(:, ind)';
    w = 2 * V(:, 1).^2;
    
    % 映射到 [0, 1] 区间
    % x_new = (b-a)/2 * x + (b+a)/2
    % w_new = (b-a)/2 * w
    x_mapped = (x + 1) / 2;
    w_mapped = (w / 2)'; % 转置为 1 x n_pts 的行向量
end