function M = calc_M_EB_angle(p, phys, sys, err, N)
% CALC_M_EB_ANGLE 计算电场与磁场夹角造成的误差(拉伸畸变)转移矩阵
% 输入:
%   p    - 粒子参数
%   phys - 共用物理量 (包含 X, Y, E_vec, t)
%   sys  - 系统参数 (包含极间距 L, func_B_vec)
%   err  - 误差参数
%   N    - 电子总数

    % 1. 初始化 3x3xN 单位矩阵
    M = repmat(eye(3), 1, 1, N);
    
    % 2. 提取局部电场 E_vec (假设 E 只与 X,Y 有关，所以在整个 Z 路径上不变)
    % E_vec 是一个 3xN 的矩阵
    E_vec = phys.E_vec; 
    
    % 3. 沿 Z 轴进行离散切片计算平均漂移速度
    n_steps = 50; % 积分步数 (可根据精度要求调整)
    z_nodes = linspace(0, sys.L, n_steps);
    
    % 预分配空间用于存储每一层的漂移速度 (3xN)
    v_drift_avg = zeros(3, N);
    
    for i = 1:n_steps
        % 获取当前 Z 切片的绝对高度向量 (1xN)
        Z_current = z_nodes(i) * ones(1, N);
        
        % 调用三维磁场函数获取当前层的 B_vec_local (3xN 矩阵)
        if isfield(sys, 'func_B_vec')
            B_vec_local = sys.func_B_vec(phys.X, phys.Y, Z_current);
        else
            B_vec_local = repmat([0; 0; sys.B], 1, N); % 降级为匀强磁场
        end
        
        % 计算当前层的 B^2 (1xN 向量)
        B_sq_local = sum(B_vec_local.^2, 1);
        
        % 核心物理：计算 E x B 漂移速度 (3xN 矩阵)
        % cross 函数会按列自动计算 10000 个电子的叉乘，极其高效！
        v_drift_local = cross(E_vec, B_vec_local) ./ B_sq_local ;
        
        % 累加求平均
        v_drift_avg = v_drift_avg + v_drift_local / n_steps;
    end
    
    % 4. 计算总漂移位移 (速度乘以时间)
    % phys.t 是之前预计算好的精确飞行时间 (1xN)
    delta_r = v_drift_avg .* phys.t;
    
    % 5. 提取 x 和 y 方向的位移量并填入转移矩阵
    % delta_r(1,:) 为 x 方向漂移，delta_r(2,:) 为 y 方向漂移
    M(1, 3, :) = delta_r(1, :);
    M(2, 3, :) = delta_r(2, :);
    
end