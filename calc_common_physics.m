function phys = calc_common_physics(p, V_init, sys)
    % 提取电子当前的绝对坐标 (1xN 向量)
    X = squeeze(V_init(1, 1, :)).';
    Y = squeeze(V_init(2, 1, :)).';
    N = length(p.U0);
    
    % 1. 获取 3xN 的电磁场向量
    if isfield(sys, 'func_B_vec')
        Z_target = 0.5 * sys.L * ones(1, N);
        phys.B_vec = sys.func_B_vec(X, Y, Z_target);
        phys.E_vec = sys.func_E_vec(X, Y);
    else
        % 默认匀强场
        phys.B_vec = repmat(sys.B_vec_0, 1, N);
        phys.E_vec = repmat(sys.E_vec_0, 1, N);
    end
    
    % 假定电压也可能具有空间分布，若没有则用常数
    if isfield(sys, 'func_U')
        phys.U_local = sys.func_U(X, Y);
    else
        phys.U_local = sys.U * ones(1, N);
    end

    if isfield(sys, 'func_L')
        phys.L_local = sys.func_L(X, Y);
    else
        phys.L_local = sys.L * ones(1, N);
    end
    
    % 2. 提取主轴分量用于 0 阶轨迹计算
    Bz = phys.B_vec(3, :); % 1xN，提取 z 方向磁场
    
    e = 1.602e-19;
    m = 9.109e-31;
    
    % 3. 核心 0 阶物理量计算 (使用 Bz 代替原来的标量 B)
    phys.v1 = sqrt(2 * e * p.U0 / m) .* cos(p.alpha);
    
    % 飞行时间 t (电场 Ez 的影响已经体现在 U_local 中)
    term_sqrt = sqrt(1 + (2 * e * phys.U_local) ./ (m * phys.v1.^2 ));
    phys.t = (m .* phys.v1 .* sys.L ./ (e .* phys.U_local )) .* (term_sqrt - 1);
    
    % 回旋角速度与相位 (仅由 Bz 决定 0 阶螺旋) [cite: 10]
    phys.omega = e .* Bz ./ m;
    phys.beta = phys.omega .* phys.t;
    
    % 回旋半径 R [cite: 19]
    phys.R = (m ./ (e .* Bz )) .* sqrt(2 * e * p.U0 / m) .* sin(p.alpha);
    
    % 计算 B 的模长平方 (用于后续 E x B / B^2 计算) 
    phys.B_norm_sq = sum(phys.B_vec.^2, 1); % 1xN 向量

    % 收束角正切 tan_gamma (用于 0 阶误差模糊度计算)
    phys.tan_gamma = (sqrt(p.U0) .* sin(p.alpha)) ./ (sqrt(p.U0 .* cos(p.alpha).^2 + phys.U_local) );
    
    % 保存坐标供后续调用
    phys.X = X;
    phys.Y = Y;
end