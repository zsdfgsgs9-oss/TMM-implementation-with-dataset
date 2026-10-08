function phys = calc_common_physics(p, V_init, sys)
    % 提取电子当前的绝对坐标 (1xN 向量)
    X = squeeze(V_init(1, 1, :)).';
    Y = squeeze(V_init(2, 1, :)).';
    N = length(p.U0);

    if isfield(sys, 'physics_model')
        model = lower(char(string(sys.physics_model)));
    elseif isfield(sys, 'use_relativistic') && sys.use_relativistic
        model = 'relativistic';
    else
        model = 'nonrelativistic';
    end
    is_relativistic = any(strcmp(model, {'relativistic', 'rel', 'sr'}));
    
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
    c = 299792458;
    E_nominal = phys.U_local ./ phys.L_local;
    
    % 3. 非相对论 0 阶物理量
    nonrel.v1 = sqrt(2 * e * p.U0 / m) .* cos(p.alpha);
    nonrel.v_perp0 = sqrt(2 * e * p.U0 / m) .* sin(p.alpha);
    nonrel.a_z = e .* phys.U_local ./ (m .* phys.L_local);
    
    % 飞行时间 t。该写法与原公式等价，但在 v1 很小时更稳定。
    nonrel.t = (-nonrel.v1 + sqrt(nonrel.v1.^2 + 2 .* nonrel.a_z .* phys.L_local)) ./ nonrel.a_z;
    
    % 回旋角速度与相位。非均匀 Bz 可选择解析或 Gauss 路径积分。
    nonrel.omega = e .* Bz ./ m;
    [nonrel.beta, nonrel.phase_diagnostics] = ...
        calc_nonrel_magnetic_phase(nonrel.v1, nonrel.a_z, nonrel.t, ...
        phys.L_local, Bz, X, Y, sys);
    
    % 回旋半径 R [cite: 19]
    nonrel.R = (m ./ (e .* Bz )) .* nonrel.v_perp0;
    nonrel.tan_gamma = (sqrt(p.U0) .* sin(p.alpha)) ./ ...
        (sqrt(p.U0 .* cos(p.alpha).^2 + phys.U_local));

    % 4. 相对论 0 阶物理量
    rel.gamma0 = 1 + (e .* p.U0) ./ (m * c^2);
    rel.p0 = m * c .* sqrt(max(rel.gamma0.^2 - 1, 0));
    rel.p_parallel0 = rel.p0 .* cos(p.alpha);
    rel.p_perp0 = rel.p0 .* sin(p.alpha);
    rel.a = sqrt(1 + (rel.p_perp0 ./ (m * c)).^2);
    rel.gamma_L = rel.gamma0 + (e .* phys.U_local) ./ (m * c^2);
    rel.p_zL = m * c .* sqrt(max(rel.gamma_L.^2 - rel.a.^2, 0));
    rel.t = (rel.p_zL - rel.p_parallel0) ./ (e .* E_nominal);
    rel.s0 = asinh(rel.p_parallel0 ./ (m * c .* rel.a));
    rel.sL = asinh(rel.p_zL ./ (m * c .* rel.a));
    rel.beta_uniform = (Bz .* c ./ E_nominal) .* (rel.sL - rel.s0);
    rel.beta = rel.beta_uniform;

    if isfield(sys, 'func_B_vec')
        n_phase_steps = 80;
        if isfield(sys, 'rel_phase_steps')
            n_phase_steps = sys.rel_phase_steps;
        end

        delta_beta_integral = zeros(1, N);
        dz = phys.L_local ./ n_phase_steps;
        for i = 1:n_phase_steps
            tau = (i - 0.5) / n_phase_steps;
            Z_current = tau .* phys.L_local;
            B_vec_local = sys.func_B_vec(X, Y, Z_current);
            Bz_local = B_vec_local(3, :);
            gamma_z = rel.gamma0 + (e .* E_nominal .* Z_current) ./ (m * c^2);
            p_z = m * c .* sqrt(max(gamma_z.^2 - rel.a.^2, 0));
            phase_density = e .* (Bz_local - Bz) ./ p_z;
            phase_density(~isfinite(phase_density)) = 0;
            delta_beta_integral = delta_beta_integral + phase_density .* dz;
        end
        rel.beta = rel.beta_uniform + delta_beta_integral;
    end

    rel.R = rel.p_perp0 ./ (e .* Bz);
    rel.tan_gamma = rel.p_perp0 ./ rel.p_zL;
    rel.tan_gamma(~isfinite(rel.tan_gamma)) = 0;
    rel.v_z0 = rel.p_parallel0 ./ (rel.gamma0 .* m);
    rel.v_zL = rel.p_zL ./ (rel.gamma_L .* m);
    rel.omega0 = e .* Bz ./ (rel.gamma0 .* m);
    rel.E_nominal = E_nominal;
    rel.c = c;
    
    % 计算 B 的模长平方 (用于后续 E x B / B^2 计算) 
    phys.B_norm_sq = sum(phys.B_vec.^2, 1); % 1xN 向量

    phys.nonrel = nonrel;
    phys.rel = rel;
    phys.model = model;
    phys.is_relativistic = is_relativistic;

    if is_relativistic
        phys.v1 = rel.v_z0;
        phys.t = rel.t;
        phys.omega = rel.omega0;
        phys.beta = rel.beta;
        phys.R = rel.R;
        phys.tan_gamma = rel.tan_gamma;
    else
        phys.v1 = nonrel.v1;
        phys.t = nonrel.t;
        phys.omega = nonrel.omega;
        phys.beta = nonrel.beta;
        phys.R = nonrel.R;
        phys.tan_gamma = nonrel.tan_gamma;
    end
    
    % 保存坐标供后续调用
    phys.X = X;
    phys.Y = Y;
end
