function M = calc_M_Br1(p, phys, sys, err, N)
% CALC_M_BR1 计算磁场Br的r方向1阶误差(旋转畸变)转移矩阵
% 输入:
%   p    - 粒子参数 (包含 U0, alpha 1xN 向量)
%   phys - 共用物理量 (包含 B_vec, t, X, Y 1xN 向量)
%   sys  - 系统参数
%   err  - 误差参数 (包含 G 或者 func_G)
%   N    - 电子总数

    % 1. 初始化 3x3xN 单位矩阵
    M = repmat(eye(3), 1, 1, N);
    
    % 2. 提取局部主磁场 Bz (1xN 向量)
    % 这里的 Bz 是电子发射位置的局部磁场
    Bz = phys.B_vec(3, :); 
    
    % 3. 提取磁场 z 方向的梯度 G
    if isfield(err, 'func_G')
        % 兼容高端用法：如果梯度本身也是空间分布的函数
        G_local = err.func_G(phys.X, phys.Y);
    elseif isfield(err, 'G')
        % 基础用法：全局统一的梯度常量
        G_local = err.G * ones(1, N);
    else
        % 如果未定义，认为无畸变
        G_local = zeros(1, N);
    end
    
    % 4. 计算旋转畸变系数 A_vec
    % 通过引入 phys.t，避免了文档源公式 [176] 中极其冗长的时间展开项
    % A = (U0 * sin^2(alpha) / Bz^3) * G^2 * t
    term1 = (p.U0 .* sin(p.alpha).^2) ./ Bz.^3 * 0.75;
    A_vec = term1 .* (G_local.^2) .* phys.t;
    
    % 5. 填入转移矩阵 [cite: 174, 175]
    % 矩阵形式为：
    % [  1  -A   0 ]
    % [  A   1   0 ]
    % [  0   0   1 ]
    M(1, 2, :) = -A_vec;
    M(2, 1, :) = A_vec;
    
end