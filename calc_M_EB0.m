function M = calc_M_EB0(p, phys, sys, err, N)
% CALC_M_EB0 计算电场Ez和磁场Bz的0阶误差(离焦)转移矩阵
% 输入:
%   p    - 粒子参数结构体 (包含 phi)
%   phys - 共用物理量结构体 (包含 B_vec, U_local, v1, tan_gamma)
%   sys  - 系统参数 (包含 L, N_turn)
%   err  - 误差参数 (包含 delta_E, delta_B)
%   N    - 电子总数

    % 1. 初始化 3x3xN 单位矩阵
    M = repmat(eye(3), 1, 1, N);
    
    % 物理常量
    e = 1.602e-19;
    m = 9.109e-31;
    
    % 2. 提取局部变量 (均为 1xN 向量)
    Bz = phys.B_vec(3, :); % 局部 z 方向主磁场
    E_nominal = phys.U_local ./ sys.L; % 标称匀强电场 E = U/L
    
    % 3. 计算电场 0 阶误差引起的焦面偏移 delta_L_E
    % 公式: δL_E = (2 * π^2 * N^2 * m) / (e * B^2) * δE
    delta_L_E = (2 * pi^2 * sys.N_turn^2 * m) ./ (e * Bz.^2 ) .* err.delta_E;
    
    % 4. 计算磁场 0 阶误差引起的焦面偏移 delta_L_B
    % 公式: δL_B = -(δB / B) * [ (4*π^2*N^2*m*E)/(e*B^2) + (2*π*N*m*v1*cosα)/(e*B^2) ]
    % 注意文档中的 sqrt(2eU0/m)*cosα 其实就是垂直初速度 v1
    term1_B = (4 * pi^2 * sys.N_turn^2 * m .* E_nominal) ./ (e * Bz.^2 );
    term2_B = (2 * pi * sys.N_turn * m .* phys.v1) ./ (e * Bz.^2 );
    delta_L_B = -(err.delta_B ./ (Bz )) .* (term1_B + term2_B);
    
    % 5. 计算靶面上的总模糊半径 delta_EB0
    % 公式: δEB0 = 2 * (δL_E + δL_B) * tan(γ)
    delta_EB0 = (delta_L_E + delta_L_B) .* phys.tan_gamma;
    
    % 6. 分解到 x 和 y 方向的位移
    % 离焦导致的扩散方向即为电子初始的方位角 phi
    dx = delta_EB0 .* cos(p.phi);
    dy = delta_EB0 .* sin(p.phi);
    
    % 7. 填入转移矩阵
    M(1, 3, :) = dx;
    M(2, 3, :) = dy;
end