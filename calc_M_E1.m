function M = calc_M_E1(p, phys, sys, err, N)
% CALC_M_E1 计算电场Ez的z方向1阶及更高阶误差(离焦)转移矩阵
% 输入:
%   p    - 粒子参数 (包含 phi)
%   phys - 共用物理量 (包含 X, Y, B_vec, tan_gamma)
%   sys  - 系统参数 (包含 N_turn)
%   err  - 误差参数 (包含 z1, z2 以及 空间误差场函数 func_dEz)
%   N    - 电子总数

    % 1. 初始化 3x3xN 转移矩阵
    M = repmat(eye(3), 1, 1, N);
    
    e = 1.602e-19;
    m = 9.109e-31;
    Bz = phys.B_vec(3, :); % 取出主磁场 (1xN 向量)
    
    % 2. 核心：向量化处理 z 方向的积分求平均 (E_average)
    if isfield(err, 'func_dEz') && isfield(err, 'z1') && isfield(err, 'z2')
        % 设定积分步数 (根据场变化的剧烈程度调整，一般 20~50 足够了)
        n_steps = 100; 
        z_nodes = linspace(err.z1, err.z2, n_steps); 
        
        % 预分配空间: n_steps 行, N 列
        dEz_matrix = zeros(n_steps, N);
        
        % 批量计算每层 Z 切片上的电场误差
        % phys.X 和 phys.Y 是一维向量 (1xN)
        for i = 1:n_steps
            % 调用外部定义的电场误差函数，传入 X, Y 和 当前层的高度 Z
            dEz_matrix(i, :) = err.func_dEz(phys.X, phys.Y, z_nodes(i));
        end
        
        % 数值积分取平均值: \int dEz dz / \delta z  ≈ mean(dEz_matrix)
        % 结果 E_average 是一个 1xN 向量
        E_average = mean(dEz_matrix, 1); 
    else
        % 如果未定义复杂空间场，则假设误差为 0
        E_average = zeros(1, N);
    end
    
    % 3. 计算焦面偏移距离 δL_E1 
    % δLE1 = 2 * (π^2 * N^2 * m) / (e * B^2) * Eaverage
    delta_LE1 = 2 * (pi^2 * sys.N_turn^2 * m) ./ (e * Bz.^2 ) .* E_average;
    
    % 4. 计算对线宽的影响 δE1 (按电子束收束角计算) 
    % δE1 = 2 * δLE1 * tan(γ)
    delta_E1 = 2 .* delta_LE1 .* phys.tan_gamma;
    
    % 5. 投影到 x 和 y 方向 
    dx = delta_E1 .* cos(p.phi);
    dy = delta_E1 .* sin(p.phi);
    
    % 6. 填入转移矩阵平移列 
    M(1, 3, :) = dx;
    M(2, 3, :) = dy;
end