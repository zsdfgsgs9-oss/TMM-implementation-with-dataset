function M = calc_M_spread(p, phys, N)
% CALC_M_SPREAD 计算发散角和能散的转移矩阵
% 输入:
%   p    - 粒子参数结构体，包含 p.phi (1xN 向量)
%   phys - 共用物理量结构体，包含 phys.R 和 phys.beta (1xN 向量)
%   N    - 电子总数
% 输出:
%   M    - 3x3xN 的转移矩阵数组

    % 1. 初始化 3x3xN 单位矩阵
    M = repmat(eye(3), 1, 1, N);
    
    % 2. 提取预计算好的物理量 (均为 1xN 向量)
    R = phys.R;       % 局部回旋半径
    beta = phys.beta; % 局部旋转相位
    phi = p.phi;      % 发射方位角
    
    % 3. 向量化计算位移项 δx 和 δy
    % 注意这里直接使用向量进行加减和三角函数运算
    dx = R .* (sin(phi + beta) - sin(phi));
    dy = R .* (cos(phi) - cos(phi + beta));
    
    % 4. 批量填入转移矩阵的平移列 (第1行第3列，第2行第3列)
    M(1, 3, :) = dx; 
    M(2, 3, :) = dy;
    
end