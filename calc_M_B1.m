function M_mat = calc_M_B1(p, phys, sys, err, N)
% CALC_M_B1 计算磁场Bz在z方向误差产生的图像缩放转移矩阵
% 输入:
%   p    - 粒子参数
%   phys - 共用物理量 (包含初始位置的 B_vec, X, Y)
%   sys  - 系统参数 (包含极间距 L 以及可能存在的 func_B_vec)
%   err  - 误差参数
%   N    - 电子总数

    % 1. 初始化 3x3xN 单位矩阵
    M_mat = repmat(eye(3), 1, 1, N);
    
    % 3. 获取硅片(靶面)处的磁场 Bz1
    % 如果系统定义了随空间变化的磁场函数，则计算 z=L 处的靶面磁场
    if isfield(sys, 'func_B_vec')
        % 假设 func_B_vec 支持传入 X, Y, Z 三个坐标
        % 如果 func_B_vec 只是 2D 的 func_B_vec(X,Y)，则需要根据你的实际函数修改
        try
            B_target = sys.func_B_vec(phys.X, phys.Y, sys.L * ones(1, N));
            B_start = sys.func_B_vec(phys.X, phys.Y, zeros(1, N));
        catch
            % 兼容旧版的 2D 函数调用
            B_target = sys.func_B_vec(phys.X, phys.Y);
            B_start = sys.func_B_vec(phys.X, phys.Y);
        end
        Bz1 = B_target(3, :);
        Bz0 = B_start(3, :);
    else
        % 如果没有定义空间函数，则认为磁场是绝对均匀的，Bz1 = Bz0
        Bz0 = phys.B_vec(3, :);
        Bz1 = phys.B_vec(3, :);
    end
    
    % 4. 计算缩放比例 M
    % 根据磁动量守恒，r^2 ∝ 1/B，因此 r ∝ 1/sqrt(B)
    % 缩放比 M = sqrt(Bz0 / Bz1)
    M_scale = sqrt(Bz0 ./ Bz1 ); 
    
    % 5. 填入转移矩阵
    % 这是一个标准的缩放矩阵，对角线元素为 M
    % [ M  0  0 ]
    % [ 0  M  0 ]
    % [ 0  0  1 ]
    M_mat(1, 1, :) = M_scale;
    M_mat(2, 2, :) = M_scale;
    
end