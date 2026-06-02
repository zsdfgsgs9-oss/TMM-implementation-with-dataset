% analyze_stripes_11_16.m
% 对文件11-16的条纹图案进行后处理
% 计算指标：质心偏移量、FWHM线宽、旋转角
% 分别对 TMM (全误差开启) 和 FEM (COMSOL) 结果进行计算并对比
%
% 关键设定:
%   - 所有误差项全部开启 (去掉 *0)
%   - 文件名中的 1mm/1μm 等为初始条纹中心 x 坐标 (标注 y 的为 y 坐标)
%   - 初始状态 CSV 需叠加对应偏移量后再送入 TMM

clear; clc; close all;

%% === 系统参数 (全误差开启) ===
sys.B0 = 1.2;
sys.L = 1e-3;
sys.U = 12.67e3;
sys.N_turn = 1;

% 磁场梯度 (原 main_process.m: err.G = 6*0)
err.G = 6;

sys.func_B_vec = @(X, Y, Z) [
    -0.5 * err.G .* X;                                        % Bx
    -0.5 * err.G .* Y;                                        % By
    sys.B0 + err.G .* (Z - 0.5 * sys.L) .* ones(size(X))     % Bz
    ];

% 电场 (原 main_process.m: Ex, Ey 后面有 *0, 现全部打开)
sys.func_E_vec = @(X, Y) [
    (sys.U / sys.L) * ones(size(X)) * 1e-3 * 0.8;            % Ex
    (sys.U / sys.L) * ones(size(X)) * 1e-3 * 0.6;            % Ey
    (sys.U / sys.L) * ones(size(X))                           % Ez
    ];

% 电场0阶误差 (原: *0)
err.delta_E = 1.267e4;
% 磁场0阶误差 (原: *0)
err.delta_B = 1.2e-3;

err.z1 = 0;
err.z2 = 1e-3;

% 电场z向高阶误差 (原: *0)
err.func_dEz = @(x, y, z) -6.335e7 * z;

err.I_beam = 3e-15;
err.R0 = 18e-9;
err.ee_method = 'gauss';

%% === 文件列表: {文件名前缀, 初始条纹中心x偏移(m), 初始条纹中心y偏移(m)} ===
file_entries = {
    '11_system_stripe_center',     0,       0;
    '12_system_stripe_1mm',        1e-3,    0;
    '13_system_stripe_1um',        1e-6,    0;
    '14_system_stripe_0p1mm',      0.1e-3,  0;
    '15_system_stripe_y1mm',       0,       1e-3;
    '16_system_stripe_0p5mm',      0.5e-3,  0;
    };

n_files = size(file_entries, 1);
results = struct();

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════════════════╗\n');
fprintf('║       条纹后处理分析 (Files 11-16)  TMM (全误差) vs FEM (COMSOL)       ║\n');
fprintf('╚══════════════════════════════════════════════════════════════════════════╝\n');

%% === 主循环 ===
for f_idx = 1:n_files
    base_name   = file_entries{f_idx, 1};
    offset_cx   = file_entries{f_idx, 2};  % 初始条纹中心 x 偏移
    offset_cy   = file_entries{f_idx, 3};  % 初始条纹中心 y 偏移

    fprintf('\n──────────────────────────────────────────────────────────────────────────\n');
    fprintf('  [%d/%d] %s\n', f_idx, n_files, base_name);
    fprintf('  初始条纹中心偏移: x = %.4f mm,  y = %.4f mm\n', offset_cx * 1e3, offset_cy * 1e3);
    fprintf('──────────────────────────────────────────────────────────────────────────\n');

    init_file  = [base_name '_initial.csv'];
    image_file = [base_name '_image.csv'];

    % =======================================================================
    % PART A: 读取 COMSOL 初始状态，叠加偏移量
    % =======================================================================
    if ~exist(init_file, 'file')
        fprintf('  [警告] 找不到初始状态文件: %s，跳过。\n', init_file);
        continue;
    end

    data_init = readmatrix(init_file);
    N_particles = size(data_init, 1);

    % 坐标转换 μm → m，并叠加偏移
    x_init = data_init(:, 2)' * 1e-6 + offset_cx;
    y_init = data_init(:, 3)' * 1e-6 + offset_cy;

    % 构造齐次坐标
    V_init_tmm = zeros(3, 1, N_particles);
    V_init_tmm(1, 1, :) = x_init;
    V_init_tmm(2, 1, :) = y_init;
    V_init_tmm(3, 1, :) = 1;

    % 物理量提取
    e_charge = 1.602e-19;
    Ep_Joule = data_init(:, 4)';
    p_tmm.U0 = Ep_Joule / e_charge;
    vx = data_init(:, 5)';
    vy = data_init(:, 6)';
    vz = data_init(:, 7)';
    v_xy = sqrt(vx.^2 + vy.^2);
    p_tmm.alpha = atan2(v_xy, vz);
    p_tmm.phi   = atan2(vy, vx);

    % 清理无效数据
    invalid_idx = isnan(p_tmm.U0) | isnan(p_tmm.alpha) | isnan(p_tmm.phi);
    if any(invalid_idx)
        V_init_tmm(:, :, invalid_idx) = [];
        p_tmm.U0(invalid_idx) = [];
        p_tmm.alpha(invalid_idx) = [];
        p_tmm.phi(invalid_idx) = [];
        N_particles = length(p_tmm.U0);
    end

    fprintf('  粒子数: %d\n', N_particles);

    % =======================================================================
    % PART B: 运行 TMM 转移矩阵流水线 (全误差开启)
    % =======================================================================
    phys = calc_common_physics(p_tmm, V_init_tmm, sys);

    M_spread = calc_M_spread(p_tmm, phys, N_particles);
    M_EB0    = calc_M_EB0(p_tmm, phys, sys, err, N_particles);
    M_Br1    = calc_M_Br1(p_tmm, phys, sys, err, N_particles);
    M_EBrho  = calc_M_EB_angle(p_tmm, phys, sys, err, N_particles);
    M_ee     = calc_M_ee(p_tmm, phys, sys, err, N_particles);
    M_E1     = calc_M_E1(p_tmm, phys, sys, err, N_particles);
    M_B1     = calc_M_B1(p_tmm, phys, sys, err, N_particles);

    M_temp1 = pagemtimes(M_EB0, M_spread);
    M_temp2 = pagemtimes(M_Br1, M_temp1);
    M_temp3 = pagemtimes(M_E1, M_temp2);
    M_temp4 = pagemtimes(M_B1, M_temp3);
    M_temp5 = pagemtimes(M_EBrho, M_temp4);
    M_total = pagemtimes(M_ee, M_temp5);

    V_final_tmm = pagemtimes(M_total, V_init_tmm);
    Pos_tmm = squeeze(V_final_tmm);
    % 减去输入位置偏移，使 TMM 结果仅包含像差/漂移效应，与 FEM 直接可比
    X_tmm = Pos_tmm(1, :) - offset_cx;
    Y_tmm = Pos_tmm(2, :) - offset_cy;

    % =======================================================================
    % PART C: 读取 COMSOL FEM 图像结果
    % =======================================================================
    if ~exist(image_file, 'file')
        fprintf('  [警告] 找不到图像文件: %s，跳过FEM分析。\n', image_file);
        X_fem = [];
        Y_fem = [];
    else
        data_img = readmatrix(image_file);
        % COMSOL 导出单位 μm，坐标翻转与 main_process.m 一致
        X_fem =  data_img(:, 1)' * 1e-6;
        Y_fem =  data_img(:, 2)' * 1e-6;
    end

    % =======================================================================
    % PART D: 条纹指标计算
    % =======================================================================
    [cx_tmm, cy_tmm, fwhm_tmm, rot_tmm, fwhm_std_tmm, fwhm_ind_tmm, ...
     fw90_tmm, fw90_std_tmm, fw90_ind_tmm, n_str_tmm] = ...
        analyze_stripe_pattern(X_tmm, Y_tmm);

    if ~isempty(X_fem)
        [cx_fem, cy_fem, fwhm_fem, rot_fem, fwhm_std_fem, fwhm_ind_fem, ...
         fw90_fem, fw90_std_fem, fw90_ind_fem, n_str_fem] = ...
            analyze_stripe_pattern(X_fem, Y_fem);
    end

    fprintf('\n  ┌─────────────────────────────────────────────────────────────────┐\n');
    fprintf('  │                      条纹量化指标结果                            │\n');
    fprintf('  ├──────────────┬──────────────────────┬──────────────────────────┤\n');
    fprintf('  │     指标      │     TMM (全误差)     │      FEM (COMSOL)        │\n');
    fprintf('  ├──────────────┼──────────────────────┼──────────────────────────┤\n');

    fprintf('  │ 识别条纹数   │  %10d 条         │', n_str_tmm);
    if ~isempty(X_fem)
        fprintf('  %10d 条         │\n', n_str_fem);
    else
        fprintf('         N/A            │\n');
    end

    fprintf('  │ FWHM 均值    │  %10.2f nm        │', fwhm_tmm * 1e9);
    if ~isempty(X_fem)
        fprintf('  %10.2f nm        │\n', fwhm_fem * 1e9);
    else
        fprintf('         N/A            │\n');
    end

    fprintf('  │ FWHM 标准差  │  %10.2f nm        │', fwhm_std_tmm * 1e9);
    if ~isempty(X_fem)
        fprintf('  %10.2f nm        │\n', fwhm_std_fem * 1e9);
    else
        fprintf('         N/A            │\n');
    end

    fprintf('  │ FW90 均值    │  %10.2f nm        │', fw90_tmm * 1e9);
    if ~isempty(X_fem)
        fprintf('  %10.2f nm        │\n', fw90_fem * 1e9);
    else
        fprintf('         N/A            │\n');
    end

    fprintf('  │ FW90 标准差  │  %10.2f nm        │', fw90_std_tmm * 1e9);
    if ~isempty(X_fem)
        fprintf('  %10.2f nm        │\n', fw90_std_fem * 1e9);
    else
        fprintf('         N/A            │\n');
    end

    fprintf('  │ 旋转角       │  %10.4f °         │', rot_tmm);
    if ~isempty(X_fem)
        fprintf('  %10.4f °         │\n', rot_fem);
    else
        fprintf('         N/A            │\n');
    end

    fprintf('  │ 质心 X       │  %10.2f nm        │', cx_tmm * 1e9);
    if ~isempty(X_fem)
        fprintf('  %10.2f nm        │\n', cx_fem * 1e9);
    else
        fprintf('         N/A            │\n');
    end

    fprintf('  │ 质心 Y       │  %10.2f nm        │', cy_tmm * 1e9);
    if ~isempty(X_fem)
        fprintf('  %10.2f nm        │\n', cy_fem * 1e9);
    else
        fprintf('         N/A            │\n');
    end

    % 质心偏移差异
    if ~isempty(X_fem)
        dx_nm = (cx_tmm - cx_fem) * 1e9;
        dy_nm = (cy_tmm - cy_fem) * 1e9;
        fprintf('  │ 质心差Δ(X,Y) │   (%8.2f, %8.2f) nm                  │\n', dx_nm, dy_nm);
    end

    fprintf('  ├──────────────┴──────────────────────┴──────────────────────────┤\n');
    fprintf('  │ 各条纹 FWHM (nm):                                                │\n');
    fprintf('  │   TMM: ');
    for s = 1:length(fwhm_ind_tmm)
        if ~isnan(fwhm_ind_tmm(s))
            fprintf(' S%d: %.2f', s, fwhm_ind_tmm(s) * 1e9);
        end
    end
    fprintf('\n');
    if ~isempty(X_fem)
        fprintf('  │   FEM: ');
        for s = 1:length(fwhm_ind_fem)
            if ~isnan(fwhm_ind_fem(s))
                fprintf(' S%d: %.2f', s, fwhm_ind_fem(s) * 1e9);
            end
        end
        fprintf('\n');
    end
    fprintf('  │ 各条纹 FW90 (nm):                                                │\n');
    fprintf('  │   TMM: ');
    for s = 1:length(fw90_ind_tmm)
        if ~isnan(fw90_ind_tmm(s))
            fprintf(' S%d: %.2f', s, fw90_ind_tmm(s) * 1e9);
        end
    end
    fprintf('\n');
    if ~isempty(X_fem)
        fprintf('  │   FEM: ');
        for s = 1:length(fw90_ind_fem)
            if ~isnan(fw90_ind_fem(s))
                fprintf(' S%d: %.2f', s, fw90_ind_fem(s) * 1e9);
            end
        end
        fprintf('\n');
    end
    fprintf('  └─────────────────────────────────────────────────────────────────┘\n');

    % 保存结果
    results(f_idx).name      = base_name;
    results(f_idx).N         = N_particles;
    results(f_idx).offset_cx = offset_cx;
    results(f_idx).offset_cy = offset_cy;
    results(f_idx).tmm.cx = cx_tmm;
    results(f_idx).tmm.cy = cy_tmm;
    results(f_idx).tmm.fwhm = fwhm_tmm;
    results(f_idx).tmm.fwhm_std = fwhm_std_tmm;
    results(f_idx).tmm.fw90 = fw90_tmm;
    results(f_idx).tmm.fw90_std = fw90_std_tmm;
    results(f_idx).tmm.rotation = rot_tmm;
    results(f_idx).tmm.n_stripes = n_str_tmm;
    if ~isempty(X_fem)
        results(f_idx).fem.cx = cx_fem;
        results(f_idx).fem.cy = cy_fem;
        results(f_idx).fem.fwhm = fwhm_fem;
        results(f_idx).fem.fwhm_std = fwhm_std_fem;
        results(f_idx).fem.fw90 = fw90_fem;
        results(f_idx).fem.fw90_std = fw90_std_fem;
        results(f_idx).fem.rotation = rot_fem;
        results(f_idx).fem.n_stripes = n_str_fem;
    end
end

%% === 汇总对比表 ===
fprintf('\n\n');
fprintf('╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n');
fprintf('║                          汇总对比表 (TMM 全误差 vs FEM)                                      ║\n');
fprintf('╠════════════╤══════════╤══════════════════╤══════════════════╤════════════════════════════════╣\n');
fprintf('║   文件     │ 初始偏移 │  FWHM均值 (nm)   │   旋转角 (°)     │  质心 X, Y (nm) TMM / FEM     ║\n');
fprintf('╟────────────┼──────────┼──────────────────┼──────────────────┼────────────────────────────────╢\n');

for f_idx = 1:n_files
    if isempty(results(f_idx).name); continue; end

    name_short = file_entries{f_idx, 1};
    if length(name_short) > 10
        name_short = [name_short(1:7) '...'];
    end

    offset_str = sprintf('(%s%s)', ...
        format_offset(results(f_idx).offset_cx, 'x'), ...
        format_offset(results(f_idx).offset_cy, 'y'));

    t = results(f_idx).tmm;
    if isfield(results(f_idx), 'fem')
        f = results(f_idx).fem;

        fwhm_str = sprintf('%7.2f / %7.2f', t.fwhm * 1e9, f.fwhm * 1e9);
        rot_str  = sprintf('%7.3f / %7.3f', t.rotation, f.rotation);
        cx_str   = sprintf('%8.1f / %8.1f', t.cx * 1e9, f.cx * 1e9);
        cy_str   = sprintf('%8.1f / %8.1f', t.cy * 1e9, f.cy * 1e9);
    else
        fwhm_str = sprintf('%7.2f /    N/A  ', t.fwhm * 1e9);
        rot_str  = sprintf('%7.3f /    N/A  ', t.rotation);
        cx_str   = sprintf('%8.1f /    N/A  ', t.cx * 1e9);
        cy_str   = sprintf('%8.1f /    N/A  ', t.cy * 1e9);
    end

    fprintf('║ %-10s │ %-8s │ %s │ %s │ %s ║\n', name_short, offset_str, fwhm_str, rot_str, cx_str);
    fprintf('║            │          │                  │                  │ %s ║\n', cy_str);
end

fprintf('╚════════════════╧══════════╧══════════════════╧══════════════════╧════════════════════════════════╝\n');

%% === 可视化: TMM vs FEM 散点图 ===
figure('Name', '条纹 TMM(全误差) vs FEM 对比 (11-16)', 'Position', [50, 50, 1400, 800]);
for f_idx = 1:n_files
    if isempty(results(f_idx).name); continue; end
    subplot(2, 3, f_idx);

    % FEM 散点
    image_file = [file_entries{f_idx, 1} '_image.csv'];
    if exist(image_file, 'file')
        data_img = readmatrix(image_file);
        X_fem_plot = data_img(:, 1)';  % μm
        Y_fem_plot = data_img(:, 2)';
        scatter(X_fem_plot, Y_fem_plot, 3, 'r', '.', 'DisplayName', 'FEM');
        hold on;
    end

    % TMM 散点 (快速重跑)
    init_file = [file_entries{f_idx, 1} '_initial.csv'];
    data_init = readmatrix(init_file);
    Np = size(data_init, 1);
    V_p = zeros(3, 1, Np);
    V_p(1, 1, :) = data_init(:, 2)' * 1e-6 + file_entries{f_idx, 2};
    V_p(2, 1, :) = data_init(:, 3)' * 1e-6 + file_entries{f_idx, 3};
    V_p(3, 1, :) = 1;

    e_chg = 1.602e-19;
    p_p.U0 = data_init(:, 4)' / e_chg;
    vxy_p = sqrt(data_init(:, 5)'.^2 + data_init(:, 6)'.^2);
    p_p.alpha = atan2(vxy_p, data_init(:, 7)');
    p_p.phi = atan2(data_init(:, 6)', data_init(:, 5)');

    phys_p = calc_common_physics(p_p, V_p, sys);
    M_s = calc_M_spread(p_p, phys_p, Np);
    M_e0 = calc_M_EB0(p_p, phys_p, sys, err, Np);
    M_br = calc_M_Br1(p_p, phys_p, sys, err, Np);
    M_er = calc_M_EB_angle(p_p, phys_p, sys, err, Np);
    M_e = calc_M_ee(p_p, phys_p, sys, err, Np);
    M_1 = calc_M_E1(p_p, phys_p, sys, err, Np);
    M_b = calc_M_B1(p_p, phys_p, sys, err, Np);

    M_t = pagemtimes(M_e, pagemtimes(M_er, pagemtimes(M_b, ...
          pagemtimes(M_1, pagemtimes(M_br, pagemtimes(M_e0, M_s))))));
    Vf = pagemtimes(M_t, V_p);
    Pf = squeeze(Vf);
    % 减去输入位置偏移
    Pf(1, :) = Pf(1, :) - file_entries{f_idx, 2};
    Pf(2, :) = Pf(2, :) - file_entries{f_idx, 3};

    scatter(Pf(1, :) * 1e6, Pf(2, :) * 1e6, 3, 'b', '.', 'DisplayName', 'TMM');
    hold off;

    short_name = regexprep(file_entries{f_idx, 1}, '^\d+_system_stripe_', '');
    title(sprintf('File %d: %s', f_idx, short_name), 'Interpreter', 'none');
    xlabel('X (μm)'); ylabel('Y (μm)');
    legend('Location', 'best');
    axis equal; grid on;
end
sgtitle('TMM 全误差 (蓝) vs FEM (红) 条纹图案对比');

%% === 量化指标条形图 ===
figure('Name', '条纹量化指标 TMM vs FEM 对比', 'Position', [100, 100, 1400, 500]);

% 简短标签
labels = cell(n_files, 1);
for f_idx = 1:n_files
    name_parts = split(file_entries{f_idx, 1}, '（');
    if length(name_parts) >= 2
        labels{f_idx} = extractBefore(name_parts{2}, '）');
    else
        labels{f_idx} = sprintf('File%d', f_idx);
    end
end

% --- 子图1: FWHM 线宽对比 ---
subplot(1, 4, 1);
fwhm_t = zeros(n_files, 1); fwhm_f = zeros(n_files, 1);
for f_idx = 1:n_files
    if ~isempty(results(f_idx).name)
        fwhm_t(f_idx) = results(f_idx).tmm.fwhm * 1e9;
        if isfield(results(f_idx), 'fem')
            fwhm_f(f_idx) = results(f_idx).fem.fwhm * 1e9;
        else; fwhm_f(f_idx) = NaN; end
    else; fwhm_t(f_idx) = NaN; fwhm_f(f_idx) = NaN; end
end
bar([fwhm_t, fwhm_f]); set(gca, 'XTickLabel', labels, 'XTick', 1:n_files);
xtickangle(30); ylabel('FWHM (nm)'); title('FWHM 线宽');
legend('TMM', 'FEM', 'Location', 'best'); grid on;

% --- 子图2: 旋转角对比 ---
subplot(1, 4, 2);
rot_t = zeros(n_files, 1); rot_f = zeros(n_files, 1);
for f_idx = 1:n_files
    if ~isempty(results(f_idx).name)
        rot_t(f_idx) = results(f_idx).tmm.rotation;
        if isfield(results(f_idx), 'fem')
            rot_f(f_idx) = results(f_idx).fem.rotation;
        else; rot_f(f_idx) = NaN; end
    else; rot_t(f_idx) = NaN; rot_f(f_idx) = NaN; end
end
bar([rot_t, rot_f]); set(gca, 'XTickLabel', labels, 'XTick', 1:n_files);
xtickangle(30); ylabel('旋转角 (°)'); title('旋转角');
legend('TMM', 'FEM', 'Location', 'best'); grid on;

% --- 子图3: 质心 X 对比 ---
subplot(1, 4, 3);
cx_t = zeros(n_files, 1); cx_f = zeros(n_files, 1);
for f_idx = 1:n_files
    if ~isempty(results(f_idx).name)
        cx_t(f_idx) = results(f_idx).tmm.cx * 1e9;
        if isfield(results(f_idx), 'fem')
            cx_f(f_idx) = results(f_idx).fem.cx * 1e9;
        else; cx_f(f_idx) = NaN; end
    else; cx_t(f_idx) = NaN; cx_f(f_idx) = NaN; end
end
bar([cx_t, cx_f]); set(gca, 'XTickLabel', labels, 'XTick', 1:n_files);
xtickangle(30); ylabel('质心 X (nm)'); title('质心 X 位置');
legend('TMM', 'FEM', 'Location', 'best'); grid on;

% --- 子图4: 质心 Y 对比 ---
subplot(1, 4, 4);
cy_t = zeros(n_files, 1); cy_f = zeros(n_files, 1);
for f_idx = 1:n_files
    if ~isempty(results(f_idx).name)
        cy_t(f_idx) = results(f_idx).tmm.cy * 1e9;
        if isfield(results(f_idx), 'fem')
            cy_f(f_idx) = results(f_idx).fem.cy * 1e9;
        else; cy_f(f_idx) = NaN; end
    else; cy_t(f_idx) = NaN; cy_f(f_idx) = NaN; end
end
bar([cy_t, cy_f]); set(gca, 'XTickLabel', labels, 'XTick', 1:n_files);
xtickangle(30); ylabel('质心 Y (nm)'); title('质心 Y 位置');
legend('TMM', 'FEM', 'Location', 'best'); grid on;

sgtitle('条纹后处理量化指标汇总: TMM (全误差) vs FEM (Files 11-16)');

fprintf('\n分析完成！\n');

%% ========================================================================
%  辅助函数: 格式化偏移量显示
% ========================================================================
function s = format_offset(val, axis_label)
    if val == 0
        s = '';
    elseif abs(val) >= 1e-3
        s = sprintf('%s=%.4gmm,', axis_label, val * 1e3);
    elseif abs(val) >= 1e-6
        s = sprintf('%s=%.4gμm,', axis_label, val * 1e6);
    else
        s = sprintf('%s=%.4gnm,', axis_label, val * 1e9);
    end
end

%% ========================================================================
%  子函数：条纹图案量化分析
%  算法：
%    1. X方向直方图投影识别条纹位置 (自适应峰值检测)
%    2. 逐条局部PCA计算旋转角
%    3. FWHM = 2.355 * std(X) 逐条计算 (假设类高斯分布)
%    4. FW90 = 包含90%电子的宽度 (5%-95%分位数)
% ========================================================================
function [cx, cy, fwhm_mean, rotation_angle_deg, fwhm_std, fwhm_list, ...
          fw90_mean, fw90_std, fw90_list, n_stripes] = ...
        analyze_stripe_pattern(X, Y)
    X = X(:)';  Y = Y(:)';  N_total = length(X);

    % 1. 整体质心
    cx = mean(X);  cy = mean(Y);
    Xc = X - cx;   Yc = Y - cy;

    % 2. X方向直方图投影 → 找条纹峰值
    n_bins = max(120, min(300, round(N_total / 120)));
    [counts, edges] = histcounts(Xc, n_bins);
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

    smooth_width = max(2, round(n_bins / 25));
    counts_sm = smoothdata(counts, 'gaussian', smooth_width);

    min_peak_dist = max(4, round(n_bins / 12));
    min_peak_height = max(counts_sm) * 0.04;
    [pks, locs] = findpeaks(counts_sm, 'MinPeakDistance', min_peak_dist, ...
        'MinPeakHeight', min_peak_height);

    % 如果太多峰 → 只保留最强7个并按间距过滤
    if length(locs) > 7
        [~, sort_idx] = sort(pks, 'descend');
        locs = locs(sort_idx(1:7));
    end
    if length(locs) > 1
        peak_positions = sort(bin_centers(locs));
        avg_spacing = (peak_positions(end) - peak_positions(1)) / (length(locs) - 1);
        min_allowed_dist = avg_spacing * 0.4;
        [~, sort_idx] = sort(pks, 'descend');
        selected = []; selected_pos = [];
        for ii = 1:length(sort_idx)
            cand_pos = bin_centers(locs(sort_idx(ii)));
            if isempty(selected_pos) || all(abs(selected_pos - cand_pos) > min_allowed_dist)
                selected = [selected, locs(sort_idx(ii))];
                selected_pos = [selected_pos, cand_pos];
            end
        end
        locs = selected;
    end
    n_stripes = length(locs);

    if n_stripes == 0
        [~, locs] = findpeaks(counts_sm, 'MinPeakDistance', 3, ...
            'MinPeakHeight', max(counts_sm) * 0.015);
        n_stripes = length(locs);
    end
    if n_stripes == 0
        fwhm_mean = NaN; fwhm_std = NaN; fwhm_list = NaN;
        fw90_mean = NaN; fw90_std = NaN; fw90_list = NaN;
        rotation_angle_deg = NaN; return;
    end

    % 3. 条纹半间距 (用于提取每条条纹的数据)
    if n_stripes > 1
        sorted_peaks = sort(bin_centers(locs));
        half_gap = min(diff(sorted_peaks)) * 0.45;
    else
        half_gap = (max(Xc) - min(Xc)) / 3;
    end

    % 4. 逐条纹分析: FWHM, FW90, 局部PCA求旋转角
    fwhm_list = zeros(1, n_stripes);
    fw90_list = zeros(1, n_stripes);
    rot_list  = zeros(1, n_stripes);

    for s = 1:n_stripes
        peak_x = bin_centers(locs(s));
        mask = abs(Xc - peak_x) < half_gap;
        Xs = Xc(mask);  Ys = Yc(mask);
        n_pts = length(Xs);

        if n_pts < 30
            fwhm_list(s) = NaN;  fw90_list(s) = NaN;  rot_list(s) = NaN;
            continue;
        end

        % --- FWHM = 2.355 * sigma (单条纹 X 分布标准差) ---
        fwhm_list(s) = 2.355 * std(Xs);

        % --- FW90 = 包含90%电子的宽度 (5%-95%分位数) ---
        Xs_sorted = sort(Xs);
        idx5  = max(1, round(0.05 * n_pts));
        idx95 = min(n_pts, round(0.95 * n_pts));
        fw90_list(s) = Xs_sorted(idx95) - Xs_sorted(idx5);

        % --- 局部PCA → 旋转角 ---
        Xs_loc = Xs - mean(Xs);  Ys_loc = Ys - mean(Ys);
        cov_loc = cov([Xs_loc; Ys_loc]');
        [V_loc, D_loc] = eig(cov_loc);
        [~, idx_sort] = sort(diag(D_loc), 'descend');
        V_loc = V_loc(:, idx_sort);
        local_angle = atan2d(V_loc(2,1), V_loc(1,1));
        rot_s = 90 - local_angle;
        if rot_s > 90; rot_s = rot_s - 180;
        elseif rot_s < -90; rot_s = rot_s + 180; end
        rot_list(s) = rot_s;
    end

    fwhm_valid = fwhm_list(~isnan(fwhm_list));
    if isempty(fwhm_valid)
        fwhm_mean = NaN; fwhm_std = NaN;
    else
        fwhm_mean = mean(fwhm_valid); fwhm_std = std(fwhm_valid);
    end

    fw90_valid = fw90_list(~isnan(fw90_list));
    if isempty(fw90_valid)
        fw90_mean = NaN; fw90_std = NaN;
    else
        fw90_mean = mean(fw90_valid); fw90_std = std(fw90_valid);
    end

    rot_valid = rot_list(~isnan(rot_list));
    if isempty(rot_valid)
        rotation_angle_deg = NaN;
    else
        rotation_angle_deg = mean(rot_valid);
    end
end
