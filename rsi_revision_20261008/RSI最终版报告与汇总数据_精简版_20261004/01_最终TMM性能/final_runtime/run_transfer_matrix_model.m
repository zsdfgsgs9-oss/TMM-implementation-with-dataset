function result = run_transfer_matrix_model(p, V_init, sys, err, num_electrons)
% RUN_TRANSFER_MATRIX_MODEL Run one full transfer-matrix simulation.

    if nargin < 5 || isempty(num_electrons)
        num_electrons = size(V_init, 3);
    end

    phys = calc_common_physics(p, V_init, sys);
    fprintf('正在构建转移矩阵 (%s)...\n', phys.model);

    use_unified_b_gradient = ~phys.is_relativistic && ...
        isfield(sys, 'b_gradient') && isfield(sys.b_gradient, 'enabled') && ...
        sys.b_gradient.enabled;

    if use_unified_b_gradient
        [M_transport, b_gradient_diagnostics] = ...
            calc_M_B_gradient(p, phys, sys, err, num_electrons);
        M_spread = M_transport;
        M_Br1 = repmat(eye(3), 1, 1, num_electrons);
        M_B1 = repmat(eye(3), 1, 1, num_electrons);
        M_EBrho = repmat(eye(3), 1, 1, num_electrons);
    else
        M_spread = calc_M_spread(p, phys, num_electrons);
        M_Br1 = calc_M_Br1(p, phys, sys, err, num_electrons);
        M_B1 = calc_M_B1(p, phys, sys, err, num_electrons);
        M_EBrho = calc_M_EB_angle(p, phys, sys, err, num_electrons);
        b_gradient_diagnostics = struct('enabled', false);
    end
    M_EB0 = calc_M_EB0(p, phys, sys, err, num_electrons);
    M_ee = calc_M_ee(p, phys, sys, err, num_electrons);
    M_E1 = calc_M_E1(p, phys, sys, err, num_electrons);

    fprintf('正在进行批量矩阵运算 (%s)...\n', phys.model);
    M_temp1 = pagemtimes(M_EB0, M_spread);
    M_temp2 = pagemtimes(M_Br1, M_temp1);
    M_temp3 = pagemtimes(M_E1, M_temp2);
    M_temp4 = pagemtimes(M_B1, M_temp3);
    M_temp5 = pagemtimes(M_EBrho, M_temp4);
    M_total = pagemtimes(M_ee, M_temp5);

    V_final = pagemtimes(M_total, V_init);
    Pos_final = squeeze(V_final);

    result.phys = phys;
    result.V_final = V_final;
    result.X_final = Pos_final(1, :);
    result.Y_final = Pos_final(2, :);
    result.metrics = summarize_beam_metrics(result.X_final, result.Y_final, ...
        squeeze(V_init(1, 1, :)).', squeeze(V_init(2, 1, :)).');
    result.diagnostics.phase = phys.nonrel.phase_diagnostics;
    result.diagnostics.b_gradient = b_gradient_diagnostics;
end
