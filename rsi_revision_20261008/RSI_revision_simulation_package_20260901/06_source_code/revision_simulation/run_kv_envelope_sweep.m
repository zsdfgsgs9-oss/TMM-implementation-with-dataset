function run_kv_envelope_sweep()
%RUN_KV_ENVELOPE_SWEEP Compare K-V and self-consistent TMM space charge.
% The round-beam edge convention is R=2*sigma_x and
% epsilon_edge=4*epsilon_rms. Results are reported after subtracting each
% model's zero-current baseline.

    this_dir = fileparts(mfilename('fullpath'));
    source_dir = fileparts(this_dir);
    package_dir = fileparts(source_dir);
    root_matlab_dir = fullfile(source_dir, 'root_matlab');
    output_dir = fullfile(this_dir, 'results');
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    addpath(root_matlab_dir);
    addpath(this_dir);

    input_file = fullfile(package_dir, '07_source_inputs', ...
        'intrinsic_pair', '01_intrinsic_error_initial.csv');
    [V_init, p] = load_comsol_initial_csv(input_file);
    x = squeeze(V_init(1, 1, :)).';
    y = squeeze(V_init(2, 1, :)).';
    x = x - mean(x);
    y = y - mean(y);
    xp = tan(p.alpha) .* cos(p.phi);
    yp = tan(p.alpha) .* sin(p.phi);
    xp = xp - mean(xp);
    yp = yp - mean(yp);

    sigma_x = sqrt(mean(x.^2));
    sigma_y = sqrt(mean(y.^2));
    eps_x = sqrt(max(mean(x.^2) * mean(xp.^2) - mean(x .* xp)^2, 0));
    eps_y = sqrt(max(mean(y.^2) * mean(yp.^2) - mean(y .* yp)^2, 0));
    R0 = sigma_x + sigma_y;
    Rprime0 = mean(x .* xp) / sigma_x + mean(y .* yp) / sigma_y;
    epsilon_edge_geom = 2 * (eps_x + eps_y);

    e = 1.602176634e-19;
    m = 9.1093837015e-31;
    eps0 = 8.8541878128e-12;
    c = 299792458;
    U = 12.67e3;
    L = 1e-3;
    B = 1.2;
    U0_axial = mean(p.U0 .* cos(p.alpha).^2);
    v0 = sqrt(2 * e * U0_axial / m);
    beta0 = v0 / c;
    epsilon_edge_norm = beta0 * epsilon_edge_geom;
    I_A = 4 * pi * eps0 * m * c^3 / e;

    sys.L = L;
    sys.U = U;
    sys.B_vec_0 = [0; 0; B];
    sys.E_vec_0 = [0; 0; U / L];
    sys.physics_model = 'nonrelativistic';
    phys = calc_common_physics(p, V_init, sys);

    zero_ode = @(z, state) envelope_rhs(z, state, 0, U0_axial, U, L, ...
        B, epsilon_edge_norm, I_A, e, m, c);
    options = odeset('RelTol', 1e-8, 'AbsTol', [1e-15, 1e-10], ...
        'MaxStep', L / 2000);
    [z_zero, state_zero] = ode15s(zero_ode, [0, L], ...
        [R0; Rprime0], options);
    W_zero = U0_axial + U .* z_zero ./ L;
    v_zero = sqrt(2 * e .* W_zero ./ m);
    a_z = e * U / (m * L);
    time_zero = (v_zero - v0) ./ a_z;
    kv_reference.time = time_zero;
    kv_reference.v_z = v_zero;
    kv_reference.R_ideal = state_zero(:, 1);
    kv_reference.omega_larmor = (e * B / (2 * m)) * ones(size(z_zero));
    kv_reference.gamma = ones(size(z_zero));

    currents = [0, 1e-12, 1e-11, 1e-10, 5e-10, 1e-9, ...
        2e-9, 1e-8, 1e-7];
    rows = table();
    kv_trajectories = cell(numel(currents), 1);
    tmm_native_trajectories = cell(numel(currents), 1);
    tmm_matched_trajectories = cell(numel(currents), 1);

    for i = 1:numel(currents)
        I = currents(i);
        ode = @(z, state) envelope_rhs(z, state, I, U0_axial, U, L, ...
            B, epsilon_edge_norm, I_A, e, m, c);
        [z, state] = ode15s(ode, [0, L], [R0; Rprime0], options);
        kv_trajectories{i} = [z, state];
        R_zero_at_z = interp1(z_zero, state_zero(:, 1), z, 'pchip');
        kv_delta_profile = state(:, 1) - R_zero_at_z;

        err_native.I_beam = I;
        err_native.R0 = R0;
        err_native.ee_method = 'self_consistent';
        err_native.ee_steps = 4000;
        err_native.ee_quiet = true;
        err_native.ee_emittance_coefficient = (epsilon_edge_norm * c)^2;
        [~, tmm_native] = calc_M_ee(p, phys, sys, err_native, 1);
        tmm_native_trajectories{i} = tmm_native;

        err_matched = err_native;
        err_matched.ee_reference_profile = kv_reference;
        [~, tmm_matched] = calc_M_ee(p, phys, sys, err_matched, 1);
        tmm_matched_trajectories{i} = tmm_matched;

        one = table(I, I / (pi * R0^2), R0 * 1e9, ...
            epsilon_edge_geom * 1e9, epsilon_edge_norm * 1e9, ...
            state(end, 1) * 1e9, 2 * state(end, 1) * 1e9, ...
            min(state(:, 1)) * 1e9, ...
            tmm_native.R_ideal(end) * 1e9, ...
            tmm_native.R_total(end) * 1e9, ...
            2 * tmm_native.legacy_delta_R_final * 1e9, ...
            2 * tmm_native.delta_R(end) * 1e9, tmm_native.scale, ...
            2 * tmm_matched.delta_R(end) * 1e9, ...
            2 * max(kv_delta_profile) * 1e9, ...
            2 * max(tmm_native.delta_R) * 1e9, ...
            2 * max(tmm_matched.delta_R) * 1e9, ...
            'VariableNames', {'I_A', 'J_initial_A_m2', 'R0_nm', ...
            'epsilon_edge_geom_nm_rad', 'epsilon_edge_norm_nm_rad', ...
            'KV_R_final_nm', 'KV_diameter_final_nm', 'KV_R_min_nm', ...
            'TMM_R_ideal_final_nm', 'TMM_R_total_final_nm', ...
            'legacy_delta_diameter_nm', 'native_delta_diameter_nm', ...
            'TMM_native_scale', 'matched_delta_diameter_nm', ...
            'KV_peak_delta_diameter_nm', 'native_peak_delta_diameter_nm', ...
            'matched_peak_delta_diameter_nm'});
        rows = [rows; one]; %#ok<AGROW>
    end

    rows.KV_delta_diameter_nm = rows.KV_diameter_final_nm - ...
        rows.KV_diameter_final_nm(1);
    rows.legacy_abs_error_nm = abs(rows.legacy_delta_diameter_nm - ...
        rows.KV_delta_diameter_nm);
    rows.native_abs_error_nm = abs(rows.native_delta_diameter_nm - ...
        rows.KV_delta_diameter_nm);
    rows.matched_abs_error_nm = abs(rows.matched_delta_diameter_nm - ...
        rows.KV_delta_diameter_nm);
    rows.native_signed_error_nm = rows.native_delta_diameter_nm - ...
        rows.KV_delta_diameter_nm;
    rows.native_peak_abs_error_nm = abs(rows.native_peak_delta_diameter_nm - ...
        rows.KV_peak_delta_diameter_nm);
    rows.matched_peak_abs_error_nm = abs(rows.matched_peak_delta_diameter_nm - ...
        rows.KV_peak_delta_diameter_nm);

    convergence = table();
    convergence_currents = [1e-9, 1e-7];
    convergence_steps = [250, 500, 1000, 2000, 4000, 8000];
    for i = 1:numel(convergence_currents)
        for j = 1:numel(convergence_steps)
            test_err = err_native;
            test_err.I_beam = convergence_currents(i);
            test_err.ee_steps = convergence_steps(j);
            [~, native_test] = calc_M_ee(p, phys, sys, test_err, 1);
            test_err.ee_reference_profile = kv_reference;
            [~, matched_test] = calc_M_ee(p, phys, sys, test_err, 1);
            one = table(convergence_currents(i), convergence_steps(j), ...
                2 * native_test.delta_R(end) * 1e9, ...
                2 * max(native_test.delta_R) * 1e9, ...
                2 * matched_test.delta_R(end) * 1e9, ...
                2 * max(matched_test.delta_R) * 1e9, ...
                'VariableNames', {'I_A', 'steps', ...
                'native_delta_diameter_nm', 'native_peak_delta_diameter_nm', ...
                'matched_delta_diameter_nm', ...
                'matched_peak_delta_diameter_nm'});
            convergence = [convergence; one]; %#ok<AGROW>
        end
    end

    output_csv = fullfile(output_dir, 'space_charge_kv_comparison.csv');
    writetable(rows, output_csv);
    writetable(convergence, fullfile(output_dir, ...
        'space_charge_kv_convergence.csv'));
    save(fullfile(output_dir, 'space_charge_kv_comparison.mat'), ...
        'rows', 'kv_trajectories', 'tmm_native_trajectories', ...
        'tmm_matched_trajectories', 'U0_axial', 'Rprime0', 'R0', ...
        'kv_reference', 'convergence');

    fprintf('\nK-V versus TMM space-charge comparison\n');
    fprintf('R0 = %.6f nm, edge emittance = %.6f nm rad\n', ...
        R0 * 1e9, epsilon_edge_geom * 1e9);
    disp(rows(:, {'I_A', 'J_initial_A_m2', 'KV_delta_diameter_nm', ...
        'legacy_delta_diameter_nm', 'native_delta_diameter_nm', ...
        'matched_delta_diameter_nm', 'KV_peak_delta_diameter_nm', ...
        'native_peak_delta_diameter_nm', 'matched_peak_delta_diameter_nm'}));
    fprintf('Saved %s\n', output_csv);
end

function derivative = envelope_rhs(z, state, I, U0, U, L, B, ...
        eps_norm, I_A, e, m, c)
    R = max(state(1), 1e-15);
    W = max(U0 + U * z / L, realmin);
    v = sqrt(2 * e * W / m);
    beta = v / c;
    pprime_over_p = 0.5 * (U / L) / W;
    k = e * B / (2 * m * v);
    eps_geom = eps_norm / beta;
    perveance = 2 * I / (I_A * beta^3);

    % For R=2*sigma in a round beam, the Sacherer term is K/R.
    derivative = [state(2); ...
        -pprime_over_p * state(2) - k^2 * R + ...
        eps_geom^2 / R^3 + perveance / R];
end
