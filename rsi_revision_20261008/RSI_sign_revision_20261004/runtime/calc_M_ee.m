function [M, diagnostics] = calc_M_ee(p, phys, ~, err, N)
%CALC_M_EE Build the 3-by-3 affine space-charge operator.
% The default model advances one shared, self-consistent beam envelope. The
% envelope displacement and velocity are internal integration variables;
% the particle state remains (x,y,1). Explicit 'sum' and 'gauss' methods
% retain the previous edge-displacement model for regression comparisons.

    eps0 = 8.8541878128e-12;
    e = 1.602176634e-19;
    m = 9.1093837015e-31;

    if nargin < 5 || isempty(N)
        N = numel(phys.t);
    end
    if ~isfield(err, 'I_beam') || ~isfield(err, 'R0')
        error('calc_M_ee requires err.I_beam and err.R0.');
    end

    method = 'self_consistent';
    if isfield(err, 'ee_method') && ~isempty(err.ee_method)
        method = lower(char(string(err.ee_method)));
    end

    switch method
        case {'self_consistent', 'self-consistent', 'selfconsistent', 'sc'}
            timer_id = tic;
            [scale, diagnostics] = solve_shared_envelope(p, phys, err, ...
                eps0, e, m);

            % Uniform round-beam self field is linear inside the beam. Its
            % action is therefore a dilation about the beam centroid.
            x_center = finite_mean(phys.X, 0);
            y_center = finite_mean(phys.Y, 0);
            M = repmat(eye(3), 1, 1, N);
            M(1, 1, :) = scale;
            M(2, 2, :) = scale;
            M(1, 3, :) = (1 - scale) * x_center;
            M(2, 3, :) = (1 - scale) * y_center;

            diagnostics.method = 'self_consistent';
            diagnostics.scale = scale;
            diagnostics.centroid = [x_center, y_center];
            diagnostics.elapsed_seconds = toc(timer_id);
            if ~isfield(err, 'ee_quiet') || ~err.ee_quiet
                fprintf(['Space charge (self-consistent envelope, %d steps) ', ...
                    'completed in %.4f s\n'], diagnostics.n_steps, ...
                    diagnostics.elapsed_seconds);
            end

        case {'sum', 'gauss', 'legacy_sum', 'legacy_gauss'}
            timer_id = tic;
            [delta_ee, diagnostics] = solve_legacy_model(phys, err, ...
                method, eps0, e, m);
            M = repmat(eye(3), 1, 1, N);

            if isfield(p, 'epsilon_ee')
                epsilon = p.epsilon_ee;
            else
                epsilon = 2 * pi * rand(1, N);
            end
            M(1, 3, :) = delta_ee .* cos(epsilon);
            M(2, 3, :) = delta_ee .* sin(epsilon);

            diagnostics.method = method;
            diagnostics.elapsed_seconds = toc(timer_id);
            if ~isfield(err, 'ee_quiet') || ~err.ee_quiet
                fprintf('Space charge (%s) completed in %.4f s\n', ...
                    method, diagnostics.elapsed_seconds);
            end

        otherwise
            error(['Unknown err.ee_method. Use ''self_consistent'', ', ...
                '''sum'', or ''gauss''.']);
    end
end

function [scale, out] = solve_shared_envelope(p, phys, err, eps0, e, m)
% Advance q'' + omega_L^2 q = a_sc(t,q) with a causal predictor-corrector.

    n_steps = 4000;
    if isfield(err, 'ee_steps')
        n_steps = max(20, round(err.ee_steps));
    elseif isfield(err, 'ee_rel_steps')
        n_steps = max(20, round(err.ee_rel_steps));
    end

    R0 = err.R0;
    I = err.I_beam;
    radius_floor = max(1e-15, abs(R0) * 1e-6);
    if isfield(err, 'ee_radius_floor')
        radius_floor = max(realmin, err.ee_radius_floor);
    end

    [time, v_z, gamma, omega_l, R_ideal, reference] = ...
        build_reference_profiles(phys, err, n_steps, radius_floor, e, m);
    dt = diff(time);
    q = zeros(size(time));
    q_dot = zeros(size(time));

    coefficient = e * I / (2 * pi * eps0 * m);
    [emittance_coefficient, emittance] = ...
        infer_emittance_coefficient(p, phys);
    if isfield(err, 'ee_emittance_coefficient')
        emittance_coefficient = max(0, err.ee_emittance_coefficient);
        emittance.source = 'err.ee_emittance_coefficient';
    end
    for k = 1:n_steps
        omega_mid = 0.5 * (omega_l(k) + omega_l(k + 1));
        a0 = space_charge_acceleration(coefficient, v_z(k), gamma(k), ...
            R_ideal(k) + q(k), radius_floor) + ...
            emittance_feedback(emittance_coefficient, R_ideal(k), ...
            q(k), radius_floor);

        [q_predict, ~] = oscillator_step(q(k), q_dot(k), ...
            a0, omega_mid, dt(k));
        a1 = space_charge_acceleration(coefficient, v_z(k + 1), ...
            gamma(k + 1), R_ideal(k + 1) + q_predict, radius_floor) + ...
            emittance_feedback(emittance_coefficient, R_ideal(k + 1), ...
            q_predict, radius_floor);

        [q(k + 1), q_dot(k + 1)] = oscillator_step(q(k), q_dot(k), ...
            0.5 * (a0 + a1), omega_mid, dt(k));

        % A negative envelope is nonphysical and indicates that the model
        % has been driven beyond its self-similar validity range.
        min_q = radius_floor - R_ideal(k + 1);
        if q(k + 1) < min_q
            q(k + 1) = min_q;
            q_dot(k + 1) = max(q_dot(k + 1), 0);
        end
    end

    R_total = max(R_ideal + q, radius_floor);
    scale = R_total(end) / max(R_ideal(end), radius_floor);

    % Previous paper model, evaluated on the same shared reference profile.
    % This isolates the effects of feedback and magnetic response.
    a_legacy = coefficient ./ ...
        (max(v_z, realmin) .* max(gamma, 1).^3 .* max(R_ideal, radius_floor));
    legacy_delta = trapz(time, (time(end) - time) .* a_legacy);

    out.time = time;
    out.R_ideal = R_ideal;
    out.delta_R = q;
    out.delta_R_dot = q_dot;
    out.R_total = R_total;
    out.v_z = v_z;
    out.gamma = gamma;
    out.omega_larmor = omega_l;
    out.legacy_delta_R_final = legacy_delta;
    out.emittance_coefficient = emittance_coefficient;
    out.emittance = emittance;
    out.n_steps = n_steps;
    out.reference = reference;
    out.hit_radius_floor = any(R_ideal + q <= radius_floor);
end

function [coefficient, out] = infer_emittance_coefficient(p, phys)
% Infer the round-beam edge emittance with bounded temporary memory.

    coefficient = 0;
    out.edge_geometric = 0;
    out.normalized = 0;
    out.source = 'unavailable';
    required = {'alpha', 'phi'};
    if ~all(isfield(p, required)) || ~isfield(phys, 'X') || ...
            ~isfield(phys, 'Y') || ~isfield(phys, 'v1')
        return
    end

    n = min([numel(p.alpha), numel(p.phi), numel(phys.X), ...
        numel(phys.Y), numel(phys.v1)]);
    if n < 2
        return
    end

    sums = zeros(1, 10);
    count = 0;
    block_size = 1e6;
    for first = 1:block_size:n
        last = min(first + block_size - 1, n);
        index = first:last;
        x = double(phys.X(index));
        y = double(phys.Y(index));
        alpha = double(p.alpha(index));
        phi = double(p.phi(index));
        xp = tan(alpha) .* cos(phi);
        yp = tan(alpha) .* sin(phi);
        valid = isfinite(x) & isfinite(y) & isfinite(xp) & isfinite(yp);
        x = x(valid); y = y(valid); xp = xp(valid); yp = yp(valid);
        if isempty(x), continue, end
        sums = sums + [sum(x), sum(y), sum(xp), sum(yp), ...
            sum(x.^2), sum(y.^2), sum(xp.^2), sum(yp.^2), ...
            sum(x .* xp), sum(y .* yp)];
        count = count + numel(x);
    end
    if count < 2
        return
    end

    means = sums(1:4) / count;
    variance_x = max(sums(5) / count - means(1)^2, 0);
    variance_y = max(sums(6) / count - means(2)^2, 0);
    variance_xp = max(sums(7) / count - means(3)^2, 0);
    variance_yp = max(sums(8) / count - means(4)^2, 0);
    covariance_x = sums(9) / count - means(1) * means(3);
    covariance_y = sums(10) / count - means(2) * means(4);
    epsilon_x = sqrt(max(variance_x * variance_xp - covariance_x^2, 0));
    epsilon_y = sqrt(max(variance_y * variance_yp - covariance_y^2, 0));
    edge_geometric = 2 * (epsilon_x + epsilon_y);

    c = 299792458;
    if isfield(phys, 'is_relativistic') && phys.is_relativistic && ...
            isfield(phys, 'rel') && isfield(phys.rel, 'p0')
        beta_gamma = finite_mean(phys.rel.p0, 0) / ...
            (9.1093837015e-31 * c);
    else
        beta_gamma = positive_rms(phys.v1, 0) / c;
    end
    normalized = beta_gamma * edge_geometric;
    coefficient = (normalized * c)^2;
    out.edge_geometric = edge_geometric;
    out.normalized = normalized;
    out.source = 'particle_second_moments';
end

function acceleration = emittance_feedback(coefficient, R_ideal, q, radius_floor)
% The prescribed ideal envelope already contains its zero-current emittance
% pressure. Only the change caused by the space-charge dilation is added.

    if coefficient == 0
        acceleration = 0;
        return
    end
    R_base = max(R_ideal, radius_floor);
    R_total = max(R_ideal + q, radius_floor);
    acceleration = coefficient * (1 / R_total^3 - 1 / R_base^3);
end

function acceleration = space_charge_acceleration(coefficient, velocity, ...
        gamma, radius, radius_floor)
% gamma^-3 includes electric/magnetic cancellation and transverse inertia.

    acceleration = coefficient / (max(abs(velocity), realmin) * ...
        max(gamma, 1)^3 * max(radius, radius_floor));
end

function [q_next, u_next] = oscillator_step(q, u, acceleration, omega, dt)
% Exact update for a harmonic oscillator with constant forcing over one step.

    theta = omega * dt;
    if abs(theta) < 1e-7
        q_next = q + u * dt + 0.5 * acceleration * dt^2;
        u_next = u + acceleration * dt;
        return
    end

    cosine = cos(theta);
    sine = sin(theta);
    q_next = q * cosine + u * sine / omega + ...
        acceleration * (1 - cosine) / omega^2;
    u_next = -omega * q * sine + u * cosine + ...
        acceleration * sine / omega;
end

function [time, v_z, gamma, omega_l, R_ideal, reference] = ...
        build_reference_profiles(phys, err, n_steps, radius_floor, e, m)
% Reduce particle-dependent kinematics to one current-weighted beam profile.

    c = 299792458;
    if isfield(err, 'ee_reference_profile')
        profile = err.ee_reference_profile;
        required = {'time', 'v_z', 'R_ideal', 'omega_larmor'};
        for k = 1:numel(required)
            if ~isfield(profile, required{k})
                error('err.ee_reference_profile is missing field %s.', required{k});
            end
        end
        source_time = profile.time(:).';
        if numel(source_time) < 2 || any(diff(source_time) <= 0)
            error('Reference-profile time must be strictly increasing.');
        end
        time = linspace(source_time(1), source_time(end), n_steps + 1);
        v_z = interp1(source_time, profile.v_z(:).', time, 'pchip');
        R_ideal = interp1(source_time, profile.R_ideal(:).', time, 'pchip');
        omega_l = interp1(source_time, profile.omega_larmor(:).', ...
            time, 'pchip');
        if isfield(profile, 'gamma')
            gamma = interp1(source_time, profile.gamma(:).', time, 'pchip');
        else
            gamma = ones(size(time));
        end
        R_ideal = max(R_ideal, radius_floor);
        reference.L = finite_mean(phys.L_local, NaN);
        reference.U = finite_mean(phys.U_local, NaN);
        reference.t_final = time(end);
        reference.omega_c = 2 * finite_mean(omega_l, 0);
        reference.R_larmor = NaN;
        reference.relativistic = any(gamma > 1 + 1e-12);
        reference.source = 'provided_profile';
        return
    end

    L_ref = finite_mean(phys.L_local, 1e-3);
    U_ref = finite_mean(phys.U_local, 0);
    E_ref = U_ref / max(L_ref, realmin);
    omega_c_ref = finite_abs_mean(phys.omega, 0);
    R_larmor_edge = finite_max_abs(phys.R, 0);
    R_larmor_rms = finite_rms(phys.R, 0);
    envelope_convention = 'rms';
    if isfield(err, 'ee_envelope_convention')
        envelope_convention = lower(char(string(err.ee_envelope_convention)));
    end

    if isfield(err, 'ee_reference_larmor_radius')
        R_larmor_edge = max(0, err.ee_reference_larmor_radius);
        R_larmor_rms = R_larmor_edge;
    end
    if isfield(err, 'ee_reference_omega_c')
        omega_c_ref = abs(err.ee_reference_omega_c);
    end

    is_relativistic = isfield(phys, 'is_relativistic') && ...
        phys.is_relativistic;
    if is_relativistic
        pz0 = finite_mean(phys.rel.p_parallel0, 0);
        a_rel = finite_mean(phys.rel.a, 1);
        gamma_final = finite_mean(phys.rel.gamma_L, 1);
        pz_final = m * c * sqrt(max(gamma_final^2 - a_rel^2, 0));
        t_final = (pz_final - pz0) / max(e * E_ref, realmin);
        time = linspace(0, max(t_final, 0), n_steps + 1);
        pz = pz0 + e * E_ref .* time;
        gamma = sqrt(a_rel^2 + (pz ./ (m * c)).^2);
        v_z = pz ./ (gamma * m);
        B_ref = finite_mean(abs(phys.B_vec(3, :)), 0);
        omega_l = e * B_ref ./ (2 * gamma * m);
        phase = cumtrapz(time, 2 * omega_l);
    else
        % The reference particle carries the ensemble mean axial energy.
        v1_ref = max(positive_rms(phys.v1, realmin), realmin);
        a_z_ref = e * U_ref / (m * max(L_ref, realmin));
        if a_z_ref > 0
            t_final = (-v1_ref + sqrt(v1_ref^2 + 2 * a_z_ref * L_ref)) / a_z_ref;
        else
            t_final = L_ref / max(v1_ref, realmin);
        end
        time = linspace(0, max(t_final, 0), n_steps + 1);
        v_z = v1_ref + a_z_ref .* time;
        gamma = ones(size(time));
        omega_l = 0.5 * omega_c_ref * ones(size(time));
        phase = omega_c_ref .* time;
    end

    if strcmp(envelope_convention, 'edge')
        % Original manuscript upper-edge construction.
        R_ideal = err.R0 + 2 * R_larmor_edge .* sin(0.5 .* phase);
    else
        % For an axisymmetric distribution with random azimuth, the edge
        % convention R=2*sigma gives this second-moment envelope.
        R_ideal = sqrt(err.R0^2 + ...
            8 * R_larmor_rms^2 .* sin(0.5 .* phase).^2);
    end
    R_ideal = max(R_ideal, radius_floor);

    reference.L = L_ref;
    reference.U = U_ref;
    reference.t_final = time(end);
    reference.omega_c = omega_c_ref;
    reference.R_larmor_edge = R_larmor_edge;
    reference.R_larmor_rms = R_larmor_rms;
    reference.envelope_convention = envelope_convention;
    reference.relativistic = is_relativistic;
    reference.source = 'native_tmm';
end

function [delta_ee, out] = solve_legacy_model(phys, err, method, eps0, e, m)
% Original non-self-consistent edge translation retained for comparison.

    I = err.I_beam;
    R0 = err.R0;
    t_flight = phys.t;
    a_z = (e .* phys.U_local) ./ (m .* phys.L_local);
    method = erase(method, 'legacy_');

    if strcmp(method, 'sum')
        n_steps = 10000;
        if isfield(err, 'ee_legacy_steps')
            n_steps = max(20, round(err.ee_legacy_steps));
        end
        tau = linspace(0, 1, n_steps).';
        t_matrix = tau * t_flight;
        v_z = phys.v1 + a_z .* t_matrix;
        radius = R0 + 2 .* phys.R .* sin(phys.omega .* t_matrix ./ 2);
        field = I ./ (2 * pi * eps0 .* max(radius, realmin) .* ...
            max(abs(v_z), realmin));
        delta_ee = (e / m) .* sum(field .* (t_flight - t_matrix), 1) .* ...
            (t_flight ./ n_steps);
        n_points = n_steps;
    else
        n_points = 50;
        [tau, weights] = get_gauss_nodes(n_points);
        t_matrix = tau * t_flight;
        v_z = phys.v1 + a_z .* t_matrix;
        radius = R0 + 2 .* phys.R .* sin(phys.omega .* t_matrix ./ 2);
        field = I ./ (2 * pi * eps0 .* max(radius, realmin) .* ...
            max(abs(v_z), realmin));
        integral = weights * ((1 - tau) .* field);
        delta_ee = (e / m) .* t_flight.^2 .* integral;
    end

    out.delta_R_final = delta_ee;
    out.n_steps = n_points;
end

function value = finite_mean(values, fallback)
    [total, count] = block_sum(values, false, false);
    if count == 0
        value = fallback;
    else
        value = total / count;
    end
end

function value = finite_abs_mean(values, fallback)
    [total, count] = block_sum(values, true, false);
    if count == 0
        value = fallback;
    else
        value = total / count;
    end
end

function value = finite_max_abs(values, fallback)
    value = -Inf;
    block_size = 1e6;
    for first = 1:block_size:numel(values)
        block = double(values(first:min(first + block_size - 1, numel(values))));
        block = abs(block(isfinite(block)));
        if ~isempty(block), value = max(value, max(block)); end
    end
    if ~isfinite(value)
        value = fallback;
    end
end

function value = finite_rms(values, fallback)
    [sum_squares, count] = block_sum(values, false, true);
    if count == 0
        value = fallback;
    else
        value = sqrt(sum_squares / count);
    end
end

function value = positive_rms(values, fallback)
    [sum_squares, count] = block_sum(values, false, true, true);
    if count == 0
        value = fallback;
    else
        value = sqrt(sum_squares / count);
    end
end

function [total, count] = block_sum(values, use_absolute, use_square, positive_only)
    if nargin < 4, positive_only = false; end
    total = 0;
    count = 0;
    block_size = 1e6;
    for first = 1:block_size:numel(values)
        block = double(values(first:min(first + block_size - 1, numel(values))));
        valid = isfinite(block);
        if positive_only, valid = valid & block > 0; end
        block = block(valid);
        if isempty(block), continue, end
        if use_absolute, block = abs(block); end
        if use_square, block = block.^2; end
        total = total + sum(block);
        count = count + numel(block);
    end
end

function [x_mapped, w_mapped] = get_gauss_nodes(N_pts)
    i = 1:N_pts-1;
    a = i ./ sqrt(4 * i.^2 - 1);
    companion = diag(a, 1) + diag(a, -1);
    [vectors, eigenvalues] = eig(companion);
    [x, order] = sort(diag(eigenvalues));
    vectors = vectors(:, order).';
    weights = 2 * vectors(:, 1).^2;
    x_mapped = (x + 1) / 2;
    w_mapped = (weights / 2).';
end
