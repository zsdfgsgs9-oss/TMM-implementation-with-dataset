function [beta, diagnostics] = calc_nonrel_magnetic_phase( ...
    vz0, a_z, flight_time, L, B_mid, X, Y, sys)
%CALC_NONREL_MAGNETIC_PHASE Nonrelativistic cyclotron phase along z(t).

    e = 1.602e-19;
    m = 9.109e-31;
    charge_to_mass = e / m;
    beta_midpoint = charge_to_mass .* B_mid .* flight_time;
    beta = beta_midpoint;
    method = "midpoint";

    if isfield(sys, 'nonrel_phase') && isfield(sys.nonrel_phase, 'method')
        method = lower(string(sys.nonrel_phase.method));
    end

    switch method
        case "midpoint"
        case "linear"
            if ~isfield(sys.nonrel_phase, 'G_Tpm')
                error('sys.nonrel_phase.G_Tpm is required for linear phase integration.');
            end
            G = sys.nonrel_phase.G_Tpm;
            path_moment = 0.5 .* vz0 .* flight_time.^2 + ...
                (1 / 6) .* a_z .* flight_time.^3 - ...
                0.5 .* L .* flight_time;
            beta = charge_to_mass .* (B_mid .* flight_time + G .* path_moment);
        case "gauss"
            if ~isfield(sys, 'func_B_vec')
                error('Gauss phase integration requires sys.func_B_vec.');
            end
            [nodes, weights] = gauss_legendre_12();
            beta = zeros(size(flight_time));
            for i = 1:numel(nodes)
                tau = 0.5 .* flight_time .* (nodes(i) + 1);
                z = vz0 .* tau + 0.5 .* a_z .* tau.^2;
                B_local = sys.func_B_vec(X, Y, z);
                beta = beta + 0.5 .* flight_time .* weights(i) .* B_local(3, :);
            end
            beta = charge_to_mass .* beta;
        otherwise
            error('Unknown nonrelativistic phase method "%s".', method);
    end

    delta = beta - beta_midpoint;
    diagnostics.method = char(method);
    diagnostics.mean_delta_rad = mean(delta);
    diagnostics.rms_delta_rad = sqrt(mean(delta.^2));
    diagnostics.max_abs_delta_rad = max(abs(delta));
end

function [nodes, weights] = gauss_legendre_12()
    nodes = [-0.981560634246719, -0.904117256370475, ...
        -0.769902674194305, -0.587317954286617, -0.367831498998180, ...
        -0.125233408511469, 0.125233408511469, 0.367831498998180, ...
        0.587317954286617, 0.769902674194305, 0.904117256370475, ...
        0.981560634246719];
    weights = [0.047175336386512, 0.106939325995318, ...
        0.160078328543346, 0.203167426723066, 0.233492536538355, ...
        0.249147045813403, 0.249147045813403, 0.233492536538355, ...
        0.203167426723066, 0.160078328543346, 0.106939325995318, ...
        0.047175336386512];
end
