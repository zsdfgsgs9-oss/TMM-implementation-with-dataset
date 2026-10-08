function M = calc_M_EB0(p, phys, sys, err, N)
%CALC_M_EB0 Audited nonrelativistic EB0 map for the frozen Case 03/04 replay.
% Keep the archived sign and one-sided displacement convention. Only the
% initial-velocity term of the B0 focal shift changes from B^2 to B.

    if phys.is_relativistic
        error('This audited override is restricted to nonrelativistic EB0.');
    end

    e = 1.602e-19;
    m = 9.109e-31;
    Bz = phys.B_vec(3, :);
    E_nominal = phys.U_local ./ phys.L_local;

    delta_L_E = (2 * pi^2 * sys.N_turn^2 * m) ./ (e * Bz.^2) .* err.delta_E;
    term1_B = (4 * pi^2 * sys.N_turn^2 * m .* E_nominal) ./ (e * Bz.^2);
    term2_B = (2 * pi * sys.N_turn * m .* phys.v1) ./ (e * Bz);
    delta_L_B = -(err.delta_B ./ Bz) .* (term1_B + term2_B);
    delta_EB0 = -(delta_L_E + delta_L_B) .* phys.tan_gamma;

    M = repmat(eye(3), 1, 1, N);
    M(1, 3, :) = delta_EB0 .* cos(p.phi);
    M(2, 3, :) = delta_EB0 .* sin(p.phi);
end
