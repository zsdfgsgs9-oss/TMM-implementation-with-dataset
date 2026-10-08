function [x_final, y_final] = nonrel_intrinsic_coordinates(x, y, U0, alpha, phi, B, U, L)
%NONREL_INTRINSIC_COORDINATES Direct form of the intrinsic TMM translation.

    e = 1.602e-19;
    m = 9.109e-31;
    v_parallel = sqrt(2 * e .* U0 / m) .* cos(alpha);
    v_perp = sqrt(2 * e .* U0 / m) .* sin(alpha);
    a_z = e * U / (m * L);
    t = (-v_parallel + sqrt(v_parallel.^2 + 2 * a_z * L)) ./ a_z;
    omega = e * B / m;
    beta = omega .* t;
    r_larmor = m .* v_perp / (e * B);

    x_final = x + r_larmor .* (sin(phi + beta) - sin(phi));
    y_final = y + r_larmor .* (cos(phi) - cos(phi + beta));
end
