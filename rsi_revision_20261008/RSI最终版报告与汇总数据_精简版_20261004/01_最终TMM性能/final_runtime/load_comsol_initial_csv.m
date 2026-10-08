function [V_init, p] = load_comsol_initial_csv(filename)
%LOAD_COMSOL_INITIAL_CSV Load COMSOL release coordinates and velocities.

    data = readmatrix(filename, 'FileType', 'text');
    data = data(all(isfinite(data(:, 1:7)), 2), :);
    if size(data, 2) < 7
        error('Expected at least seven numeric columns in %s.', filename);
    end

    e = 1.602e-19;
    x = data(:, 2).' * 1e-6;
    y = data(:, 3).' * 1e-6;
    vx = data(:, 5).';
    vy = data(:, 6).';
    vz = data(:, 7).';

    V_init = zeros(3, 1, size(data, 1));
    V_init(1, 1, :) = x;
    V_init(2, 1, :) = y;
    V_init(3, 1, :) = 1;

    p.id = data(:, 1).';
    p.U0 = data(:, 4).' / e;
    p.alpha = atan2(hypot(vx, vy), abs(vz));
    p.phi = atan2(vy, vx);
    p.vx = vx;
    p.vy = vy;
    p.vz = vz;
end
