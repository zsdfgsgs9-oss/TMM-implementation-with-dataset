function [V_init, p] = load_comsol_initial_csv(filename)
%LOAD_COMSOL_INITIAL_CSV Load COMSOL release coordinates and velocities.

    data = readmatrix(filename, 'FileType', 'text');
    if size(data,2)==9
        % Indexed export: row, pidx, x_um, y_um, z_um, K_J, vx, vy, vz.
        assert(all(abs(data(isfinite(data(:,1)),1)-data(isfinite(data(:,1)),2))<1e-8));
        data=data(:,[2 3 4 6 7 8 9]);
    elseif size(data,2)~=7
        error('Expected the documented seven- or nine-column launch export.');
    end
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
    assert(all(abs(p.id-round(p.id))<1e-8) && numel(unique(p.id))==numel(p.id));
    p.U0 = data(:, 4).' / e;
    p.alpha = atan2(hypot(vx, vy), abs(vz));
    p.phi = atan2(vy, vx);
    p.vx = vx;
    p.vy = vy;
    p.vz = vz;
end
