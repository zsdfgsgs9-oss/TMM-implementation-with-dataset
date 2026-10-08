function run_mc_convergence()
%RUN_MC_CONVERGENCE Run the 6 particle counts by 10 seed TMM study.

    this_dir = fileparts(mfilename('fullpath'));
    root_dir = fileparts(this_dir);
    output_dir = fullfile(this_dir, 'results');
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    addpath(root_dir);
    addpath(this_dir);

    counts = [1e3, 1e4, 1e5, 1e6, 1e7, 5e7];
    seeds = 42:51;
    chunk_size = 1e6;
    rows = table();
    for i = 1:numel(counts)
        for j = 1:numel(seeds)
            n = counts(i);
            seed = seeds(j);
            fprintf('MC N=%g, seed=%d\n', n, seed);
            timer = tic;
            metrics = run_one_sample(n, seed, chunk_size);
            one = table(n, seed, metrics.cx * 1e9, metrics.cy * 1e9, ...
                metrics.rms_radius * 1e9, metrics.FW50 * 1e9, ...
                metrics.FW90 * 1e9, metrics.FW99 * 1e9, toc(timer), ...
                'VariableNames', {'N', 'seed', 'cx_nm', 'cy_nm', 'rms_nm', ...
                'FW50_nm', 'FW90_nm', 'FW99_nm', 'runtime_s'});
            rows = [rows; one]; %#ok<AGROW>
            writetable(rows, fullfile(output_dir, 'mc_seed_summary.csv'));
        end
    end

    summary = summarize_seeds(rows);
    writetable(summary, fullfile(output_dir, 'mc_convergence_summary.csv'));
    save(fullfile(output_dir, 'mc_workspace.mat'), 'rows', 'summary');
end

function metrics = run_one_sample(n, seed, chunk_size)
    if n <= 1e6
        rng(seed, 'twister');
        [x, y, U0, alpha, phi] = generate_particles(n);
        [xf, yf] = nonrel_intrinsic_coordinates(x, y, U0, alpha, phi, ...
            1.2, 12.67e3, 1e-3);
        metrics = summarize_beam_metrics(xf, yf);
        return;
    end

    rng(seed, 'twister');
    sx = 0; sy = 0; sr2 = 0; max_r0 = 0;
    done = 0;
    while done < n
        k = min(chunk_size, n - done);
        [x, y, U0, alpha, phi] = generate_particles(k);
        [xf, yf] = nonrel_intrinsic_coordinates(x, y, U0, alpha, phi, ...
            1.2, 12.67e3, 1e-3);
        sx = sx + sum(xf); sy = sy + sum(yf);
        sr2 = sr2 + sum(xf.^2 + yf.^2);
        max_r0 = max(max_r0, max(hypot(xf, yf)));
        done = done + k;
    end
    cx = sx / n; cy = sy / n;
    rms_radius = sqrt(max(sr2 / n - cx^2 - cy^2, 0));

    bin_width = 5e-12;
    r_max = max_r0 + hypot(cx, cy) + bin_width;
    edges = 0:bin_width:(ceil(r_max / bin_width) * bin_width + bin_width);
    counts = zeros(1, numel(edges) - 1);
    rng(seed, 'twister');
    done = 0;
    while done < n
        k = min(chunk_size, n - done);
        [x, y, U0, alpha, phi] = generate_particles(k);
        [xf, yf] = nonrel_intrinsic_coordinates(x, y, U0, alpha, phi, ...
            1.2, 12.67e3, 1e-3);
        counts = counts + histcounts(hypot(xf - cx, yf - cy), edges);
        done = done + k;
    end

    metrics.cx = cx;
    metrics.cy = cy;
    metrics.rms_radius = rms_radius;
    metrics.r_50 = histogram_quantile(counts, edges, 0.50);
    metrics.r_90 = histogram_quantile(counts, edges, 0.90);
    metrics.r_99 = histogram_quantile(counts, edges, 0.99);
    metrics.FW50 = 2 * metrics.r_50;
    metrics.FW90 = 2 * metrics.r_90;
    metrics.FW99 = 2 * metrics.r_99;
end

function [x, y, U0, alpha, phi] = generate_particles(n)
    radius = 18e-9 * sqrt(rand(1, n));
    theta = 2 * pi * rand(1, n);
    x = radius .* cos(theta);
    y = radius .* sin(theta);
    U0 = 0.5 * rand(1, n);
    phi = 2 * pi * rand(1, n);
    alpha = asin(sqrt(sin(deg2rad(10.7))^2 * rand(1, n)));
end

function q = histogram_quantile(counts, edges, probability)
    cumulative = cumsum(counts);
    target = probability * cumulative(end);
    idx = find(cumulative >= target, 1, 'first');
    before = 0;
    if idx > 1, before = cumulative(idx - 1); end
    fraction = (target - before) / max(counts(idx), 1);
    q = edges(idx) + fraction * (edges(idx + 1) - edges(idx));
end

function summary = summarize_seeds(rows)
    unique_n = unique(rows.N);
    summary = table();
    names = {'FW50_nm', 'FW90_nm', 'FW99_nm'};
    tcrit = 2.262157162740992;
    for i = 1:numel(unique_n)
        subset = rows(rows.N == unique_n(i), :);
        one.N = unique_n(i);
        for j = 1:numel(names)
            values = subset.(names{j});
            base = erase(names{j}, '_nm');
            one.([base '_mean_nm']) = mean(values);
            one.([base '_std_nm']) = std(values);
            half = tcrit * std(values) / sqrt(numel(values));
            one.([base '_ci95_low_nm']) = mean(values) - half;
            one.([base '_ci95_high_nm']) = mean(values) + half;
        end
        summary = [summary; struct2table(one)]; %#ok<AGROW>
    end
end
