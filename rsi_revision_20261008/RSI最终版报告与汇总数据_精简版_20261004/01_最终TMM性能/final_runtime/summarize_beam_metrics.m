function metrics = summarize_beam_metrics(X, Y, X_ref, Y_ref)
% SUMMARIZE_BEAM_METRICS Basic centroid and spot-size metrics.

    metrics.cx = mean(X);
    metrics.cy = mean(Y);

    r = sqrt((X - metrics.cx).^2 + (Y - metrics.cy).^2);
    metrics.rms_radius = sqrt(mean(r.^2));

    r_sorted = sort(r);
    n_total = length(r_sorted);
    metrics.r_50 = r_sorted(max(1, round(0.50 * n_total)));
    metrics.r_90 = r_sorted(max(1, round(0.90 * n_total)));
    metrics.r_99 = r_sorted(max(1, round(0.99 * n_total)));
    metrics.FW50 = 2 * metrics.r_50;
    metrics.FW90 = 2 * metrics.r_90;
    metrics.FW99 = 2 * metrics.r_99;

    if nargin >= 4 && ~isempty(X_ref) && ~isempty(Y_ref)
        metrics.rms_blur_X = std(X - X_ref);
        metrics.rms_blur_Y = std(Y - Y_ref);
    else
        metrics.rms_blur_X = NaN;
        metrics.rms_blur_Y = NaN;
    end
end
