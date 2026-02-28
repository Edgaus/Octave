function [x_max, x_min] = der_functi(sample)
    % 1. Load the data
    file_name = ['Experimental\', sample, '.txt'];
    data = load(file_name);
    x = data(:, 1);
    y_raw = data(:, 2);

    % 2. Split the data at X = 300
    split_idx = find(x >= 300, 1);
    y_part1 = y_raw(1:split_idx-1);
    y_part2 = y_raw(split_idx:end);

    % 3. Apply heavy smoothing ONLY to the noisy part (>300)
    % A window size of 25 averages over ~6 nm of data, killing the noise
    window_size = 25;
    y_part2_smooth = movmean(y_part2, window_size);

    % 4. Recombine
    y_combined = [y_part1; y_part2_smooth];

    % 5. Safe Shift and Thresholds
    y_offset = min(y_combined);
    y_shifted = y_combined - y_offset + 1;

    min_h = 0.005;
    min_dist = 40; % 40 points = 10 nm minimum separation

    % 6. Find Peaks and Valleys
    [~, max_idx] = findpeaks(y_shifted, "MinPeakHeight", min_h, "MinPeakDistance", min_dist);
    [~, min_idx] = findpeaks(-y_shifted + max(y_shifted), "MinPeakHeight", min_h, "MinPeakDistance", min_dist);

    x_max = x(max_idx);
    x_min = x(min_idx);

    % 7. Diagnostic Plot (Pops up automatically so you can check it)
    figure(99); clf;
    plot(x, y_raw, 'color', [0.8 0.8 0.8], 'DisplayName', 'Raw Noise'); hold on;
    plot(x, y_combined, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Smoothed Signal');
    plot(x_max, y_combined(max_idx), 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'Found Maxima');
    plot(x_min, y_combined(min_idx), 'bv', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Found Minima');
    plot([300 300], ylim(), 'g--', 'DisplayName', 'Smoothing Boundary');
    title(['Extracted True Peaks for ', sample]);
    xlim([225, 390]);
    legend();
end
