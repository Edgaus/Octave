function [x_max, x_min] = der_functi(sample)
    % 1. Load the data
    file_name = ['Experimental\', sample, '.txt'];
    data = load(file_name);
    x_full = data(:, 1);
    y_raw_full = data(:, 2);

    % --- NEW: Crop the data to the region of interest (225 - 390 nm) ---
    % This prevents the massive noise spikes below 225nm from being counted as peaks.
    region_idx = find(x_full >= 225 & x_full <= 390);
    x = x_full(region_idx);
    y_raw = y_raw_full(region_idx);

    % 2. Split the data at X = 300
    split_idx = find(x >= 300, 1);
    y_part1 = y_raw(1:split_idx-1);
    y_part2 = y_raw(split_idx:end);

    % 3. Apply heavy smoothing ONLY to the noisy part (>300)
    window_size = 25;
    y_part2_smooth = movmean(y_part2, window_size);

    % 4. Recombine
    y_combined = [y_part1; y_part2_smooth];

    % 5. Safe Shift and Thresholds
    y_offset = min(y_combined);
    y_shifted = y_combined - y_offset + 1;

    % Slightly stricter thresholds to ensure exactly 2 peaks
    min_h = 0.01;  % Ignore tiny ripples
    min_dist = 60; % 60 points = 15 nm minimum separation

    % 6. Find Peaks and Valleys
    [~, max_idx] = findpeaks(y_shifted, "MinPeakHeight", min_h, "MinPeakDistance", min_dist);
    [~, min_idx] = findpeaks(-y_shifted + max(y_shifted), "MinPeakHeight", min_h, "MinPeakDistance", min_dist);

    x_max = x(max_idx);
    x_min = x(min_idx);

    % Optional: Uncomment these lines to visually verify it only found 2 of each!
    % figure(99); clf;
    % plot(x, y_combined, 'k-'); hold on;
    % plot(x_max, y_combined(max_idx), 'r^', 'MarkerSize', 10);
    % plot(x_min, y_combined(min_idx), 'bv', 'MarkerSize', 10);
end
