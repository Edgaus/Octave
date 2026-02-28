function [x_max, x_min] = der_functi( name_sample  )
  pkg load signal
  sam = name_sample

  file_name = ['Experimental\', sam,'.txt'  ]
% Load the data (assuming two columns: Time and Value)
  data = load(file_name);

x_full = data(:, 1);
y_full = data(:, 2);

% --- 2. Extract Region of Interest (225 < X < 395) ---
region_idx = find(x_full > 225 & x_full < 395);
x = x_full(region_idx);
y_raw = y_full(region_idx);

% --- 3. Split the Data at X = 300 ---
split_idx = find(x >= 300, 1);
x_part1 = x(1:split_idx-1);
y_part1 = y_raw(1:split_idx-1);
x_part2 = x(split_idx:end);
y_part2 = y_raw(split_idx:end);

% --- 4. Smooth ONLY Part 2 ---
window_size = 5;
y_part2_smooth = movmean(y_part2, window_size);

% --- 5. Recombine the Data ---
y_combined = [y_part1; y_part2_smooth];

% --- 6. Safe Shift for Peak Detection ---
y_offset = min(y_combined);
y_shifted = y_combined - y_offset + 1;

% --- 7. Set Thresholds ---
min_h = 0.02;
min_dist = 10;

% --- 8. Find Maxima ---
[pks_shifted, max_idx] = findpeaks(y_shifted, "MinPeakHeight", min_h, "MinPeakDistance", min_dist);
y_max = y_combined(max_idx);
x_max = x(max_idx); % <--- HERE ARE YOUR SEPARATED MAX X VALUES

% --- 9. Find Minima ---
[vls_shifted, min_idx] = findpeaks(-y_shifted + max(y_shifted), "MinPeakHeight", min_h, "MinPeakDistance", min_dist);
y_min = y_combined(min_idx);
x_min = x(min_idx); % <--- HERE ARE YOUR SEPARATED MIN X VALUES

% ==========================================
% NEW: SEPARATING, PRINTING, AND SAVING
% ==========================================



% --- Plotting (Same as before) ---
figure(1);
plot(x, y_combined, 'k-', 'LineWidth', 1.5); hold on;
plot(x_max, y_max, 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'Maxima');
plot(x_min, y_min, 'bv', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Minima');
plot([300 300], ylim(), 'g--', 'LineWidth', 2, 'DisplayName', 'Smoothing Boundary');
title('Separated Extrema for 225 < X < 395');
legend('Location', 'northeastoutside');
grid on;
xlim([220, 400]);
end
