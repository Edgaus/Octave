function [total_error, R_opt] = target_function_peaks(thick, n_aln, n_si, lambda_nm, exp_x_max, exp_x_min)

    % --- A. Generate Simulated Curve ---
    R = zeros(1, length(lambda_nm));
    for m = 1:length(lambda_nm)
        n_structure = [1.0, n_aln(m), n_si(m)];
        R(m) = ComplexFresnelMatrix(lambda_nm(m), n_structure, thick, deg2rad(45), 'm');
    end
    R_opt = R;

    % --- B. Find Peaks in Simulated Curve dynamically ---
    y_offset = min(R);
    y_shifted = R - y_offset + 1;

    % Thresholds for simulated data (simulations are usually smooth, so we can use low thresholds)
    min_h = 0.005;
    min_dist = 40; % Wait ~8.8 nm before finding another peak to avoid fake noise peaks

    [~, sim_max_idx] = findpeaks(y_shifted, "MinPeakHeight", min_h, "MinPeakDistance", min_dist);
    [~, sim_min_idx] = findpeaks(-y_shifted + max(y_shifted), "MinPeakHeight", min_h, "MinPeakDistance", min_dist);

    sim_x_max = lambda_nm(sim_max_idx);
    sim_x_min = lambda_nm(sim_min_idx);

    % --- C. Apply User's Custom Formula Error ---
    total_error = 0;

    % Safety Catch: if simulation flatlines
    if isempty(sim_x_max) || isempty(sim_x_min)
        total_error = 1e6;
        return;
    end

    % Compare Maxima
    for i = 1:length(exp_x_max)
        [~, closest_idx] = min(abs(sim_x_max - exp_x_max(i)));
        matched_sim_x = sim_x_max(closest_idx);
        % sum( abs(sim_x / raw_x - 1) )
        total_error = total_error + ( abs(matched_sim_x / exp_x_max(i) - 1) ).^2;
    end

    % Compare Minima
    for i = 1:length(exp_x_min)
        [~, closest_idx] = min(abs(sim_x_min - exp_x_min(i)));
        matched_sim_x = sim_x_min(closest_idx);
        total_error = total_error +( abs(matched_sim_x / exp_x_min(i) - 1)).^2;
    end

    % --- D. Apply Mismatch Penalty ---
    % If the thickness guess creates 5 peaks but the experiment only has 3,
    % we MUST punish it so it doesn't choose that thickness.
    if length(sim_x_max) ~= length(exp_x_max) || length(sim_x_min) ~= length(exp_x_min)
        total_error = total_error + 10; % Massive penalty
    end

end
