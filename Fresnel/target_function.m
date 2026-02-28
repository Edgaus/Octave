function [total_error, R_opt] = target_function(thick, n_aln, n_si, lambda_nm_its, true_x_max, true_x_min)

    R = zeros(1, length(lambda_nm_its));

    for m = 1:length(lambda_nm_its)
        % Structure: [Air, AlN, Silicon]
        n_structure = [1.0, n_aln(m), n_si(m)];

        % REMOVED the *100 so the output is a fraction (0 to 1) matching m668.txt
        R(m) = ComplexFresnelMatrix(lambda_nm_its(m), n_structure, thick, deg2rad(45), 'm');
    end

    R_opt = R; % Assign to output variable!
    y_sim = R;

    y_offset = min(y_sim);
    y_shifted = y_sim - y_offset + 1;

    % You may need to tweak min_h if simulated peaks are very shallow!
    min_h = 0.005;
    min_dist = 10;

    [~, sim_max_idx] = findpeaks(y_shifted, "MinPeakHeight", min_h, "MinPeakDistance", min_dist);
    [~, sim_min_idx] = findpeaks(-y_shifted + max(y_shifted), "MinPeakHeight", min_h, "MinPeakDistance", min_dist);

    sim_x_max = lambda_nm_its(sim_max_idx);
    sim_x_min = lambda_nm_its(sim_min_idx);

    total_error = 0;

    if isempty(sim_x_max) || isempty(sim_x_min)
        total_error = 1e6;
        return;
    end

    for i = 1:length(true_x_max)
        [~, closest_idx] = min(abs(sim_x_max - true_x_max(i)));
        matched_sim_x = sim_x_max(closest_idx);
        total_error = total_error + abs(matched_sim_x / true_x_max(i) - 1);
    end

    for i = 1:length(true_x_min)
        [~, closest_idx] = min(abs(sim_x_min - true_x_min(i)));
        matched_sim_x = sim_x_min(closest_idx);
        total_error = total_error + abs(matched_sim_x / true_x_min(i) - 1);
    end

    % I temporarily reduced this penalty so the optimizer doesn't get instantly
    % rejected if a minor simulated peak is missing.
    if length(sim_x_max) ~= length(true_x_max) || length(sim_x_min) ~= length(true_x_min)
        total_error = total_error + 0.5;
    end

end
