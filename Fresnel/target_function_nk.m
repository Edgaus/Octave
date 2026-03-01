function total_error = target_function_nk(x, Energy_eV, n_target, k_target)

    % 1. Unpack all 8 guesses
    E0_g      = x(1); Gamma0_g  = x(2); A0_g      = x(3);
    E1_g      = x(4); Gamma1_g  = x(5); B1_g      = x(6);
    eps_inf_g = x(7); A_ex_g    = x(8); % <--- Unpack A_ex

    % 2. MASSIVE penalty wall to strictly enforce physics
    if E0_g < 5.8 || E0_g > 6.2|| A0_g < 0 || B1_g < 0 || eps_inf_g < 0 || Gamma0_g < 0 || Gamma1_g < 0 || A_ex_g < 0
        total_error = 1e12;
        return;
    end

    % 3. Pass A_ex_g into the adachi function!
    [n_sim, k_sim] = adachi_ext(Energy_eV, E0_g, Gamma0_g, A0_g, E1_g, Gamma1_g, B1_g, eps_inf_g, A_ex_g);

    % ========================================================
    % 4. CALCULATE POINT-BY-POINT ERRORS (AS ARRAYS)
    % ========================================================
    error_n_array = abs((n_sim ./ n_target(:)') - 1);
    error_k_array = abs(k_sim - k_target(:)');

    % ========================================================
    % 5. BUILD THE DYNAMIC WEIGHT ARRAY
    % ========================================================
    lambda_current = 1240 ./ Energy_eV;
    W = ones(1, length(lambda_current));

    % Find the bandgap region
    bandgap_region = find(lambda_current >= 190 & lambda_current <= 225);

    % Drop the weight from 100 down to 10.
    % This is enough to guide the optimizer without blinding it to the flat tail!
    W(bandgap_region) = 100;

    % ========================================================
    % 6. APPLY WEIGHTS AND SUM FOR FINAL SCORE
    % ========================================================
    weighted_error_n = sum(W .* error_n_array);
    weighted_error_k = sum(W .* error_k_array);

    % Keep the global k multiplier active
    global_k_multiplier = 100;

    total_error = weighted_error_n + (global_k_multiplier * weighted_error_k);

end
