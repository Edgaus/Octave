function total_error = target_function_nk(x, Energy_eV, n_target, k_target)

    % 1. Unpack the 7 guesses
    E0_g      = x(1); A0_g      = x(2);
    E1_g      = x(3);; B1_g      = x(4);
    eps_inf_g = x(5);

    % 2. THE PENALTY WALL
    % We MUST enforce eps_inf > 2.0 so the n curve doesn't collapse to zero!
    % We MUST enforce Gamma > 0 so the math doesn't crash into a singularity!
    if E0_g < 5 || E0_g > 6.3 || A0_g < 0  ||  B1_g < 0
        total_error = 1e18;
        return;
    end

    % 3. Pass 7 guesses into the math engine
    [n_sim, k_sim] = adachi_ext(Energy_eV, E0_g, A0_g, E1_g, B1_g, eps_inf_g);

    % 4. CALCULATE POINT-BY-POINT ERRORS
    error_n_array = (abs((n_sim ./ n_target(:)') - 1)).^2;
    error_k_array = (abs(k_sim - k_target(:)')).^2;

    % 5. BUILD THE DYNAMIC WEIGHT ARRAY
    lambda_current = 1240 ./ Energy_eV;
    W = ones(1, length(lambda_current));

    bandgap_region = find(lambda_current >= 190 & lambda_current <= 225);
    W(bandgap_region) = 1000;

    % 6. APPLY WEIGHTS AND SUM FOR FINAL SCORE
    weighted_error_n = sum(W .* error_n_array);
    weighted_error_k =  sum(W .* error_k_array);


    total_error = (  weighted_error_n + weighted_error_k);
end
