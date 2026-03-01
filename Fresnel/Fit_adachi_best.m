function [total_error, R_opt, n_out, k_out] = Fit_adachi_best(x, Ener_region, thick, y_exp_a, lambda_region)

    % 1. Unpack guesses
    E0_g        = x(1);
    Gamma0_g    = x(2);
    A0_g        = x(3);
    E1_g        = x(4);
    Gamma1_g    = x(5);
    B1_g        = x(6);
    eps_inf_g   = x(7);

    % 2. Generate n and k (Calling adachi_full correctly)
    [n_aln, k_aln] = adachi(Ener_region, E0_g, Gamma0_g, A0_g, E1_g, Gamma1_g, B1_g, eps_inf_g);
    n_complex_aln = n_aln - 1i * k_aln;

    % 3. Run Fresnel Transfer Matrix
    si_n2 = importdata('Data\Si_n2.txt');
    si_k2 = importdata('Data\Si_k2.txt');
    n_si = interp1(si_n2(:,1), si_n2(:,2), lambda_region, 'pchip') - 1i*interp1(si_k2(:,1), si_k2(:,2), lambda_region, 'pchip');

    R = zeros(1, length(lambda_region));
    for m = 1:length(lambda_region)
        n_structure = [1.0, n_complex_aln(m), n_si(m)];
        R(m) = ComplexFresnelMatrix(lambda_region(m), n_structure, thick, deg2rad(45), 'm');
    end

    % 4. Custom Weighted Error with RELAXED boundaries
    % Give the optimizer room to search!
    if E0_g < 5.5 || E0_g > 6.5
        total_error = 1e6;
    elseif A0_g < 0 || B1_g < 0 || eps_inf_g < 0 || Gamma0_g < 0 || Gamma1_g < 0
        total_error = 1e6;
    else
        W = ones(1, length(lambda_region));
        uv_region = find(lambda_region <= 225);
        W(uv_region) = 100;

        y_exp_row = y_exp_a(:)';
        relative_error = abs( (R ./ y_exp_row) - 1 );
        total_error = sum( W .* relative_error );
    end

    R_opt = R;
    n_out = n_aln;
    k_out = k_aln;

end
