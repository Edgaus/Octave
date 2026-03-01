function [total_error, R_opt, n_out, k_out] = Fit_adachi_best(x, lambda_nm, Energy_eV, y_exp)

    % 1. Unpack the optimizer's current guesses
    thick   = x(1);
    A0_g    = x(2);
    B1_g    = x(3);
    eps_g   = x(4);
    eps_g   = x(5);
    eps_g   = x(6);
    eps_g   = x(7);








    % Fixed physical parameters for AlN (you can add these to 'x' later if you want!)
    E0_fixed = 6.20; Gamma0_fixed = 0.15;
    E1_fixed = 8.70; Gamma1_fixed = 0.60;
    A_ex_fixed = 0.5; Eb_fixed = 0.060; Gamma_ex_fixed = 0.08;

    % 2. Generate n and k using your perfectly written Adachi function
    [n_aln, k_aln] = adachi(Energy_eV, E0_fixed, Gamma0_fixed, A0_g, E1_fixed, Gamma1_fixed, B1_g, A_ex_fixed, Eb_fixed, Gamma_ex_fixed, eps_g);

    % Combine into complex array
    n_complex_aln = n_aln - 1i * k_aln;

    % 3. Run the Fresnel Transfer Matrix
    % Silicon substrate is fixed. We interpolate it here for safety.
    si_n2 = importdata('Data\Si_n2.txt'); si_k2 = importdata('Data\Si_k2.txt');
    n_si = interp1(si_n2(:,1), si_n2(:,2), lambda_nm, 'pchip') - 1i*interp1(si_k2(:,1), si_k2(:,2), lambda_nm, 'pchip');

    R = zeros(1, length(lambda_nm));
    for m = 1:length(lambda_nm)
        n_structure = [1.0, n_complex_aln(m), n_si(m)];
        R(m) = ComplexFresnelMatrix(lambda_nm(m), n_structure, thick, deg2rad(45), 'm');
    end








% 4. Calculate the Custom Weighted Error
    % We still keep the massive penalty if it guesses outside our physical limits!
    if thick < 100 || thick > 300
        total_error = 1e6;
    elseif E0_g < 6.0 || E0_g > 6.4
        total_error = 1e6;
    elseif A0_g < 0 || B1_g < 0 || eps_g < 0
        total_error = 1e6;
    else
        % --- A. Build the Weight Array ---
        % Create an array of 1s the same size as lambda_nm
        W = ones(1, length(lambda_nm));

        % Find the indices where the wavelength is in the UV bandgap region (<= 225 nm)
        % and multiply their weight by 100
        uv_region = find(lambda_nm <= 225);
        W(uv_region) = 100;

        % --- B. Apply the Relative Absolute Error Formula ---
        % Force y_exp to be a row vector so the matrix math aligns perfectly with R
        y_exp_row = y_exp(:)';

        % Your exact formula: abs( R_sim / R_exp - 1 )
        relative_error = abs( (R ./ y_exp_row) - 1 );

        % Multiply by the weights and sum it all up!
        total_error = sum( W .* relative_error );
    end

end
