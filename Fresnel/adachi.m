function [n, k] = adachi_full(E, E0, Gamma0, A0, E1, Gamma1, B1, A_ex, Eb, Gamma_ex, eps_inf)

    % ----------------------------------------------------
    % 1. E0 Term (3D Critical Point: Free Electron-Hole)
    % ----------------------------------------------------
    shi_0 = (E + 1i*Gamma0) ./ E0;

    % Using sqrt() instead of .^(1/2) is slightly faster in Octave
    eps_0 = A0 .* (E0^(-1.5)) .* (shi_0.^(-2)) .* ( 2 - sqrt(1+shi_0) - sqrt(1-shi_0) );

    % ----------------------------------------------------
    % 2. Exciton Term (Discrete Bound State)
    % ----------------------------------------------------
    E_ex = E0 - Eb; % The exciton resonance sits just below the bandgap
    eps_ex = A_ex ./ (E_ex - E - 1i*Gamma_ex);

    % ----------------------------------------------------
    % 3. E1 Term (2D Logarithmic Critical Point)
    % ----------------------------------------------------
    shi_1 = (E + 1i*Gamma1) ./ E1;
    eps_1 = -B1 .* (shi_1.^(-2)) .* log( 1 - shi_1.^2 );

    % ----------------------------------------------------
    % 4. Total Complex Dielectric Function
    % ----------------------------------------------------
    eps = eps_inf + eps_0 + eps_1

    % ----------------------------------------------------
    % 5. Extract n and k robustly
    % ----------------------------------------------------
    e_mag = abs(eps);
    e_real = real(eps);

    n = sqrt( (e_mag + e_real) ./ 2 );
    k = sqrt( (e_mag - e_real) ./ 2 );

end
