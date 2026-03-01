function [n, k] = adachi_ext(E, E0, Gamma0, A0, E1, Gamma1, B1, eps_inf, A_ex)
    % E0 Term
    shi_0 = (E + 1i*Gamma0) ./ E0;
    eps_0 = A0 .* (E0^(-1.5)) .* (shi_0.^(-2)) .* ( 2 - sqrt(1+shi_0) - sqrt(1-shi_0) );

    % Exciton Term (Now using the free parameter A_ex!)
    Eb = 0.060;
    Gamma_ex = 0.08;
    E_ex = E0 - Eb;
    eps_ex = A_ex ./ (E_ex - E - 1i*Gamma_ex);

    % E1 Term
    shi_1 = (E + 1i*Gamma1) ./ E1;
    eps_1 = -B1 .* (shi_1.^(-2)) .* log( 1 - shi_1.^2 );

    % Total
    eps = eps_inf + eps_0 + eps_ex + eps_1;

    e_mag = abs(eps);
    e_real = real(eps);

    n = sqrt( (e_mag + e_real) ./ 2 );
    k = sqrt( (e_mag - e_real) ./ 2 );
end
