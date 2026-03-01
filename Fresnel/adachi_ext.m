function [n, k] = adachi_ext(E, E0, Gamma0, A0, E1, Gamma1, B1, eps_inf)

    % 1. Fundamental Gap (E0 Term)
    shi_0 = (E + 1i*Gamma0) ./ E0;
    eps_0 = A0 .* (E0^(-1.5)) .* (shi_0.^(-2)) .* ( 2 - sqrt(1+shi_0) - sqrt(1-shi_0) );

    % 2. Logarithmic Background (E1 Term)
    shi_1 = (E + 1i*Gamma1) ./ E1;
    eps_1 = -B1 .* (shi_1.^(-2)) .* log( 1 - shi_1.^2 );

    % Total Dielectric Function (Exciton removed)
    eps = eps_inf + eps_0 + eps_1;

    e_mag = abs(eps);
    e_real = real(eps);

    n = sqrt( (e_mag + e_real) ./ 2 );
    k = sqrt( (e_mag - e_real) ./ 2 );
end
