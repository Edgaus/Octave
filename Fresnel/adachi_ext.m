function [n, k] = adachi_ext(E, E0, A0, E1, B1, eps_inf)

    % 1. Fundamental Gap (E0 Term)
    shi_0 = (E ) ./ E0;
    eps_0 = A0 .* (E0^(-1.5)) .* (shi_0.^(-2)) .* ( 2 - sqrt(1+shi_0) - sqrt(1-shi_0) );

    % 2. Spin-Orbit Split-Off Gap
    Delta_so = 0.020; % Usually around 50 meV
    E_so = E0 + Delta_so;
    shi_so = (E ) ./ E_so;

    % FIXED: The Spin-Orbit transition strength is theoretically 0.5 * A0
    eps_so = (0.5 * A0) .* (E_so^(-1.5)) .* (shi_so.^(-2)) .* ( 2 - sqrt(1+shi_so) - sqrt(1-shi_so) );

    % 3. Logarithmic Background (E1 Term)
    shi_1 = (E ) ./ E1;
    eps_1 = -B1 .* (shi_1.^(-2)) .* log( 1 - shi_1.^2 );

    % Total Dielectric Function
    eps = eps_inf + eps_0 + eps_so + eps_1;

    e_mag = abs(eps);
    e_real = real(eps);

    n = sqrt( (e_mag + e_real) ./ 2 );
    k = sqrt( (e_mag - e_real) ./ 2 );
end
