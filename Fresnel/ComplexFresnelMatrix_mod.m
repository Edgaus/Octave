function R = ComplexFresnelMatrix_mod(wave, n, thickness, theta0, polar)
    % Initialize the output array R to be the same size as the input 'wave' array
    num_waves = length(wave);
    R = zeros(size(wave));

    % Loop over every wavelength provided
    for k = 1:num_waves

        % 1. Extract the current wavelength
        w = wave(k);

        % 2. Extract the refractive index profile for THIS wavelength
        % If 'n' is a Cell Array (array of arrays): n_current = n{k}
        % If 'n' is a 2D Matrix: n_current = n(:, k)
        if iscell(n)
            n_current = n{k};
        elseif size(n, 2) == num_waves
            n_current = n(:, k);
        else
            n_current = n; % Fallback if n is constant across all wavelengths
        end

        % 3. Reset the Transfer Matrices for this specific wavelength
        M_P = eye(2);
        M_S = eye(2);

        % 4. Iterate through the layers
        for m = 2:length(n_current)-1

            theta_j = asin( (n_current(1)/n_current(m))*sin(theta0) );

            % Use the specific wavelength 'w' and specific index 'n_current'
            beta = (2*pi/w) * thickness(m-1) * n_current(m) * cos(theta_j);

            Mp = [ cos(beta)                              1i*cos(theta_j)*sin(beta)/n_current(m) ; ...
                   1i*n_current(m)*sin(beta)/cos(theta_j) cos(beta)];
            M_P = M_P * Mp;

            Ms = [ cos(beta)                                1i*sin(beta)/(n_current(m)*cos(theta_j)) ; ...
                   1i*n_current(m)*sin(beta)*cos(theta_j)   cos(beta)];
            M_S = M_S * Ms;

        end

        % 5. Calculate substrate angles and amplitudes
        theta_j_1 =  asin( (n_current(1)/n_current(end))*sin(theta0) );

        A_P = (1/(2*cos(theta0))) * [1 cos(theta0); 1 -cos(theta0) ] * M_P * [cos(theta_j_1) 0; n_current(end) 0];

        A_S = (1/(2*cos(theta0))) * [ cos(theta0) 1; cos(theta0) -1] * M_S * [1 0; cos(theta_j_1)*n_current(end) 0];

        r_P = A_P(2,1)/A_P(1,1);
        r_S = A_S(2,1)/A_S(1,1);

        % 6. Store the result in the output array R at index k
        if polar == 's'
            R(k) = abs(r_S)^2;
        elseif polar == 'p'
            R(k) = abs(r_P)^2;
        elseif polar == 'm'
            R(k) = 0.5 * (abs(r_P)^2 + abs(r_S)^2);
        end

    end % End of wavelength loop

end
