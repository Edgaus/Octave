function R = ComplexFresnelMatrix(wave, n, thickness, theta0, polar)

    M_P = eye(2);
    M_S = eye(2);
    for m = 2:length(n)-1

        theta_j = asin( (n(1)/n(m))*sin(theta0) );

        beta = (2*pi/wave)* thickness(m-1)*n(m)*cos(theta_j);


        Mp = [ cos(beta) 1i*cos(theta_j)*sin(beta)/n(m) ; 1i*n(m)*sin(beta)/cos(theta_j) cos(beta)];
        M_P = M_P*Mp;

        Ms = [ cos(beta) 1i*sin(beta)/(n(m)*cos(theta_j)) ; 1i*n(m)*sin(beta)*cos(theta_j) cos(beta)];
        M_S = M_S*Ms;

    end

    theta_j_1 =  asin( (n(1)/n(end))*sin(theta0) );


    A_P = (1/(2*cos(theta0))) * [1 cos(theta0); 1 -cos(theta0) ] * M_P * [cos(theta_j_1) 0; n(end) 0];

    A_S = (1/(2*cos(theta0))) * [ cos(theta0) 1; cos(theta0) -1] * M_S * [1 0; cos(theta_j_1)*n(end) 0];

    r_P = A_P(2,1)/A_P(1,1);
    r_S = A_S(2,1)/A_S(1,1);

    if polar=='s'
        R = abs(r_S)^2;

    elseif polar=='p'
        R = abs(r_P)^2;
    elseif polar =='m'
        R = 1/2  *(abs(r_P)^2 + abs(r_S)^2);
    end



end
