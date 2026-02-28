function [Eigval, Eigfun] = Finite_Diferences(P,Mass, grid)


    const=3.8151;

    n = length(grid);
    L = zeros(1,n);
    L(1,1) = sqrt(grid(1,1));
    L(1,n) = sqrt(grid(1,n));
    for i =2:n-1
        L(1,i) = sqrt( (grid(1,i)+grid(1,i-1) )/2 );
    end



    for i=1:n-1
        offd1(i)= -const/(     (  (Mass(i) +Mass(i+1))*0.5  )    *grid(i)*L(i)^2       ) ; %off-diagonal + matrix element
    end

    for i=2:n
       offd2(i-1)= -const/(  (  (Mass(i) +Mass(i-1))*0.5  ) *grid(i-1)*L(i)^2 ) ; %off-diagonal + matrix element
    end

    for i = 2:n-1
        cent(1,i) = -(  offd2(i-1) +offd1(i) ) + P(1,i);
    end

    cent(1,1) = cent(1,2);
    cent(1,n) = cent(1,n-1);



    A = diag(cent) + diag(offd1,1) + diag(offd2,-1);

    [ Eigfun,Eigval] = eig(A);

    Eigval = diag(Eigval)*1000;

end





