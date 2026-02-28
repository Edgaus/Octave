function [V,M,grid,x_grid] = grid_var_octave( Widths, Masses, potential )
    witdh = (Widths(1,1)+Widths(1,2))/2;
    upbound = 32;
    lowbound = 2;
    number_elments = 3;
    n =2 ;

    x=zeros(1,number_elments+1);
    x(1,1) = upbound;
    factor = zeros(1, number_elments+1);
    factor(1,1) = 1;

    for i = 2: number_elments
        factor(1,i) = 1 / ((1/n)^(i-1));
        x(1,i) = upbound/factor(1,i);
    end

    x(1,number_elments+1) = lowbound;
    factor(1, number_elments+1) = upbound/lowbound;



    fact = floor( witdh / ((number_elments+1)*upbound));
    factor = fact*factor;

    sobrante = witdh-(number_elments+1)*upbound*fact;

    m = 2 ;% fracción de elementos

    factor(end-1) = factor(end-1) + floor(sobrante/m/x(end-1)) ;
    factor(end) = factor(end) + floor(((m-1)/m)*sobrante/x(end));

    grid = [];
    for i=1:number_elments+1
        grid = [grid, repmat(x(1,i),1,factor(1,i))];
    end


    x_grid = zeros(1,length(grid));

    flip_grid = flip(grid);
    x_grid(1,1) = flip_grid(1,1);

    for i=2:length(grid)
        x_grid(1,i) = flip_grid(1,i)+x_grid(1,i-1);
    end

    x_grid = [ -1*flip(x_grid),0,x_grid ];
    grid = [grid, flip(grid)];

    funpar = @(z)   heaviside(z+Widths(1,1)/2)-heaviside(z-Widths(1,1)/2);
    V =[];
    M = [];
    for i =1:length(x_grid)
        V(i,1) = potential*(1-funpar(x_grid(1,i)));
        M(i,1) = Masses(1,2) +(Masses(1,1)-Masses(1,2))*( funpar(x_grid(1,i))  );
    end




end



