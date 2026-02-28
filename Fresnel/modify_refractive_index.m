function n_index_mod = modify_refractive_index(Gamma_eV, E_g_eV, n_index, lambda)
  % MODIFY_REFRACTIVE_INDEX_FUNC Adds a custom function curve to n_index.

  n_index_mod = n_index;

  #energy_lambda = 1240./lambda
  energy_lambda = linspace(-10,10,100)

  #Gamma_x = Gamma_eV
  Gamma _x= 2

  #E0x = E_g_eV
  E0x = 0

  x = (energy_lambda -E0x)./Gamma_x

  A = 1

  f = @(y)      (     ( y.^2 +1 ).^(-1/2)    )*sqrt(   sqrt(   y.^2+1   )     +y     )

  der = A*(1/2)*(Gamma^(-1/2))*( f(-x) )


  plot( x, der)


end
