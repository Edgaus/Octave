
    #energy_lambda = 1240./lambda
    energy_lambda = linspace(-20,20,1000)

    #Gamma_x = Gamma_eV
    Gamma_x= 1

    #E0x = E_g_eV
    E0x = 0

    x = (energy_lambda -E0x)./Gamma_x

    A = 1

    f = @(y) ( (y.^2 + 1).^(-1/2) ) .* sqrt( sqrt(y.^2 + 1) + y );


    dedw = A*1i*(1/2)*(   Gamma_x^(-1/2)   )*(   f(x)  - 1i*f(-x)  )

   dedg = A*1i*(1/2)*(   Gamma_x^(-1/2)   )*(   f(-x)  + 1i*f(x)  )

    real_dedw = real( dedw )
    imag_dedw = imag( dedw )

    plot( energy_lambda, real_dedw , energy_lambda ,  imag_dedw  )



  %% Guardar los resultados en un archivo de texto

  % 1. Convertir los arreglos en vectores columna (verticales) usando (:)
  Wavelength_col = energy_lambda(:);
  Reflectance_col = imag_dedw(:);



  % 2. Unir ambas columnas en una sola matriz
  datos_finales = [Wavelength_col, Reflectance_col];
  fileID = fopen('Simulation\imag_dedw.txt', 'w');

% 4. Escribir los encabezados separados por una tabulación (\t)
  fprintf(fileID, 'Energy\tDialectric_function\n');
  fprintf(fileID, '%e\t%f\n', datos_finales');
% 6. Cerrar el archivo (¡Muy importante para que los datos se guarden en el disco!)
  fclose(fileID);
  disp('Los datos se han guardado exitosamente');







