clear
clearvars

A=1

FP = @(x) A*sqrt(sqrt(x.^2+1)+x)./sqrt(x.^2+1);
FM = @(x)  A*sqrt(sqrt(x.^2+1)-x)./sqrt(x.^2+1);
FPN = @(x) -A*sqrt(sqrt(x.^2+1)+x)./sqrt(x.^2+1);


Gamma = 1

x = linspace(-10, 10, 100)

xc = 0

W = (x-xc)/Gamma

dedw = 1i*( 1/2)*( Gamma^(-1/2) )*( FP(W)-1i*FM( W))
dedG = 1i*( 1/2)*( Gamma^(-1/2) )*( FM(W)+1i*FP( W))

plot( x, real( dedw), x, imag( dedw  ))

plot( x, real( dedG), x, imag( dedG  ))


Wavelength_col = x(:);
der_dw =  real ( dedw (:) );
dei_dw = imag( dedw (:) );

der_dG =  real ( dedG (:) );
dei_dG = imag( dedG (:) );




% 2. Unir ambas columnas en una sola matriz
datos_finales = [Wavelength_col, der_dw, dei_dw,  der_dG, dei_dG  ];

% 3. Abrir (o crear) el archivo de texto en modo escritura ('w')
% Puedes cambiar 'Resultados_Reflectancia.txt' por el nombre que prefieras
fileID = fopen('Simulation\F_W_behaviour.txt', 'w');

% 4. Escribir los encabezados separados por una tabulación (\t)
fprintf(fileID, 'Wavelength\tder_dw\tdei_dw\tder_dG\tdei_dG\n');

% 5. Escribir los datos numéricos
% Usamos %e para Wavelength (porque está en metros, ej. 2.00e-07)
% y %f para Reflectance (números decimales normales)
% OJO: Es necesario transponer la matriz (datos_finales') dentro de fprintf
% porque Octave lee los datos columna por columna al escribirlos.
fprintf(fileID, '%e\t%f\t%f\t%f\t%f\n', datos_finales');

% 6. Cerrar el archivo (¡Muy importante para que los datos se guarden en el disco!)
fclose(fileID);


