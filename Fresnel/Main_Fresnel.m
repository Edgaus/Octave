clear
clearvars
pkg load tablicious

%% Lectura de archivos iniciales


%% Film metrics



fim_spectra =  importdata('Experimental\m668.txt')



% 3. Extract your columns into two arrays
Wavelength_film = fim_spectra(:, 1)*1e-9 ;
film_spe = fim_spectra(:, 2)*100;






%% AlN



index_refrac_aln =  importdata('Data\AlN_n2.txt')



% 3. Extract your columns into two arrays
Wavelength_index_aln = index_refrac_aln(:, 1)*1e-9 ;
aln_index = index_refrac_aln(:, 2);


extintion_coeff =  importdata('Data\AlN_k2.txt')



% 3. Extract your columns into two arrays
Wavelength_extin_aln = extintion_coeff(:, 1)*1e-9;
aln_extin = extintion_coeff(:, 2);



%% Si



index_refrac_si =  importdata('Data\Si_n2.txt')



% 3. Extract your columns into two arrays
Wavelength_index_si = index_refrac_si(:, 1)*1e-9 ;
si_index = index_refrac_si(:, 2);


extintion_coeff_si =  importdata('Data\Si_k2.txt')



% 3. Extract your columns into two arrays
Wavelength_extin_si = extintion_coeff_si(:, 1)*1e-9 ;
si_extin = extintion_coeff_si(:, 2);







%%  Interpolate complex indices

lambda_nm = linspace(150,400,750)*1e-9 ;  % [nm]


n0 = 1.0;  % air


n_aln = interp1(Wavelength_index_aln,  aln_index,  lambda_nm,'pchip') ...
     - 1i*interp1(Wavelength_extin_aln,  aln_extin,  lambda_nm,'pchip');

n_si = interp1(Wavelength_index_si,  si_index,  lambda_nm,'pchip') ...
     - 1i*interp1(Wavelength_extin_si,  si_extin,  lambda_nm,'pchip');




%% Initial Conditions


R = zeros(1,length(lambda_nm));


d = [188.19660]*1e-9;     %m905


theta_inc =45;


%% Initial Code



for m=1:length(lambda_nm)
      n_structure = [n0,n_aln(m), n_si(m)];
      R(m) = ComplexFresnelMatrix(lambda_nm(m), n_structure, d, deg2rad(theta_inc), 'm' )*100;
end



figure;
plot(lambda_nm,R, Wavelength_film, film_spe );



%% Guardar los resultados en un archivo de texto

% 1. Convertir los arreglos en vectores columna (verticales) usando (:)
Wavelength_col = lambda_nm(:)*1e9;
Reflectance_col = R(:);

% 2. Unir ambas columnas en una sola matriz
datos_finales = [Wavelength_col, Reflectance_col];

% 3. Abrir (o crear) el archivo de texto en modo escritura ('w')
% Puedes cambiar 'Resultados_Reflectancia.txt' por el nombre que prefieras
fileID = fopen('Simulation\sim.txt', 'w');

% 4. Escribir los encabezados separados por una tabulación (\t)
fprintf(fileID, 'Wavelength\tReflectance\n');

% 5. Escribir los datos numéricos
% Usamos %e para Wavelength (porque está en metros, ej. 2.00e-07)
% y %f para Reflectance (números decimales normales)
% OJO: Es necesario transponer la matriz (datos_finales') dentro de fprintf
% porque Octave lee los datos columna por columna al escribirlos.
fprintf(fileID, '%e\t%f\n', datos_finales');

% 6. Cerrar el archivo (¡Muy importante para que los datos se guarden en el disco!)
fclose(fileID);

disp('Los datos se han guardado exitosamente en Resultados_Reflectancia.txt');

