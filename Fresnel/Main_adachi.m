1; % <--- DO NOT REMOVE. Required for Octave scripts with local functions.

pkg load signal;
clear all;

% ========================================================
% 1. LOAD AND ALIGN DATA
% ========================================================
sample = 'm669';
file_name = ['Experimental\', sample, '.txt'];
data = load(file_name);

x_raw = data(:, 1); % Wavelength in nm
y_raw = data(:, 2); % Reflectance (0 to 1)

% Safely extract the exact region of interest DIRECTLY from the data
region_idx = find(x_raw >= 200 & x_raw <= 400.75);
lambda_nm = x_raw(region_idx);
y_exp = y_raw(region_idx);

Energy_eV = 1240 ./ lambda_nm;

% ========================================================
% 2. INITIALIZE PARAMETERS
% ========================================================
if strcmp(sample, 'm668')
    thickness = 162.31601;
elseif strcmp(sample, 'm669')
    thickness = 211.43622;
elseif strcmp(sample, 'm672')
    thickness = 193.77740;
elseif strcmp(sample, 'm687')
    thickness = 204.78842;
elseif strcmp(sample, 'm691')
    thickness = 176.29674;
else
    thickness = 200; % Fallback
end

% Guesses: [E0, A0, E1, B1, eps_inf]
x_guess = [6, 30.0, 10, 2.0, 4.0,-2];

options = optimset('Display', 'iter', 'MaxIter', 2000, 'MaxFunEvals', 4000, 'TolX', 1e-4);

disp('Starting Master Optical Optimization...');

% ========================================================
% 3. RUN FMINSEARCH
% ========================================================
[best_params, final_error] = fminsearch(@(x) Fit_adachi_best(x, Energy_eV, thickness, y_exp, lambda_nm), x_guess, options);

% Unpack the winning parameters
best_E0       = best_params(1);
best_A0       = best_params(2);
best_E1       = best_params(3);
best_B1       = best_params(4);
best_eps_inf  = best_params(5);

fprintf('\n=== OPTIMIZATION COMPLETE ===\n');
fprintf('E0:         %.3f eV\n', best_E0 );
fprintf('A0:         %.2f\n', best_A0);
fprintf('E1:         %.2f\n', best_E1);
fprintf('B1:         %.2f\n', best_B1);
fprintf('eps_inf:    %.2f\n', best_eps_inf);

% ========================================================
% 4. EXTRACT AND PLOT
% ========================================================
[~, R_simulated, n_final, k_final] = Fit_adachi_best(best_params, Energy_eV, thickness, y_exp, lambda_nm);

% Save Reflectance
% ========================================================
Wavelength_col = lambda_nm(:);
Reflectance_sim_col = R_simulated(:);
Reflectance_exp_col = y_exp(:);

file_name1 = ['Simulation\', sample, 'Reflec_sim_.txt'];
datos_finales1 = [Wavelength_col, Reflectance_sim_col, Reflectance_exp_col];

% Ensure the 'Simulation' folder exists in your directory!
fileID = fopen(file_name1, 'w');

% CORRECCIÓN: Le agregamos el texto del encabezado
fprintf(fileID, 'Wavelength(nm)\tR_simulated\tR_experimental\n');

% Guardamos los números
fprintf(fileID, '%e\t%f\t%f\n', datos_finales1');
fclose(fileID);

disp(['Los datos se han guardado exitosamente en ', file_name1]);



% Save Index
% ========================================================
Wavelength_col = lambda_nm(:);
n_sim_col = n_final(:);
k_sim_col = k_final(:);

file_name2 = ['Simulation\', sample, 'nk_sim_.txt'];
datos_finales2 = [1240./Wavelength_col, n_sim_col, k_sim_col];

% CORRECCIÓN 1: Usamos file_name2 para no borrar el archivo anterior
fileID = fopen(file_name2, 'w');

% CORRECCIÓN 2: Reemplazamos el fprintf vacío por los encabezados de columna
%fprintf(fileID, 'Wavelength(nm)\tn\tk\n');

% Guardamos la matriz de datos
fprintf(fileID, '%e\t%f\t%f\n', datos_finales2');
fclose(fileID);

% CORRECCIÓN 3: Actualizamos el mensaje para que muestre el archivo correcto
disp(['Los datos n y k se han guardado exitosamente en ', file_name2]);






% Plot Reflectance
figure(1); clf;
plot(lambda_nm, y_exp, 'k-', 'LineWidth', 1.5, 'DisplayName', ['Experiment (', sample, ')']); hold on;
plot(lambda_nm, R_simulated, 'r--', 'LineWidth', 2, 'DisplayName', 'Adachi + Fresnel Sim');
xlabel('Wavelength (nm)'); ylabel('Reflectance');
title('Reflectance Match'); legend(); grid on; xlim([187, 400]);

% Plot n and k
figure(2); clf;
[ax, h1, h2] = plotyy(lambda_nm, n_final, lambda_nm, k_final);
set(h1, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b');
set(h2, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r');
ylabel(ax(1), 'Refractive Index (n)'); ylabel(ax(2), 'Extinction Coefficient (k)');
xlabel('Wavelength (nm)'); title('Extracted Optical Constants for this Film'); grid on;


