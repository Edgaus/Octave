1; % <--- DO NOT REMOVE. Tells Octave this is a script with local functions.

pkg load optim;
clear all;

% ========================================================
% 1. LOAD YOUR EXPERIMENTAL n AND k DATA
% ========================================================
data_n = importdata('Data\AlN_n2.txt');
data_k = importdata('Data\AlN_k2.txt');

% Create the exact wavelength grid you will use in your reflection script
lambda_nm = linspace(187, 400, 500);
Energy_eV = 1240 ./ lambda_nm;

% Interpolate the data onto this common grid
n_exp = interp1(data_n(:,1), data_n(:,2), lambda_nm, 'pchip');
k_exp = interp1(data_k(:,1), data_k(:,2), lambda_nm, 'pchip');


% ========================================================
% 2. INITIALIZE THE OPTIMIZER
% ========================================================
% Guesses: [E0, Gamma0, A0, E1, Gamma1, B1, eps_inf, A_ex]
x_guess = [6.1, 0.15, 35.0, 8.70, 0.60, 2.0, 1.5, 1.5]; % A_ex is the 8th parameter

options = optimset('Display', 'iter', 'MaxIter', 3000, 'MaxFunEvals', 5000, 'TolX', 1e-4);

disp('Starting the Adachi Pre-Fit Optimization...');

% ========================================================
% 3. RUN FMINSEARCH
% ========================================================
[best_params, final_error] = fminsearch(@(x) target_function_nk(x, Energy_eV, n_exp, k_exp), x_guess, options);

% Unpack ALL 8 winning parameters correctly
best_E0       = best_params(1);
best_Gamma0   = best_params(2);
best_A0       = best_params(3);
best_E1       = best_params(4);
best_Gamma1   = best_params(5);
best_B1       = best_params(6);
best_eps_inf  = best_params(7);
best_A_ex     = best_params(8); % <--- Unpacked successfully!

fprintf('\n=== PRE-FIT OPTIMIZATION COMPLETE ===\n');
fprintf('Copy these values into your main reflection script x_guess!\n\n');
fprintf('E0:         %.4f eV\n', best_E0 );
fprintf('Gamma0:     %.4f\n', best_Gamma0);
fprintf('A0:         %.4f\n', best_A0);
fprintf('E1:         %.4f eV\n', best_E1);
fprintf('Gamma1:     %.4f\n', best_Gamma1);
fprintf('B1:         %.4f\n', best_B1);
fprintf('eps_inf:    %.4f\n', best_eps_inf);
fprintf('A_ex:       %.4f\n', best_A_ex); % <--- Printed successfully!


% ========================================================
% 4. PLOT THE FIT MATCH
% ========================================================
% ========================================================
% 4. PLOT THE FIT MATCH
% ========================================================
% Pass all 8 parameters into the math engine
[n_sim, k_sim] = adachi_ext(Energy_eV, best_E0, best_Gamma0, best_A0, best_E1, best_Gamma1, best_B1, best_eps_inf, best_A_ex);

figure(1); clf;

subplot(2,1,1);
plot(lambda_nm, n_exp, 'k-', 'LineWidth', 2, 'DisplayName', 'Target n Data'); hold on;
plot(lambda_nm, n_sim, 'b--', 'LineWidth', 2, 'DisplayName', 'Adachi Fit n');
ylabel('Refractive Index (n)'); title('Adachi Pre-Fit Match'); legend('Location', 'northeast'); grid on; xlim([187, 400]);

subplot(2,1,2);
plot(lambda_nm, k_exp, 'k-', 'LineWidth', 2, 'DisplayName', 'Target k Data'); hold on;
plot(lambda_nm, k_sim, 'r--', 'LineWidth', 2, 'DisplayName', 'Adachi Fit k');
xlabel('Wavelength (nm)'); ylabel('Extinction Coeff (k)'); legend('Location', 'northeast'); grid on; xlim([187, 400]);chi_ext(Energy_eV, best_E0, best_Gamma0, best_A0, best_E1, best_Gamma1, best_B1, best_eps_inf, best_A_ex);


% ========================================================
% 5. SAVE THE OPTIMIZED DATA TO TEXT FILE
% ========================================================
nk_final_data = [lambda_nm(:), n_sim(:), k_sim(:)];
fileID = fopen('Simulation\Optimized_nk_PreFit.txt', 'w');
fprintf(fileID, 'Wavelength(nm)\tn\tk\n');
fprintf(fileID, '%e\t%f\t%f\n', nk_final_data');
fclose(fileID);
disp('The optimized n and k curves have been successfully saved to Simulation\Optimized_nk_PreFit.txt');
