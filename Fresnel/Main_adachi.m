pkg load signal;
clear all;

% 1. Define Energy Array
lambda_nm = linspace(190,400, 210);
Energy_eV = 1240 ./ lambda_nm;


sample = 'm668';
file_name = ['Experimental\', sample, '.txt'];
data = load(file_name);
x_full = data(:, 1); % Wavelength in nm
y_full = data(:, 2); % Reflectance (0 to 1)

% Use the thickness optimization previous codes

if sample == 'm668'
  thickness  =  162.31601;
elseif sample =='m669'
  thickness  = 211.43622;
elseif sample =='m669'
  thickness  = 193.77740;
elseif sample =='m669'
  thickness  = 204.78842;
elseif sample =='m669'
  thickness  = 176.29674;
else
  thickness  = 0;
end

% Guesses: [E0, Gamma0, A0, E1, Gamma1, B1, eps_inf]
x_guess =   [5.8, 0.025, 35.0,8.70,0.60, 2.0, 1.5]

% Optimizer options (MaxFunEvals allows it to run longer to find the perfect fit)
options = optimset('Display', 'iter', 'MaxIter', 2000, 'MaxFunEvals', 4000, 'TolX', 1e-4);

disp('Starting Master Optical Optimization...');

% Run fminsearch!
[best_params, final_error] = fminsearch(@(x_guess) Fit_adachi_best(x_guess, Energy_eV, y_exp), x_guess, options);



% Unpack the winning parameters
best_E0  = best_params(1);
best_Gamma0     = best_params(2);
best_A0     = best_params(3);
best_E1    = best_params(4);
best_Gamma1    = best_params(5);
best_B1    = best_params(6);
best_eps_inf   = best_params(7);





fprintf('\n=== OPTIMIZATION COMPLETE ===\n');
fprintf('E0:  %.2f nm\n', best_E0 );
fprintf('A0:         %.2f\n', best_A0);
fprintf('Gamma0:         %.2f\n', best_Gamma0);
fprintf('E1:         %.2f\n', best_E1);
fprintf('Gamma1:    %.2f\n', best_Gamma1);
fprintf('eps_B1:    %.2f\n', best_B1);
fprintf('eps_inf:    %.2f\n', best_eps_inf);








% ========================================================
% 3. PLOT THE FINAL MATCH
% ========================================================
% Run the target function one last time to extract the simulated curve
[~, R_simulated, n_final, k_final] = master_target_function(best_params, lambda_nm, Energy_eV, y_exp);

figure(1); clf;
plot(lambda_nm, y_exp, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Experiment (m668)'); hold on;
plot(lambda_nm, R_simulated, 'r--', 'LineWidth', 2, 'DisplayName', 'Adachi + Fresnel Sim');
xlabel('Wavelength (nm)'); ylabel('Reflectance');
title('Reflectance Match'); legend(); grid on; xlim([225, 390]);

% Plot the extracted n and k!
figure(2); clf;
[ax, h1, h2] = plotyy(lambda_nm, n_final, lambda_nm, k_final);
set(h1, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b');
set(h2, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r');
ylabel(ax(1), 'Refractive Index (n)'); ylabel(ax(2), 'Extinction Coefficient (k)');
xlabel('Wavelength (nm)'); title('Extracted Optical Constants for this Film'); grid on;
