1; % Octave script flag

pkg load optim;
clear all;

% ========================================================
% 1. LOAD YOUR EXPERIMENTAL n AND k DATA
% ========================================================
data_n = importdata('Data\AlN_n2.txt');
data_k = importdata('Data\AlN_k2.txt');

lambda_nm = linspace(150, 1100, 951);
Energy_eV = 1240 ./ lambda_nm;

n_exp = data_n(:,2);
k_exp = data_k(:,2);

% ========================================================
% 2. INITIALIZE THE OPTIMIZER
% ========================================================
% Guesses: [E0, A0, E1, B1, eps_inf]
x_guess = [5.2, 30.0, 10, 2.0, 4.0,-2];

options = optimset('Display', 'iter', 'MaxIter', 4000, 'MaxFunEvals', 6000, 'TolX', 1e-5);

disp('Starting the Complete Adachi Fit (Spin-Orbit + Exciton)...');

% ========================================================
% 3. RUN FMINSEARCH
% ========================================================
[best_params, final_error] = fminsearch(@(x) target_function_nk(x, Energy_eV, n_exp, k_exp), x_guess, options);

% Unpack the 8 winning parameters
best_E0       = best_params(1);
best_A0       = best_params(2);
best_E1       = best_params(3);
best_B1       = best_params(4);
best_eps_inf  = best_params(5);


fprintf('\n=== PRE-FIT OPTIMIZATION COMPLETE ===\n');
fprintf('E0:         %.4f eV\n', best_E0 );
fprintf('A0:         %.4f\n', best_A0);
fprintf('E1:         %.4f eV\n', best_E1);
fprintf('B1:         %.4f\n', best_B1);
fprintf('eps_inf:    %.4f\n', best_eps_inf);


% ========================================================
% 4. PLOT THE FIT MATCH
% ========================================================
[n_sim, k_sim] = adachi_ext(Energy_eV, best_E0, best_A0, best_E1, best_B1, best_eps_inf);

figure(1); clf;

subplot(2,1,1);
plot(lambda_nm, n_exp, 'k-', 'LineWidth', 2, 'DisplayName', 'Target n Data'); hold on;
plot(lambda_nm, n_sim, 'b--', 'LineWidth', 2, 'DisplayName', 'Adachi Fit n');
ylabel('Refractive Index (n)'); title('Adachi Fit (Spin-Orbit + Exciton)'); legend('Location', 'northeast'); grid on; xlim([150, 1100]);

subplot(2,1,2);
plot(lambda_nm, k_exp, 'k-', 'LineWidth', 2, 'DisplayName', 'Target k Data'); hold on;
plot(lambda_nm, k_sim, 'r--', 'LineWidth', 2, 'DisplayName', 'Adachi Fit k');
xlabel('Wavelength (nm)'); ylabel('Extinction Coeff (k)'); legend('Location', 'northeast'); grid on; xlim([150, 1100]);

% ========================================================
% 5. SAVE THE OPTIMIZED DATA
% ========================================================
nk_final_data = [lambda_nm(:), n_sim(:), k_sim(:)];
fileID = fopen('Simulation\Optimized_nk_PreFit.txt', 'w');
fprintf(fileID, 'Wavelength(nm)\tn\tk\n');
fprintf(fileID, '%e\t%f\t%f\n', nk_final_data');
fclose(fileID);
disp('Saved successfully to Simulation\Optimized_nk_PreFit.txt');
