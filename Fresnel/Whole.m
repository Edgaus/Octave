clearvars
clear
pkg load signal


sample = 'm691';
file_name = ['Experimental\', sample, '.txt'];
limit_1 = 250
limit_2 = 390


% Load Raw Data
data = load(file_name);
x_full = data(:, 1);
y_full = data(:, 2);

% --- 1. GET EXPERIMENTAL PEAKS ---
% This is done ONLY ONCE before the optimization starts
disp('Extracting Experimental Peaks...');
[x_max_exp, x_min_exp] = der_functi(sample,limit_1, limit_2);

disp('Experimental Maxima found at:'); disp(x_max_exp');
disp('Experimental Minima found at:'); disp(x_min_exp');

%% --- 2. LOAD OPTICAL CONSTANTS ---
index_refrac_aln = importdata('Data\AlN_n2.txt');
Wavelength_index_aln = index_refrac_aln(:, 1);
aln_index = index_refrac_aln(:, 2);

extintion_coeff = importdata('Data\AlN_k2.txt');
Wavelength_extin_aln = extintion_coeff(:, 1);
aln_extin = extintion_coeff(:, 2);

index_refrac_si = importdata('Data\Si_n2.txt');
Wavelength_index_si = index_refrac_si(:, 1);
si_index = index_refrac_si(:, 2);

extintion_coeff_si = importdata('Data\Si_k2.txt');
Wavelength_extin_si = extintion_coeff_si(:, 1);
si_extin = extintion_coeff_si(:, 2);

% Interpolate onto standard wavelength grid
lambda_nm = linspace(limit_1, limit_2, 750);  % [nm]

n_aln = interp1(Wavelength_index_aln, aln_index, lambda_nm, 'pchip') ...
     - 1i*interp1(Wavelength_extin_aln, aln_extin, lambda_nm, 'pchip');

n_si = interp1(Wavelength_index_si, si_index, lambda_nm, 'pchip') ...
     - 1i*interp1(Wavelength_extin_si, si_extin, lambda_nm, 'pchip');

%% --- 3. RUN OPTIMIZATION ---
min_thickness = 150;
max_thickness = 230;
options = optimset('Display', 'iter', 'TolX', 1e-9);

disp('Starting Peak-Matching Optimization...');

% Call the fminbnd using our custom peak-matching function
[best_thickness, final_error] = fminbnd(@(t) target_function(t, n_aln, n_si, lambda_nm, x_max_exp, x_min_exp), ...
                         min_thickness, max_thickness, options);

fprintf('\nSuccess! The optimal thickness is: %.5f nm\n', best_thickness);

%% --- 4. PLOT FINAL RESULTS ---
% Run the function one last time to get the Simulated Array for plotting
[~, R_best] = target_function(best_thickness, n_aln, n_si, lambda_nm, x_max_exp, x_min_exp);

figure(1); clf;
plot(x_full, y_full, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Experimental Raw');
hold on;
plot(lambda_nm, R_best, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Simulated (%.2f nm)', best_thickness));

% Plot true experimental peaks as markers
y_max_exp = interp1(x_full, y_full, x_max_exp);
y_min_exp = interp1(x_full, y_full, x_min_exp);
plot(x_max_exp, y_max_exp, 'b^', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Exp Max');
plot(x_min_exp, y_min_exp, 'gv', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Exp Min');

xlim([225, 390]);
xlabel('Wavelength (nm)'); ylabel('Reflectance');
title('Peak-to-Peak Optimization Match');
legend('Location', 'northeast');
grid on;
