% Main_Fresnel
pkg load signal

sample = 'm672';
file_name = ['Experimental\', sample, '.txt'];
% Load the data
data = load(file_name);
x_full = data(:, 1);
y_full = data(:, 2);

% Assuming this function correctly gives you arrays of X coordinates
[x_max_exp, x_min_exp] = der_functi(sample);

%% AlN Data
index_refrac_aln =  importdata('Data\AlN_n2.txt');
Wavelength_index_aln = index_refrac_aln(:, 1); % Keep in nm!
aln_index = index_refrac_aln(:, 2);

extintion_coeff =  importdata('Data\AlN_k2.txt');
Wavelength_extin_aln = extintion_coeff(:, 1); % Keep in nm!
aln_extin = extintion_coeff(:, 2);

%% Si Data
index_refrac_si =  importdata('Data\Si_n2.txt');
Wavelength_index_si = index_refrac_si(:, 1); % REMOVED *1e-9 to match lambda_nm
si_index = index_refrac_si(:, 2);

extintion_coeff_si =  importdata('Data\Si_k2.txt');
Wavelength_extin_si = extintion_coeff_si(:, 1); % REMOVED *1e-9
si_extin = extintion_coeff_si(:, 2);

%%  Interpolate complex indices
lambda_nm = linspace(225, 390, 750);  % [nm]

% Interpolate AlN and Si
n_aln = interp1(Wavelength_index_aln, aln_index, lambda_nm, 'pchip') ...
     - 1i*interp1(Wavelength_extin_aln, aln_extin, lambda_nm, 'pchip');

n_si = interp1(Wavelength_index_si, si_index, lambda_nm, 'pchip') ...
     - 1i*interp1(Wavelength_extin_si, si_extin, lambda_nm, 'pchip');

% Safety Check!
if any(isnan(n_aln)) || any(isnan(n_si))
    error('CRITICAL ERROR: Interpolation produced NaNs. Check the wavelength ranges in your .txt files! They must cover 225 to 390 nm.');
end

%% Optimization
min_thickness = 150;
max_thickness = 250;

options = optimset('Display', 'iter', 'TolX', 1e-4);

disp('Starting Optimization...');
% fminbnd returns [best_x, best_error_value]
[best_thickness, final_error_val] = fminbnd(@(t) target_function(t, n_aln, n_si, lambda_nm, x_max_exp, x_min_exp), ...
                         min_thickness, max_thickness, options);

fprintf('\nSuccess! The optimal thickness is: %.5f nm\n', best_thickness);

%% Final Run & Plotting
% Call the function one last time with the best thickness to get the actual Simulated Array (R_opt)
[~, R_best] = target_function(best_thickness, n_aln, n_si, lambda_nm, x_max_exp, x_min_exp);

figure(2);
clf;
plot(x_full, y_full, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Experimental (m668)');
hold on;
plot(lambda_nm, R_best, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Simulated (%.1f nm)', best_thickness));

% Plot true peaks to visualize alignment
y_max_exp = interp1(x_full, y_full, x_max_exp);
y_min_exp = interp1(x_full, y_full, x_min_exp);
plot(x_max_exp, y_max_exp, 'b^', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Exp Max');
plot(x_min_exp, y_min_exp, 'gv', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Exp Min');

xlim([225, 390]);
xlabel('Wavelength (nm)');
ylabel('Reflectance');
title('Thin Film Transfer Matrix Optimization Match');
legend('Location', 'northeast');
grid on;
