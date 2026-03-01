pkg load signal;
clear all;

% 1. Define Energy Array
Energy = linspace(3, 9, 500);

% 2. Define Parameters
E0      = 6.20;
Gamma0  = 0.15;
A0      = 35.0;
A_ex    = 0.5;
Eb      = 0.060;
Gamma_ex= 0.08;
E1      = 8.70;
Gamma1  = 0.60;
B1      = 2.0;
eps_inf = 1.5;

% 3. Calculate n and k (Now calling 'adachi' to match your filename)
[n_full, k_full] = adachi(Energy, E0, Gamma0, A0, E1, Gamma1, B1, A_ex, Eb, Gamma_ex, eps_inf);

% 4. Plotting with Octave's 'plotyy'
figure(1); clf;

% plotyy creates the dual axes.
% ax holds the left/right axes. h1/h2 hold the line data.
[ax, h1, h2] = plotyy(Energy, n_full, Energy, k_full);

% Format the lines
set(h1, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b');
set(h2, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r');

% Format the axes colors
set(ax(1), 'ycolor', 'b');
set(ax(2), 'ycolor', 'r');

% Add Labels
ylabel(ax(1), 'Refractive Index (n)');
ylabel(ax(2), 'Extinction Coefficient (k)');
xlabel('Photon Energy (eV)');
title('Full Adachi Model (E0 + Exciton + E1) for AlN');
grid on;

xlim(ax(1), [3, 9]);
xlim(ax(2), [3, 9]);
