%% %%%%%%%%%%% AlInN Band Gap %%%%%%%%%%%%

% 1. Constants
AlN_gap = 6.2; 
InN_gap = 0.7; 
x_In = [0.2, 0.21, 0.22];


%2 Papers data
b0 = 14.3;  %eV
x0 = 0.01;
n = 4; %
s = 0.122;

%3 Functions
bowing_AlInN = bowing_function(x_In, b0, x0, n, s); 
vl_band_gap_AlInN = @(x) (AlN_gap.*(1-x) + InN_gap.*x) - ...
                   (bowing_function(x, b0, x0, n, s) .* x .* (1-x));

%4 Values
gap_AlInN = vl_band_gap_AlInN(x_In)

%% %%%%%%%%%%%% Campos Pz Ps%%%%%%%%%%%%%%%%%%


% AlN constants

Ps_AlN = -0.081;  % C/m2

a_AlN = 3.112; % Armostrongs
e31_AlN = -0.58; % C/m2 
e33_AlN = 1.55;  % C/m2

C33_AlN = 389; % GPa
C13_AlN = 99; % GPa

% InN constants

Ps_InN = -0.042;  % C/m2

a_InN = 3.533; % Armostrongs
e31_InN = -0.57; % C/m2 
e33_InN = 0.97;  % C/m2

C33_InN = 182; % GPa
C13_InN = 121; % GPa

% AlInN constants
vegard_law = @(x, AlN_const, InN_const)  (1-x).*AlN_const +InN_const.*x ; 

a_AlInN = vegard_law( x_In, a_AlN, a_InN  );

e31_AlInN = vegard_law( x_In, e31_AlN, e31_InN  );
e33_AlInN = vegard_law( x_In, e33_AlN, e33_InN  );

C33_AlInN = vegard_law( x_In, C33_AlN, C33_InN  );
C13_AlInN = vegard_law( x_In, C13_AlN, C13_InN  );

pref_pz = e31_AlInN - e33_AlInN.*( C13_AlInN./C33_AlInN );
strain = (a_AlInN - a_AlN)/a_AlN ; % Para capa AlN relajada y AlInN tensa

% Calculo campos polarización

Ps = vegard_law( x_In, Ps_AlN, Ps_InN  );

residual_stress = 9;  %  
Pz = -2*pref_pz.*strain*(residual_stress)./100;

fprintf(' \n Los campos Pizoelectricos :\n');
fprintf('%-12.5g ', Pz);

fprintf('\n \n Los campos Espontaneos:\n');
fprintf('%-12.5g ', Ps);

%%% Contruccion de las bandas de energía %%%%

thickness_well = 10e-9;
thickness_barrier = 25e-9;


conduction_band_edge_profile = @(x)    ( 1 - heaviside(x-225))*AlN_gap + ...
                    (heaviside(x-225) - heaviside(x-235))*gap_AlInN(1) + ...
                    ( heaviside(x-235) - heaviside(x-260))*AlN_gap + ...
                    (heaviside(x-260) - heaviside(x-270))*gap_AlInN(2) + ...
                    ( heaviside(x-270) - heaviside(x-295))*AlN_gap + ...
                    ( heaviside(x-295) - heaviside(x-305))*gap_AlInN(3) + ...
                    ( heaviside(x-305) - heaviside(x-331))*AlN_gap;

x = linspace(0,330, 330);
plot(x, band_pont(x))
