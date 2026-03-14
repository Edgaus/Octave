clear
clearvars

sample = 'm668'
filename = ['Experimental\',sample,'_deltaRR.txt'];
% 1. Read the file into a single structure variable
rawData = importdata(filename, "\t", 1);

filename12 = 'Data\AlN_k2.txt'
k_raw = importdata(filename12, "\t", 1);
numMatrix12 = rawData.data;
k_AlN = numMatrix12(:,2)

% 2. Extract the numeric data matrix from the structure
numMatrix = rawData.data;

% 3. Assign each column to your specific variables
% The colon (:) means "take all rows" for that specific column
energy_x = numMatrix(:, 1);
R0       = numMatrix(:, 2);
Rt       = numMatrix(:, 3);
DeltaRR  = numMatrix(:, 4);


% importar los datos del indice de refraccion del Si y del AlN
filename = 'Si_n2.txt';
% NOTE: If importdata fails in Octave, change the next line to: B = load(filename);
[B, ~] = importdata(filename);



#filename = ['Simulation\', sample, 'nk_sim_.txt'];
filename = ['AlN-n.dat'];

% NOTE: If importdata fails here, change to: AlNn = dlmread(filename, '', 1, 0);
[AlNn, delimiterOut] = importdata(filename);

Eg = 6

gm = 0.12;        % valor del ensanchamiento en eV
dW = 1*(4.5E-4);   % cambio del punto critico (dw/dT)*dT  dT=1;
dG = 0.1*(15E-4);  % cambio del ensanchamiento  (dG/dt)*dT   dT=1;
lg = 1240/Eg;
rs = 0.1;

for kk = 1:1000
    l = 149.9 + kk*rs;
    nn0(kk) = interp1(AlNn(:,1), AlNn(:,2),l );
    kk0(kk) = interp1(AlNn(:,1), AlNn(:,3), l);
    nsi(kk) = interp1(B(:,1), B(:,2), l);
    wl(kk) = l;
end



% OCTAVE FIX: Replace MATLAB's smooth(x, 0.1) with movmean(x, 21)
span = 20;
n0s = movmean(nn0, span) ;
k0s = movmean(kk0, span);
nsis = movmean(nsi, span);




dW = 5*(4.5E-4);
dG=0.1*(15E-4); %cambio del ensanchamiento  (dG/dt)*dT   dT=1;



n = n0s;
k = k0s;
aln_0s = n-1i.*k

A = 1

FP = @(x) A*sqrt(sqrt(x.^2+1)+x)./sqrt(x.^2+1);
FM = @(x)  A*sqrt(sqrt(x.^2+1)-x)./sqrt(x.^2+1);
FPN = @(x) -A*sqrt(sqrt(x.^2+1)+x)./sqrt(x.^2+1);

xc =(1240./wl - 1240/lg)/gm;





der= FM(xc)*dW - FPN(xc)*dG;
dei= FP(xc)*dW + FM(xc)*dG;


dn_aln = ( n./ (  2*(n.^2 + k.^2) )      ).*der + ( k./(2*(n.^2 + k.^2)) ).*dei;

dk_aln =  (- k./ ( 2* (n.^2 + k.^2) )  ).*der + ( n./(2*(n.^2 + k.^2)) ).*dei;



% --- 4. Create the Bridge and Call the Modifier ---
n_aln_mod = dn_aln+n;
k_aln_mod = dk_aln+k;

aln_mod = n_aln_mod - 1i*k_aln_mod

plot(wl, dk_aln )


% 1. Convertir los arreglos en vectores columna (verticales) usando (:)
Wavelength_col =1240./ wl(:);
Extintion_col = k0s(:);
Extintion_col_mod = k_aln_mod(:);
Difference_col = 100.* dk_aln(:);


% 2. Unir ambas columnas en una sola matriz
datos_finales = [Wavelength_col, Extintion_col,Extintion_col_mod,  Difference_col ];

% 3. Abrir (o crear) el archivo de texto en modo escritura ('w')
% Puedes cambiar 'Resultados_Reflectancia.txt' por el nombre que prefieras

path = ['Simulation\extintion_Simulation.txt' ]

fileID = fopen(path, 'w');

% 4. Escribir los encabezados separados por una tabulación (\t)
fprintf(fileID, 'Wavelength\tComparation\n');

% 5. Escribir los datos numéricos

% porque Octave lee los datos columna por columna al escribirlos.
fprintf(fileID, '%e\t%f\t%f\t%f\n', datos_finales');

% 6. Cerrar el archivo (¡Muy importante para que los datos se guarden en el disco!)
fclose(fileID);



























#plot(1240./wl,  dk_aln, 1240./wl, kk0.*100)



%% Initial Conditions
R = zeros(1,length(wl));
R_per = zeros(1,length(wl));
DeltaR = zeros(1,length(wl));


d =[162.31601]*1e-9;     %m905

n0 = 1
theta_inc =45;


%% Initial Code


for m=1:length(wl)
      n_structure = [n0,   aln_0s(m) , nsi(m)];
      n_structure_mod = [n0,  aln_mod(m), nsi(m)];

      R(m) = ComplexFresnelMatrix(wl(m)*1e-9, n_structure, d, deg2rad(theta_inc), 'm' );
      R_mod(m) = ComplexFresnelMatrix(wl(m)*1e-9, n_structure_mod, d, deg2rad(theta_inc), 'm' );
      DeltaR(m) =2*(R_mod(m)-R(m)) /(R_mod(m) +R(m) +0.1) ;
end

#plot(1240./wl,  100.*DeltaR, energy_x, DeltaRR)
%% Guardar los resultados en un archivo de texto

% 1. Convertir los arreglos en vectores columna (verticales) usando (:)
Wavelength_col =1240./ wl(:);
Reflectance_col = R(:);
Reflectance_col_mod = R_mod(:);
Difference_col = 100.* DeltaR(:);


% 2. Unir ambas columnas en una sola matriz
datos_finales = [Wavelength_col, Reflectance_col  ];

% 3. Abrir (o crear) el archivo de texto en modo escritura ('w')
% Puedes cambiar 'Resultados_Reflectancia.txt' por el nombre que prefieras

path = ['Simulation\',sample,'_Simulation.txt' ]

fileID = fopen(path, 'w');

% 4. Escribir los encabezados separados por una tabulación (\t)
fprintf(fileID, 'Wavelength\tComparation\n');

% 5. Escribir los datos numéricos

% porque Octave lee los datos columna por columna al escribirlos.
fprintf(fileID, '%e\t%f\t%f\t%f\n', datos_finales');

% 6. Cerrar el archivo (¡Muy importante para que los datos se guarden en el disco!)
fclose(fileID);

