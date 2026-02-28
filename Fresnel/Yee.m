%%% Parametros de la simulacion
d = 200;  %% Espesor de la capa de GaN (Note: Using AlN data below)
t0 = 45;  %% angulo de incidencia respecto a la normal
n0 = 1;   % la primer interfaz es vacio

% importar los datos del indice de refraccion del Si y del AlN
filename = 'Si_n2.txt';
% NOTE: If importdata fails in Octave, change the next line to: B = load(filename);
[B, ~] = importdata(filename);

filename = 'AlN-n.dat';
% NOTE: If importdata fails here, change to: AlNn = dlmread(filename, '', 1, 0);
[AlNn, delimiterOut] = importdata(filename);

% Element-wise math added to anonymous functions
FP = @(x) sqrt(sqrt(x.^2+1)+x)./sqrt(x.^2+1);
FM = @(x) sqrt(sqrt(x.^2+1)-x)./sqrt(x.^2+1);
FPN = @(x) -1*sqrt(sqrt(x.^2+1)+x)./sqrt(x.^2+1);

gm = 0.150;        % valor del ensanchamiento en eV
dW = 5*(4.5E-4);   % cambio del punto critico (dw/dT)*dT  dT=1;
dG = 0.1*(15E-4);  % cambio del ensanchamiento  (dG/dt)*dT   dT=1;
lg = 207.5;
rs = 0.1;

for kk = 1:201
    l = 199.9 + kk*rs;
    nn0(kk) = interp1(AlNn(:,1), AlNn(:,2), l+204-lg);
    kk0(kk) = interp1(AlNn(:,1), AlNn(:,3), l+204-lg);
    nsi(kk) = interp1(B(:,1), B(:,2), l);
    wl(kk) = l;
end

% OCTAVE FIX: Replace MATLAB's smooth(x, 0.1) with movmean(x, 21)
span = 21;
n0s = movmean(nn0, span) - 0.2;
k0s = movmean(kk0, span);
nsis = movmean(nsi, span);

% calculo de la parte real y compleja del indice de refraccion
for kk = 1:201
    er(kk) = n0s(kk)^2 - k0s(kk)^2;
    ei(kk) = 2*n0s(kk)*k0s(kk);
end

% calculo en el cambio del indice de la parte real y e imaginaria
for kk = 1:201
    l = wl(kk);
    % primero calculo el valor de x
    xc = (1240/l - 1240/lg)/gm;
    der(kk) = FM(xc)*dW - FPN(xc)*dG;
    dei(kk) = FP(xc)*dW + FM(xc)*dG;
end

% caclulo de los cambios del indice de refraccion y coeficiente de extincion
for kk = 1:201
    l = wl(kk);
    dn(kk) = 1/( 4*n0s(kk) ) * ( (er(kk)/sqrt(er(kk)^2+ei(kk)^2)+1)*der(kk) + ei(kk)/sqrt(er(kk)^2+ei(kk)^2)*dei(kk) );

    if l < lg
        dk(kk) = 1/(4*k0s(kk)) * ((er(kk)/sqrt(er(kk)^2+ei(kk)^2)-1)*der(kk) + ei(kk)/sqrt(er(kk)^2+ei(kk)^2)*dei(kk));
    else
        dk(kk) = 0.5/sqrt(er(kk))*dei(kk) - 0.25*ei(kk)/(er(kk))^(3/2)*der(kk);
    end
end


new_n1 = zeros( 1, length( n0s  ))

% calculo el valor de la reflectancia previa (CON variacion)
for kk = 1:201
    l = wl(kk);
    n1 = n0s(kk) + dn(kk);

    n1c = k0s(kk) + dk(kk);
    new_n1(kk) = n1c;
    n2 = nsis(kk);

    % Se calcula los angulos usando ley de Snell
    t1 = asind(n0*sind(t0)/n1);
    t2 = asind(n1*sind(t1)/n2);
    k1 = 2*pi*n1/l;
    k1c = 2*pi*abs(n1c)/l;

    % Polarizacion tipo s
    rs01 = (sind(t1)*cosd(t0)-sind(t0)*cosd(t1))/(sind(t1)*cosd(t0)+sind(t0)*cosd(t1));
    rs10 = -rs01;
    rs12 = (sind(t2)*cosd(t1)-sind(t1)*cosd(t2))/(sind(t2)*cosd(t1)+sind(t1)*cosd(t2));
    ts01 = 2*sind(t1)*cosd(t0)/(sind(t1)*cosd(t0)+sind(t0)*cosd(t1));
    ts10 = 2*sind(t0)*cosd(t1)/(sind(t1)*cosd(t0)+sind(t0)*cosd(t1));
    rstot = rs01 + (ts01*ts10*rs12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)))/(1-rs10*rs12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)));

    % Polarizacion tipo p
    rp01 = (sind(t1)*cosd(t1)-sind(t0)*cosd(t0))/(sind(t1)*cosd(t1)+sind(t0)*cosd(t0));
    tp01 = 2*sind(t1)*cosd(t0)/(sind(t1)*cosd(t1)+sind(t0)*cosd(t0));
    rp10 = -rp01;
    tp10 = 2*sind(t0)*cosd(t1)/(sind(t0)*cosd(t0)+sind(t1)*cosd(t1));
    rp12 = (sind(t2)*cosd(t2)-sind(t1)*cosd(t1))/(sind(t2)*cosd(t2)+sind(t1)*cosd(t1));
    rptot = rp01 + (tp01*tp10*rp12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)))/(1-rp10*rp12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)));

    aa(kk,1) = l;
    aa(kk,2) = abs(rstot)*abs(rstot);
    aa(kk,3) = angle(rstot);
    aa(kk,4) = abs(rptot)*abs(rptot);
    aa(kk,5) = angle(rptot);
    aa(kk,6) = 0.5*(aa(kk,2)+aa(kk,4));
end




% calculo el valor de la reflectancia previa (SIN variacion)
for kk = 1:201
    l = wl(kk);
    n1 = n0s(kk);
    n1c = k0s(kk);
    n2 = nsis(kk);

    % Se calcula los angulos usando ley de Snell
    t1 = asind(n0*sind(t0)/n1);
    t2 = asind(n1*sind(t1)/n2);
    k1 = 2*pi*n1/l;
    k1c = 2*pi*abs(n1c)/l;

    % Polarizacion tipo s
    rs01 = (sind(t1)*cosd(t0)-sind(t0)*cosd(t1))/(sind(t1)*cosd(t0)+sind(t0)*cosd(t1));
    rs10 = -rs01;
    rs12 = (sind(t2)*cosd(t1)-sind(t1)*cosd(t2))/(sind(t2)*cosd(t1)+sind(t1)*cosd(t2));
    ts01 = 2*sind(t1)*cosd(t0)/(sind(t1)*cosd(t0)+sind(t0)*cosd(t1));
    ts10 = 2*sind(t0)*cosd(t1)/(sind(t1)*cosd(t0)+sind(t0)*cosd(t1));
    rstot = rs01 + (ts01*ts10*rs12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)))/(1-rs10*rs12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)));

    % Polarizacion tipo p
    rp01 = (sind(t1)*cosd(t1)-sind(t0)*cosd(t0))/(sind(t1)*cosd(t1)+sind(t0)*cosd(t0));
    tp01 = 2*sind(t1)*cosd(t0)/(sind(t1)*cosd(t1)+sind(t0)*cosd(t0));
    rp10 = -rp01;
    tp10 = 2*sind(t0)*cosd(t1)/(sind(t0)*cosd(t0)+sind(t1)*cosd(t1));
    rp12 = (sind(t2)*cosd(t2)-sind(t1)*cosd(t1))/(sind(t2)*cosd(t2)+sind(t1)*cosd(t1));
    rptot = rp01 + (tp01*tp10*rp12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)))/(1-rp10*rp12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)));

    bb(kk,1) = l;
    bb(kk,2) = abs(rstot)*abs(rstot);
    bb(kk,3) = angle(rstot);
    bb(kk,4) = abs(rptot)*abs(rptot);
    bb(kk,5) = angle(rptot);
    bb(kk,6) = 0.5*(bb(kk,2)+bb(kk,4));
end
plot(1240./bb(:,1),2*(aa(:,6)-bb(:,6))./(aa(:,6)+bb(:,6)+0.1))

%% Guardar los resultados en un archivo de texto

% 1. Convertir los arreglos en vectores columna (verticales) usando (:)
Wavelength_col = wl(:);
Reflectance_col = new_n1(:);




% 2. Unir ambas columnas en una sola matriz
datos_finales = [Wavelength_col, Reflectance_col];

% 3. Abrir (o crear) el archivo de texto en modo escritura ('w')
% Puedes cambiar 'Resultados_Reflectancia.txt' por el nombre que prefieras
fileID = fopen('Simulation\Yee_newk.txt', 'w');

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

