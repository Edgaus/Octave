
%%% Parametros de la simulacion
d=200; %%Espesor de la capa de GaN
t0=45; %% angulo de incidencia respecto a la normal
n0=1; % la primer iterfaz es vacio



% importar los datos del indice de refraccion del Si y del GaN
filename = 'Si_n2.txt';
[B,~]=importdata(filename);
%B(:,1)=1000*B(:,1);

%filename = 'GaNH-n.txt';
filename = 'AlN-n.dat';

[AlNn, delimiterOut]=importdata(filename);
%GaNH-n(:,1)==1240./GaNH-n(:,1);

%filename='GaNH-ext.txt';
% filename='AlN-k.txt';
% [AlNk,delimiterOut] = importdata(filename);
%GaNH-ext(:,1)==1240./GaNH-ext(:,1);
A  =2

FP=@(x) A*sqrt(sqrt(x^2+1)+x)/sqrt(x^2+1);
FM=@(x) A*sqrt(sqrt(x^2+1)-x)/sqrt(x^2+1);
FPN=@(x) -A*sqrt(sqrt(x^2+1)+x)/sqrt(x^2+1);


gm=0.150;  %valor del ensanchamiento en eV
dW=5*(4.5E-4);%cambio del punto critico (dw/dT)*dT  dT=1;
dG=0.1*(15E-4); %cambio del ensanchamiento  (dG/dt)*dT   dT=1;
lg=207.5;
rs=0.1;
for kk=1:201;
    l=199.9+kk*rs;
nn0(kk)=interp1(AlNn(:,1),AlNn(:,2),l+204-lg);
kk0(kk)=interp1(AlNn(:,1),AlNn(:,3),l+204-lg);
nsi(kk)=interp1(B(:,1),B(:,2),l);
wl(kk)=l;
%n1=2.2;recorr
%n1=interp1(GaNH-n(:,1),GaNH-n(:,2),l)+10E-7;  % aqui se agrega el cambio de
%deltan

end
n0s=smooth(nn0,0.1)-0.2;
k0s=smooth(kk0,0.1);
nsis=smooth(nsi,0.1);
%calculo de la parte real y compleja del indice de refraccion

for kk=1:201;
    l=199.9+kk*rs;
er(kk)=n0s(kk)^2-k0s(kk)^2;
ei(kk)=2*n0s(kk)*k0s(kk);

end

%calculo en el cambio del indice de la parte real y e imaginaria de la funcion dielectrica
for kk=1:201;
    l=199.9+kk*rs;
%primero calculo el valor de x
 xc=(1240/l-1240/lg)/gm;
 der(kk)=FM(xc)*dW-FPN(xc)*dG;
 dei(kk)=FP(xc)*dW+FM(xc)*dG;



end
%caclulo de los cambios del indice de refraccion y coeficiente de
%refraccion
for kk=1:201;
    l=199.9+kk*rs;
%primero calculo el valor de x
 dn(kk)=1/(  4*n0s(kk)    )    *(          (er(kk)/sqrt(er(kk)^2+ei(kk)^2)+1)*der(kk)+ei(kk)/sqrt(er(kk)^2+ei(kk)^2)*dei(kk)  );







 if l<lg
 dk(kk)=1/(4*k0s(kk))*((er(kk)/sqrt(er(kk)^2+ei(kk)^2)-1)*der(kk)+ei(kk)/sqrt(er(kk)^2+ei(kk)^2)*dei(kk));
 else
     dk(kk)=0.5/sqrt(er(kk))*dei(kk)-0.25*ei(kk)/(er(kk))^(3/2)*der(kk);
 end
end




%calculo el valor de la reflectancia previa
for kk=1:201;
    l=wl(kk);
%primero calculo el valor de x
n1=n0s(kk)+dn(kk);
%n1=interp1(GaNn(:,1),GaNn(:,2),l)+10E-7;  % aqui se agrega el cambio de
%deltan
n1c=k0s(kk)+dk(kk);
n2=nsis(kk);






%Se calcula los angulos usando ley de Snell
t1=asind(n0*sind(t0)/n1);
t2=asind(n1*sind(t1)/n2);
k1=2*pi*n1/l;
k1c=2*pi*abs(n1c)/l;

%Polarizacion tipo s
rs01=(sind(t1)*cosd(t0)-sind(t0)*cosd(t1))/(sind(t1)*cosd(t0)+sind(t0)*cosd(t1));
rs10=-rs01;
rs12=(sind(t2)*cosd(t1)-sind(t1)*cosd(t2))/(sind(t2)*cosd(t1)+sind(t1)*cosd(t2));
ts01=2*sind(t1)*cosd(t0)/(sind(t1)*cosd(t0)+sind(t0)*cosd(t1));
ts10=2*sind(t0)*cosd(t1)/(sind(t1)*cosd(t0)+sind(t0)*cosd(t1));
rstot=rs01+(ts01*ts10*rs12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)))/(1-rs10*rs12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)));

%Polarizacion tipo p
rp01=(sind(t1)*cosd(t1)-sind(t0)*cosd(t0))/(sind(t1)*cosd(t1)+sind(t0)*cosd(t0));
tp01=2*sind(t1)*cosd(t0)/(sind(t1)*cosd(t1)+sind(t0)*cosd(t0));
rp10=-rp01;
tp10=2*sind(t0)*cosd(t1)/(sind(t0)*cosd(t0)+sind(t1)*cosd(t1));
rp12=(sind(t2)*cosd(t2)-sind(t1)*cosd(t1))/(sind(t2)*cosd(t2)+sind(t1)*cosd(t1));
rptot=rp01+(tp01*tp10*rp12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)))/(1-rp10*rp12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)));

aa(kk,1)=l; %longitud de onda
aa(kk,2)=abs(rstot).*abs(rstot);  %magnitud de la reflectancia tipo s
aa(kk,3)=angle(rstot);            %fase reflectancia tipo s
aa(kk,4)=abs(rptot).*abs(rptot);  %magnitud de la reflectancia tipo p
aa(kk,5)=angle(rptot);            %fase de la reflectancia
aa(kk,6)=0.5*(aa(kk,2)+aa(kk,4));  %reflectancia para luz no polarizada
end

%calculo el valor de la reflectancia previa
for kk=1:201;
    l=wl(kk);
%primero calculo el valor de x
n1=n0s(kk);
%n1=interp1(GaNn(:,1),GaNn(:,2),l)+10E-7;  % aqui se agrega el cambio de
%deltan
n1c=k0s(kk);
n2=nsis(kk);

%Se calcula los angulos usando ley de Snell
t1=asind(n0*sind(t0)/n1);
t2=asind(n1*sind(t1)/n2);
k1=2*pi*n1/l;
k1c=2*pi*abs(n1c)/l;

%Polarizacion tipo s
rs01=(sind(t1)*cosd(t0)-sind(t0)*cosd(t1))/(sind(t1)*cosd(t0)+sind(t0)*cosd(t1));
rs10=-rs01;
rs12=(sind(t2)*cosd(t1)-sind(t1)*cosd(t2))/(sind(t2)*cosd(t1)+sind(t1)*cosd(t2));
ts01=2*sind(t1)*cosd(t0)/(sind(t1)*cosd(t0)+sind(t0)*cosd(t1));
ts10=2*sind(t0)*cosd(t1)/(sind(t1)*cosd(t0)+sind(t0)*cosd(t1));
rstot=rs01+(ts01*ts10*rs12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)))/(1-rs10*rs12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)));

%Polarizacion tipo p
rp01=(sind(t1)*cosd(t1)-sind(t0)*cosd(t0))/(sind(t1)*cosd(t1)+sind(t0)*cosd(t0));
tp01=2*sind(t1)*cosd(t0)/(sind(t1)*cosd(t1)+sind(t0)*cosd(t0));
rp10=-rp01;
tp10=2*sind(t0)*cosd(t1)/(sind(t0)*cosd(t0)+sind(t1)*cosd(t1));
rp12=(sind(t2)*cosd(t2)-sind(t1)*cosd(t1))/(sind(t2)*cosd(t2)+sind(t1)*cosd(t1));
rptot=rp01+(tp01*tp10*rp12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)))/(1-rp10*rp12*exp(2i*k1*d*cosd(t1)-2*k1c*d*cosd(t1)));

bb(kk,1)=l; %longitud de onda
bb(kk,2)=abs(rstot).*abs(rstot);  %magnitud de la reflectancia tipo s
bb(kk,3)=angle(rstot);            %fase reflectancia tipo s
bb(kk,4)=abs(rptot).*abs(rptot);  %magnitud de la reflectancia tipo p
bb(kk,5)=angle(rptot);            %fase de la reflectancia
bb(kk,6)=0.5*(bb(kk,2)+bb(kk,4));  %reflectancia para luz no polarizada
end

%figure('Color', 'w', 'Position', [100 100 800 400]);
title('My Plot Title');
subplot(3,2,1)
plot(1240./bb(:,1),bb(:,6))
legend('Reflectancia','Location','SW')
subplot(3,2,2)
plot(1240./bb(:,1),dn)
legend('dn','Location','NW')
subplot(3,2,4)
plot(1240./bb(:,1),dk)
legend('dk','Location','NW')
subplot(3,2,3)
plot(1240./bb(:,1),(aa(:,6)-bb(:,6)));
legend('dR','Location','SW')
subplot(3,2,5)
plot(1240./bb(:,1),2*(aa(:,6)-bb(:,6))./(aa(:,6)+bb(:,6)+0.1))
legend('dR/R','Location','SW')
subplot(3,2,6)
plot(1240./bb(:,1),n0s);
legend('n(w)','Location','NW');


%AQUI PRETENDO EXPORTAR LOS ARCHIVOS
res(:,1)=1240./bb(:,1);
res(:,2)=n0s;
res(:,3)=k0s;
res(:,4)=dn;
res(:,5)=dk;
res(:,6)=aa(:,6);
res(:,7)=bb(:,6);

% fileID = fopen('res.dat','w');
%  fprintf(fileID,'%f %f %f %f %f %f %f\n',res');
%   fclose(fileID);
% %
