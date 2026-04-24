clear all;
hold off

xc=0.2;
Eg_AlN=6.2;
Eg_InN=0.7;
Eg_b=2.75;
Eg_AlInN=xc*Eg_InN+(1-xc)*Eg_AlN-xc*(1-xc)*Eg_b;
BC=0.7*(Eg_AlN-Eg_AlInN);
BV=0.3*(Eg_AlN-Eg_AlInN);
buffer=2000;
Lw=100;
Lb=250;
PSP_AlN=-0.081;
eps_AlN=8.5;
eps_InN=15.3;
eps_AlInN=xc*eps_InN+(1-xc)*eps_AlN;


%las polarizaciones en cada tegio
p(1)=PSP_AlN;
l(1)=buffer+Lb;
eps(1)=eps_AlN;
p(2)=polt_AlInN(0.20);
l(2)=Lw;
xc=0.2;
eps(2)=xc*eps_InN+(1-xc)*eps_AlN;
p(3)=PSP_AlN;
l(3)=Lb;
eps(3)=eps_AlN;
p(4)=polt_AlInN(0.21);
l(4)=Lw;
xc=0.21;
eps(4)=xc*eps_InN+(1-xc)*eps_AlN;
p(5)=PSP_AlN;
l(5)=Lb;
eps(5)=eps_AlN;
p(6)=polt_AlInN(0.22);
l(6)=Lw;
xc=0.22;
eps(6)=xc*eps_InN+(1-xc)*eps_AlN;
p(7)=PSP_AlN;
l(7)=Lb;
eps(7)=eps_AlN;

%masas efectivas

norm=0;
for i=1:7
norm=norm+l(i)/eps(i);
end
for i=1:7
    F(i)=0;
    for j=1:7
        F(i)=F(i)+(p(j)-p(i))*l(j)/(eps(j)*8.85E-12)/(eps(i)*norm);
    end
end
ver=0;
for i=1:7
    ver=ver+F(i)*l(i);
end
%construccion del potencia

%buffer
cont=0;
    [mel,mhh]=meffec(0);
[ dEc,dEv,gap ] = offsets( 0 );
for i=1:buffer+Lb
    cont=cont+1;

%el potencial electrico
veoff(cont)=dEc;
vhhoff(cont)=-1*gap+dEv;
if i==1
v(cont)=-F(1)*1E-10;
else
v(cont)=v(cont-1)-F(1)*1E-10;
end
meel(cont)=mel;
mehh(cont)=mhh;
end
fprintf('el contador es %d \n',cont)
%primer pozo
  [mel,mhh]=meffec(0.2);
[ dEc,dEv,gap ] = offsets( 0.2 );
for i=1:Lw
    cont=cont+1;

  veoff(cont)=-1*dEc;
vhhoff(cont)=-1*6.2+dEv;
%el potencial electrico
v(cont)=v(cont-1)-F(2)*1E-10;


meel(cont)=mel;
mehh(cont)=mhh;
end

%primer barrera
   [mel,mhh]=meffec(0);
[ dEc,dEv,gap ] = offsets( 0);
for i=1:Lb
    cont=cont+1;

 
veoff(cont)=-1*dEc;
vhhoff(cont)=-1*gap+dEv;
%el potencial electrico
v(cont)=v(cont-1)-F(3)*1E-10;


meel(cont)=mel;
mehh(cont)=mhh;
end
%segudo pozo
 [mel,mhh]=meffec(0.21);
[ dEc,dEv,gap ] = offsets( 0.21 );
for i=1:Lw
    cont=cont+1;

   
veoff(cont)=-1*dEc;
vhhoff(cont)=-6.2+dEv;
%el potencial electrico
v(cont)=v(cont-1)-F(4)*1E-10;


meel(cont)=mel;
mehh(cont)=mhh;
end

%segunda barrera
[mel,mhh]=meffec(0);
[ dEc,dEv,gap ] = offsets( 0);
for i=1:Lb
    cont=cont+1;

    
veoff(cont)=-1*dEc;
vhhoff(cont)=-1*gap+dEv;
%el potencial electrico
v(cont)=v(cont-1)-F(5)*1E-10;


meel(cont)=mel;
mehh(cont)=mhh;
end

%tercer pozo
 [mel,mhh]=meffec(0.22);
[ dEc,dEv,gap ] = offsets( 0.22 );
for i=1:Lw
    cont=cont+1;

   
veoff(cont)=-1*dEc;
vhhoff(cont)=-1*6.2+dEv;
%el potencial electrico
v(cont)=v(cont-1)-F(6)*1E-10;


meel(cont)=mel;
mehh(cont)=mhh;
end


%ultima barrera
[mel,mhh]=meffec(0);
[ dEc,dEv,gap ] = offsets( 0);
for i=1:Lb
    cont=cont+1;

    
veoff(cont)=-1*dEc;
vhhoff(cont)=-1*gap+dEv;
%el potencial electrico
v(cont)=v(cont-1)-F(7)*1E-10;


meel(cont)=mel;
mehh(cont)=mhh;
end

plot(veoff-v)
hold on
plot(vhhoff-v)
[e,ff]=matricial(meel',veoff'-v',1);
[ehh,ffhh]=matricial(mehh',-vhhoff'+v',1);

%calculo de la energia del primer pozo
 [melw,mhhw]=meffec(0.2);
[melb,mhhb]=meffec(0);
[ dEc,dEv,gap ] = offsets( 0.2 );
%constantes que requiero
hbar=6.582E-16;  %eV-s
m0=0.511E6/(3E8)^2;  %eV/(m/s)^2
alpaw=(2*melw*m0*abs(F(2))/hbar^2)^(1/3);
alpab=(2*melb*m0*abs(F(3))/hbar^2)^(1/3);

alpawhh=(2*mhhw*m0*abs(F(2))/hbar^2)^(1/3);
alpabhh=(2*mhhb*m0*abs(F(3))/hbar^2)^(1/3);%construccion del potencia
for i=1:Lb
ws0(i)=dEc;
ws1(i)=abs(F(2))*(100E-10)-(250-i)*abs(F(1))*1E-10;
ws2(i)=dEv;
ws3(i)=-abs(F(2))*(100E-10)+(250-i)*abs(F(1))*1E-10;

end
for i=Lb+1:Lb+Lw
ws0(i)=0;
ws1(i)=abs(F(2))*(1E-10)*(Lw+Lb-i);
ws2(i)=0;
ws3(i)=-abs(F(2))*(1E-10)*(Lw+Lb-i);

end
for i=Lb+Lw+1:Lb+Lw+Lb
ws0(i)=dEc;
ws1(i)=abs(F(3))*(1E-10)*(i-Lw-Lb);
ws2(i)=dEv;
ws3(i)=-abs(F(3))*(1E-10)*(i-Lw-Lb);

end

for i=1:2000
s=i*.001;
gg(i,1)=s;
gg(i,2)=airy(0,alpaw*s/abs(F(2)))*airy(1,alpab*(dEc-s)/abs(F(3)))*alpab*melw/(melb*alpaw)+airy(1,alpaw*s/abs(F(2)))*airy(0,alpab*(dEc-s)/abs(F(3)));
end
plot(gg(:,1),gg(:,2)*1E3)