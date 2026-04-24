function [ dEc,dEv,gap ] = offsets( xc )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Eg_AlN=6.2;
Eg_InN=0.7;
Eg_b=2.75;
gap=xc*Eg_InN+(1-xc)*Eg_AlN-xc*(1-xc)*Eg_b;
dEc=0.7*(Eg_AlN-gap);
dEv=0.3*(Eg_AlN-gap);
end

