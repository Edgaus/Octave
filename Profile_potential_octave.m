%% Profile of Band Conduction

GaAs_gap = 1.424;  %eV

x_Al = 0.3;  %Al concentration in AlGaAs

if x_Al<0.45
  AlGaAs_gap = 1.424 + 1.247*x_Al;
else
  AlGaAs_gap = 1.900 + 0.125*x + 0.143*x^2;
end

AlGaAs_gap


width_layers = [ 150, 350, 400, 5000 ]

gap_material = [ GaAs_gap, AlGaAs_gap, AlGaAs_gap, GaAs_gap ]

Points = 100;

x = linspace(0, 5000, Points );


P = zeros(1, Points)

for i = 1: Points
  if x(i)<150
    P( 1,i) = GaAs_gap ;
  elseif x(i)< 350
    P(1,i) = AlGaAs_gap ;
  elseif x(i) < 400
    P(1,i) = AlGaAs_gap ;
  else
    P(1,i) = GaAs_gap;

  end
 end


 carrier_density = [ -1E18, -1E18, 1E13, 1E16 ];









