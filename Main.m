% Main
format longg
Pot = 0.23; % Value of the potential in the barrier to solve (eV), Pot=0 inside the well
widht_well = 56; % Width of the well (Armstrong)
widht_barrier = 12*56 ; % Width of barrier (Armstrong)
mass_well = 0.067; % Effective mass of the free electron in the well (GaAs)
mass_barrier = 0.092; % Effective mass of the hole in the barrier (GaAs)

mesh_well = 4;
mesh_barrier = 4;
points_well = widht_well/mesh_well;
points_barrier = widht_barrier/mesh_barrier;

x_left = -(widht_well+widht_barrier)/2 : mesh_barrier : -widht_well/2 ;
x_center = -widht_well/2: mesh_well : widht_well/2;
x_right = widht_well/2 : mesh_barrier : (widht_well+widht_barrier)/2;

x = unique([x_left, x_center, x_right]);

h = diff(x);

Masses = ones(1, points_barrier+points_well) * mass_barrier;

Potential = ones(1, points_barrier+points_well) * Pot;

for i = points_barrier/2+1: points_barrier/2+ points_well
  Masses(1,i) = mass_well;
  Potential(1,i) = 0;
end


 [Eigevalues, Eigenfunctions] = Finite_diferences_octave(Potential, Masses, h );
 [Energys, I] = sort(Eigevalues);
 [Eigen, I] = sort(Eigenfunctions);

 Eigen(1,:);
 Energys(1,:)
 Energys(2,:)

