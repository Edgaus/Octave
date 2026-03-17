% Vamos a resolver el caso unidimensional de un átomo con tight binding
% Acorde a los calculos analiticos, uno lllega a la expresión del hamiltoniano general a primeros vecinos-
clear all


e11 = 12 % Energia del estado base para el atomo hidrogeno. Se puede considerar que ya asimilo el shift

t = 2  % factor de acoplamiento (traslape) entre entre el átomo n-ésimo con n-1 y n+1

n = 1000 % Número de celdas unitarias i.e. de átomos en nuestro cristal unidemsional


H0 = eye(n).*e11;  % Diagonal

upper_diag = [zeros(n,1), eye(n,n-1)] ;   % Superdiagonal
down_diag = [zeros(1,n) ; eye(n-1,n)];

hamiltonian = H0 -  t.*(down_diag +upper_diag);  % Creates upper diagonal matrix

[ V,D  ] = eig( hamiltonian);
eigenvalues = diag(D);

k = linspace(0,n, n);

plot(k, eigenvalues)
