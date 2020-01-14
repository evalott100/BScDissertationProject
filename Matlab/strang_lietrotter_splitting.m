%{
	Eva Lott
	University of Dundee BSc project
	'Time splitting spectral methods for Schrodinger equations in the
	semiclassical'
	
	12/01/20
	strang_lietrotter_splitting.m
	Functions to demonstrate and compute Lie-Trotter and Strang splitting
	schemes
%}


%{
	Approximates function d/dx = Ax = (B + C)x using strang splitting
	delt    timestep
	B, C    matrices generated from split
	nPoints       Time to compute the approximated function at
%}

clc
clear


syms x(t) y(t)
A = [1 2; -1 1];
X = [x; y];
odeSystem = diff(X) == A*X;
ini = X(0) == [2; -1];
sols = dsolve(odeSystem,ini);

A = [1 2; -1 1];
C = [0 2; 0 1];
B = [1 0; -1 0];
X = [x; y];
ini =  [2; -1];
delt = 1;
T0 = 0;
T1 = 10;


ltAp = lieTrotter(delt, B, C, ini, T0, T1);

clf
fplot(sols.x,'b')
hold on
fplot(sols.y,'r')
xlabel('t')

n=1
while (T1 - T0)/delt > n*delt
	hold on
	var = ltAp(:,n)
	plot(n*delt,var(1),'ob')
	hold on
	plot(n*delt,var(2),'or')
	n = n+1
end

grid on
xlim([T0-2 T1])
ylim([-250, 250])

function ret = lieTrotter (delt, B, C, ini, T0, T1)
	numPoints = (T1 - T0)/delt;
	
	n = 1;
	x_ini(:,n) = ini;

	while numPoints >= n*delt
		x(:,n) = x_ini(:,n);
		x(:,n+1) = expm(B*delt)*x(:,n);
		y(:,n) = x(:,n+1);
		y(:,n+1) = expm(C*delt)*y(:,n);
		x_ini(:,n+1) = y(:,n+1);
		n = n+1;
	end
	ret = x;
end

%{
	Approximates function d/dx = Ax = (B + C)x using strang splitting
	delt    timestep
	B, C    matrices generated from split
	nPoints       Time to compute the approximated function at
%}
function ret = strang (delt, nPoints, B, C)
end