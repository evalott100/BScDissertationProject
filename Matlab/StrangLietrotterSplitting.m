%{
	Eva Lott
	University of Dundee BSc project
	'Time splitting spectral methods for Schrodinger equations in the
	semiclassical'
	
	12/01/20
	Functions to demonstrate and compute Lie-Trotter and Strang splitting
	schemes
%}
clf
clc
clear

%plotExample(0.1,'st')
%plotExample(0.5,'lt')
plotCommutativeExample(1)

%{
	Approximates function dx/dt = Ax = (B + C)x using Lie-Trotter splitting
	delt    timestep
	B, C    matrices generated from split
	ini     initial condition matrices
	T0      start point in the domain
	T1      end point in the domain
%}
function ret = lieTrotter (delt, B, C, ini, T0, T1)
	n = 0;
	x_ini(:,1) = ini;
	
	domain = (T1 - T0)/delt;
	while domain >= n*delt
		x(:,n+1) = x_ini(:,n+1);
		x(:,n+2) = expm(B*delt)*x(:,n+1);
		y(:,n+1) = x(:,n+2);
		y(:,n+2) = expm(C*delt)*y(:,n+1);
		x_ini(:,n+2) = y(:,n+2);
		n = n+1;
	end
	ret = x;
end

%{
	Approximates function dx/dt = Ax = (B + C)x using Strang splitting
	delt    timestep
	B, C    matrices generated from split
	ini     initial condition ma(delt, nPoints, B, C)trices
	T0      start point in the domain
	T1      end point in the domain
%}
function ret = strang (delt, B, C, ini, T0, T1)
	n = 0;
	x_ini(:,1) = ini;
	
	domain = (T1 - T0)/delt;
	while domain >= n*delt
		x(:,n+1) = x_ini(:,n+1);
		x(:,n+2) = expm(0.5*B*delt)*x(:,n+1);
		
		y(:,n+1) = x(:,n+2);
		y(:,n+2) = expm(C*delt)*y(:,n+1);
		
		x(:,n+2) = expm(0.5*B*delt)*y(:,n+2);
		
		x_ini(:,n+2) = x(:,n+2);
		n = n+1;
	end
	ret = x;
end

%{
	Plots lie trotter example with given timestep for approximation
%}
function plotExample (delt, method)
	syms x(t) y(t)
	A = [1 2; -1 1];
	X = [x; y];
	odeSystem = diff(X) == A*X;
	ini = X(0) == [2; -1];
	sols = dsolve(odeSystem,ini);

	
	
	A = [1 2; -1 1];
	C = [0 2; 0 1];
	B = [1 0; -1 0];
	expm(A*delt) - expm(B*delt)*expm(C*delt)
	
	X = [x; y];
	ini =  [2; -1];
	T0 = 0;
	T1 = 10;

	if strcmp(method,'lt')
		approx = lieTrotter(delt, B, C, ini, T0, T1);
		title('Lie-Trotter approximation')
	elseif strcmp(method,'st')
		approx = strang(delt, B, C, ini, T0, T1);
		title('Strang approximation')
	end
	
	fplot(sols.x,'b')
	hold on
	fplot(sols.y,'r')
	xlabel('t')

	n = 0;
	domain = (T1 - T0)/delt;
	while domain >= n*delt
		hold on
		var = approx(:,n+1);
		plot(n*delt,var(1),'ob')
		hold on
		plot(n*delt,var(2),'or')
		n = n+1;
	end

	grid on
	xlim([T0-2 T1])
	ylim([-250, 250])
end

function plotCommutativeExample (delt)
	syms x(t) y(t)
	A = [1 2; -1 1];
	X = [x; y];
	odeSystem = diff(X) == A*X;
	ini = X(0) == [2; -1];
	sols = dsolve(odeSystem,ini);

	A = [1 2; -1 1];
	C = [1 0; 0 1];
	B = [0 2; -1 0];
	X = [x; y];
	ini =  [2; -1];
	T0 = 0;
	T1 = 10;
	expm(A*delt) - expm(B*delt)*expm(C*delt)


	approx = lieTrotter(delt, B, C, ini, T0, T1);
	
	fplot(sols.x,'b')
	hold on
	fplot(sols.y,'r')
	xlabel('t')

	n = 0;
	domain = (T1 - T0)/delt;
	while domain >= n*delt
		hold on
		var = approx(:,n+1);
		plot(n*delt,var(1),'ob')
		hold on
		plot(n*delt,var(2),'or')
		n = n+1;
	end

	grid on
	xlim([T0-2 T1])
	ylim([-250, 250])
end
