%{
	Eva Lott
	University of Dundee BSc project
	'Time splitting spectral methods for Schrodinger equations in the
	semiclassical reigime'
	
	03/02/20
	schrodinger.m
	Aprroximate splitted Schrodinger equation with Fourier and ODE
	solutions
%}
clc; clear;
syms V(t) U0(x)


% Choose epsilon between 10^-3 and 10^-2
eps = 10^(-2.5);

% Time Grid size
tM = 100;
% Time domain
t0 = 0;
t1 = 1;
% Timestep
tH = (t1 - t0)/tM;

% Space Grid size
xM = 100;
% Space domain
x0 = 0;
x1 = 1;
%'Space'step
xH = (x1 - x0)/xM;

% The potential in use is specific to the problem
V(x) = (x^2)/2;

%initial value of u(x,T0), i.e the first U*
u0(x) = exp(-25*(x - 0.5)^2)*exp(1i*(1+x)/eps);


u = Strang(eps, tM, t0, t1, tH, xM, x0, x1, xH, V, u0);
posDensity = PositionDensity(u)
vec = PlotAtConstantTime(50,posDensity,xM,xH)
PlotEntireDomain(posDensity)

%{
	Approximates a Schrodinger equation using strang
%}
function ret = Strang (eps, tM, t0, t1, tH, xM, x0, x1, xH, V, u0)
	% populate solution space
	u = zeros(xM,tM);
	for j = 0 : xM-1
		u(j+1,1) = u0(x0 + j*xH);
	end
	
	% For each timestep (column of U)
	for tIter = 0 : xM-1
		tIter
		uSt = zeros(xM);
		% For each space step (column of U)
		for sIter = 0 : xM-1
			sIter
			% Solve the ODE section of the split
			uSt(sIter+1) = exp(-1i*V(x0 + (sIter+1)*xH)*tH/(2*eps))*u(tIter+1,sIter+1);
			
			% Solve the fourier part of the split
			sumval = 0;
			for l = -xM/2 : (xM/2 - 1)
				mul = (2*pi*l)/(x1 - x0);
				
				% Calculate fourier coeffiec using ODE init cond
				uStFourier = 0;
				for j = 0 : xM - 1
					uStFourier = uStFourier + uSt(j+1)*exp(-1i*mul*(j+1)*xH);
				end
				
				sumval = sumval + exp(-1i*eps*tH*(mul^2)/2) * uStFourier * exp(1i*mul*((sIter+1)*xH));
			end
			% New Ust from the fourier solution
			uSt(sIter+1) = (1/xM)*sumval;
			u(tIter+1, sIter+1) = exp(-1i*V(sIter+1)*tH/(2*eps))*uSt(sIter+1);
		end			
	end
	ret = u
end

%{
	Take a matrix of complex values and apply the norm squared to every
	element inside
%}
function ret = PositionDensity(u)
	[rows, columns] = size(u);
	ret = zeros(rows, columns);
	for i = 1 :  rows
		for j = 1 : columns
			ret(i,j) = norm(u(i,j),2)^2;
		end	
	end

end

function ret = PlotAtConstantTime (timepoint, u, xM, xH)

	vec = u(:,timepoint);
	
	xlabel('x')

	for n = 0 : xM - 1
		hold on
		plot(xH*(n+1),vec(n+1),'ob')
	end

	grid on
	ret = vec
end

function ret = PlotEntireDomain (u,tH,xH)
	surf(u)
end



































