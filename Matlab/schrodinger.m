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
clf; clc; clear;
syms V(t) u0(x)



% Choose vEpsilon between 10^-3 and 10^-2
vEps = 10^(-2.5);

% Time Grid size
tM = 500;
% Time domain
t0 = 0;
t1 = 0.8;
% Timestep
tH = (t1 - t0)/tM;

% Space Grid size
xM = 500;
% Space domain
x0 = 0;
x1 = 1;
%'Space'step
xH = (x1 - x0)/xM;

% The potential in use is specific to the problem
V(x) = 10;

%initial value of u(x,T0), i.e the first U*
u0(x) = exp(-25*((x-0.5)^2))*exp(1i*(-1/5)*log(exp(5*(x - 0.5))+exp(-5*(x - 0.5)))/vEps);


u = Strang(vEps, tM, t0, t1, tH, xM, x0, x1, xH, V, u0);
%posDensity = PositionDensity(u)
%PlotAtConstantTime(135,posDensity,xM,xH,tH);
%PlotEntireDomain(posDensity);

%{
	Approximates a Schrodinger equation using strang
%}
function ret = Strang (vEps, tM, t0, t1, tH, xM, x0, x1, xH, V, u0)
	% populate solution space
	u = zeros(xM,tM);
	x = zeros(xM);

	for j = 1 : xM
		x(j) = x0 + (j-1)*xH;
		u(j,1) = u0(x(j));
	end
	
	% For each timestep (column of U)
	for n = 1 : tM-1
		uSt1 = zeros(xM);
		uSt2 = zeros(xM);
		% Solve the ODE section of the split
		for j = 1 : xM
			uSt1(j) = exp(-1i*V(x(j))*tH/(2*vEps))*u(j,n);
		end
		
		% Solve the fourier part of the split
		for j = 1 : xM			
			sumval = 0;
			for l = -xM/2 : (xM/2 - 1)
				% mu_l value
				mul = (2*pi*l)/(x1 - x0);
				
				% Calculate fourier coeffiec using ODE init cond
				uStFourier = 0;
				for j2 = 1 : xM
					uStFourier = uStFourier + uSt1(j2)*exp(-1i*mul*(x(j2) - x0));
				end
				
				sumval = sumval + exp(-1i*vEps*tH*(mul^2)/2) * uStFourier * exp(1i*mul*(x(j) - x0));
			end
			% New Ust from the fourier solution
			uSt2(j) = (1/xM)*sumval;
		end
		
		% Solution for the next timestep
		for j = 1 : xM
			[j,xM,n+1,tM]
			u(j, n+1) = exp(-1i*V(x(j))*tH/(2*vEps))*uSt2(j);
		end
	end		
	ret = u;
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

function ret = PlotAtConstantTime (timepoint, u, xM, xH, tH)
	t = timepoint*tH
	figure(1);
	vec = u(:,timepoint);
	
	xlabel('x')

	for n = 0 : xM - 1
		hold on
		xlabel('x')
		ylabel('y')
		
		plot(xH*(n),vec(n+1),'ob')
		
	end

	grid on
end

function ret = PlotEntireDomain (u,tH,xH)
	figure(2);
	surf(u)
	%view(2)
end



































