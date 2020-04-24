%{
	Eva Lott
	University of Dundee BSc project
	'Time splitting spectral methods for Schrodinger equations in the
	semiclassical reigime'
	
	18/04/20
	Aprroximate splitted Schrodinger equation with slow Fourier and ODE
	solutions
%}
clf; clc; %clear;


% TODO: Make space step proportional to eps
% Choose vEpsilon between 10^-3 and 10^-2
vEps = 10^-3;

% Time Grid size
tM = 500;
% Time domain
t0 = 0;
t1 = 1;
% Timestep
tH = (t1 - t0)/tM;

% Space Grid size
xM =100;
% Space domain
x0 = 0;
x1 = 1;
%'Space'step
xH = (x1 - x0)/xM;

% The potential in use is specific to the problem
syms V(x)
V(x) = (x^2)/2;
%V(x) = 10
% Evaluate the potential at each space step
Vx = zeros(xM);
for j = 1 : xM
	Vx(j) = V(x0 + (j-1)*xH);
end

%initial value of u(x,T0), i.e the first U*
syms u0(x)
u0(x) = exp(-25*(x - 1/2)^2)*exp(1i/vEps * (1 + x));
%u0(x) = exp(-25*((x-0.5)^2))*exp(1i*(-1/5)*log(exp(5*(x - 0.5))+exp(-5*(x - 0.5)))/vEps);
U0 = zeros(xM);
for j = 1 : xM
	U0(j) = u0(x0 + (j-1)*xH);
end


u = Strang(vEps, tM, t0, t1, tH, xM, x0, x1, xH, Vx, U0);
%u = LieTrotter(vEps, tM, t0, t1, tH, xM, x0, x1, xH, Vx, U0);
posDensity = PositionDensity(u);
%currDensity = CurrentDensity(u,xH, vEps);
PlotAtConstantTime(270,posDensity,xM,xH,tH);
%PlotAtConstantTime(270,currDensity,xM,xH,tH);
%PlotEntireDomain(posDensity);
%PlotEntireDomain(currDensity);


%{
	Approximates a Schrodinger equation using strang
%}
function ret = Strang (vEps, tM, t0, t1, tH, xM, x0, x1, xH, Vx, U0)
	% Populate solution space
	u = zeros(xM,tM);

	% Solution in discretised space for the first timestep t = t0
	for j = 1 : xM
		u(j,1) = U0(j);
	end
	
	% For each timestep (column of U)
	for n = 1 : tM-1
		uSt1 = zeros(xM);
		uSt2 = zeros(xM);
		% Solve the ODE section of the split
		for j = 1 : xM
			uSt1(j) = exp(-1i*Vx(j)*tH/(2*vEps))*u(j,n);
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
					uStFourier = uStFourier + uSt1(j2)*exp(-1i*mul*((j2-1)*xH));
				end 
				
				sumval = sumval + exp(-1i*vEps*tH*(mul^2)/2) * uStFourier * exp(1i*mul*((j-1)*xH));
			end
			% New Ust from the fourier solution
			uSt2(j) = (1/xM)*sumval;
		end
		
		% Solution for the next timestep
		for j = 1 : xM
			[j,xM,n+1,tM]
			u(j, n+1) = exp(-1i*Vx(j)*tH/(2*vEps))*uSt2(j);
		end
	end		
	ret = u;
end


%{
	Approximates a Schrodinger equation using lie trotter splitting
%}
function ret = LieTrotter (vEps, tM, t0, t1, tH, xM, x0, x1, xH, Vx, U0)
	% Populate solution space
	u = zeros(xM,tM);

	% Solution in discretised space for the first timestep t = t0
	for j = 1 : xM
		u(j,1) = U0(j);
	end
	
	% For each timestep (column of U)
	for n = 1 : tM-1
		uSt1 = zeros(xM);
		
		% Solve the fourier part of the split
		for j = 1 : xM			
			sumval = 0;
			for l = -xM/2 : (xM/2 - 1)
				% mu_l value
				mul = (2*pi*l)/(x1 - x0);
				
				% Calculate fourier coeffiec using init cond
				uStFourier = zeros(1,xM);
				for j2 = 1 : xM
					uStFourier(j2) = u(j2,n)*exp(-1i*mul*((j2-1)*xH));
				end

				sumval = sumval + exp(-1i*vEps*tH*(mul^2)/2) * sum(uStFourier) * exp(1i*mul*((j-1)*xH));
			end
			% Ust from the fourier solution
			uSt1(j) = (1/xM)*sumval;
		end
		
		% Solution for the next timestep
		for j = 1 : xM
			[j,xM,n+1,tM]
			u(j, n+1) = exp(-1i*Vx(j)*tH/vEps)*uSt1(j);
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
	expVal = zeros(rows, columns);
	for i = 1 :  rows
		for j = 1 : columns
			expVal(i,j) = norm(u(i,j),2)^2;
		end	
	end
	ret = expVal;
end


%{
	Evaluate discrete wavefunction for current density
%}
function ret = CurrentDensity(u, xH, vEps)
	[xM, tM] = size(u);
	currDens = zeros(xM, tM);
	
	
	%{
	% approximate the derivative of u every point
	for n = 1 :  tM
		for j = 2 : xM - 1
			derivu = (u(j+1,n) - u(j-1,n))/(2*xH);
	
			% calculate the current density at each point
			currDens(j,n) = vEps*imag(conj(u(j,n))*derivu);
		end
	end
	
	%}
	% For each timestep we calculate spatial dft, and differentiate idft
	for n = 1 : tM
		uFourier = zeros(1,xM);
		for k = 1 : xM
			for j = 1 : xM
				% The (k - 1 - xM/2) term is just so k ranges from -xM/2 ... xM/2 - 1
				uFourier(k) = uFourier(k) + u(j,n)*exp(-1i*pi*2*(k - 1 - xM/2)*(j-1)/xM);  
			end
		end
		
		for j = 1 : xM
			Du = 0;
			for k = 1 : xM
				Du = Du + (1i*2*pi*(k - 1 - xM/2)/xM)*uFourier(k)*exp(1i*2*pi*(k - 1 - xM/2)*(j-1)/xM);
			end
			%I arbitrarily comment out the 1/xM factor
			%Weirdly, this gives results that look similar to Crank-Nicholson
			%Du = (1/xM) * Du;
			currDens(j,n) = vEps*imag(conj(u(j,n))*Du);
		end
	end


	ret = currDens
end


function ret = PlotAtConstantTime (timepoint, u, xM, xH, tH)
	t = (timepoint)*tH
	figure(1);
	vec = u(:,timepoint);
	
	for n = 1 : xM - 2
	%for n = 0 : xM - 1
		hold on
		plot(xH*(n),vec(n+1),'or')
		
	end
	ylim([-2,2])
	xlabel('x')
	ylabel('position density')
	grid on
end


function ret = PlotEntireDomain (u,tH,xH)
	
	figure(2);


	surf(u)
	xlabel('timestep n  : (t_n = t_0 + n \times tH)')
	ylabel('spacestep m  : (x_m = x_0 + m \times xH)')
	zlabel('current density')
	
	%view(2)
end
