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
	T       Time to compute the approximated function at
%}
function ret = lieTrotter (delt, T, B, C)
	while 
end

%{
	Approximates function d/dx = Ax = (B + C)x using strang splitting
	delt    timestep
	B, C    matrices generated from split
	T       Time to compute the approximated function at
%}
function ret = strang (delt, T, B, C)
end