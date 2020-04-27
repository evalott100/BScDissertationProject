%{
	Eva Lott
	University of Dundee BSc project
	'Time splitting spectral methods for Schrodinger equations in the
	semiclassical'
	
	21/04/20
	Solve the heat equation at a specific time, using the method from
    "A brief introduction to pseudo-spectral methods:application to
    diffusion problems"
	https://cel.archives-ouvertes.fr/cel-01256472v2/document
%}

l  = 1.0;
N  = 256;
dx = 2*l/N;
x  = (1-N/2:N/2)*dx;
nu = 0.01;
T  = 5.0;
dk = pi/l;
k  = [0:N/2 1-N/2:-1]*dk;
k2 = k.^2;
u0 = sech(10.0*x).^2;
u0_hat = fft(u0);
T  = real(ifft(exp(-nu*k2*T).*u0_hat));
plot(x,T)
xlabel('x')


