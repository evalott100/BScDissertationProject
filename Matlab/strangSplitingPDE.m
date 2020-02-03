%{
	Eva Lott
	University of Dundee BSc project
	'Time splitting spectral methods for Schrodinger equations in the
	semiclassical'
	
	12/01/20
	strangSplittingPDE.m
	Example taken from
	https://folk.ntnu.no/andreas/OpSplit/Chapter2/Example2_2.html

	Demonstrating a splitting method in a linear transport problem.
%}
clf
clc
clear

N     = 128;  % Number of grid points. NB, N must be even!
T     = 0.25; % Final time
dt    = T/6;  % Time step in splitting
x     = linspace(-pi,pi,N); y=x;
[X,Y] = meshgrid(x,y);
r1    = sqrt((X-1).^2+(Y-1).^2);
r2    = abs(X+1.5)+abs(Y+1.5);
u0    = max(0.0,1-0.5*r1) + 1.5*(r2<0.8);


surf(x,y,u0), axis([-pi pi -pi pi 0 1.5]),
shading interp; light, view(-12,24), title('Inital data')

Nt=ceil(T/dt); dt=T/Nt;
u=zeros(length(x),length(y),Nt+1);
u(:,:,1)=u0;
for i=1:Nt
	u1         = heat(u(:,:,i),dt,1);   % Diffusing in the x-direction
	u(:,:,i+1) = heat(u1,dt,2);         % Diffusing in the y-direction
end

surf(x,y,u(:,:,Nt+1)); axis([-pi pi -pi pi 0 1.5])
shading interp; light, view(-12,24), title('Solution')


function u=heat(u,t,dir)
%%
% Numerically solves the heat equation (for a time |t|)
% 
% $$u_t = u_{xx}$$
% 
% with initial data |u| and periodic boundary conditions by a spectral 
% differentiation method. |u| must be a matrix and |x| is a coordinate in
% the column direction if |dir=1| and in the row direction if |dir==2|.
% Default is |dir=1|. 



%% Initial setup
if (dir>2),
	error('Only 2D diffusions allowed');
end;
transpose=0;
S=size(u);
N=S(1);
if (N==1),
	u=u';
	transpose=1;
end;
if nargin<3,
	dir=1;
end;
if (dir==2)
	u=u';
	transpose=1;
end;
S=size(u);
N=S(1);
if (mod(N,2)~=0),
	error('Not even length');
end;
%%
% Using a spectral scheme to solve the 1D heat equation for 
% each column of u. 
n=[0:N/2-1 0 -N/2+1:-1]';
e=ones(1,S(2));
uhat=fft(u);
uhat=(exp(-t.*n.^2)*e).*uhat;
u=real(ifft(uhat));
if (transpose)
	u=u';
end;
end