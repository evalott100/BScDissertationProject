clc
clear


% https://uk.mathworks.com/help/symbolic/solve-a-system-of-differential-equations.html
syms x(t) y(t)
A = [1 2; -1 1];
X = [x; y];
odeSystem = diff(X) == A*X;
ini = X(0) == [2; -1];

sols = dsolve(odeSystem,ini);
sols.x
sols.y
T0 = 0;
T1 = 10;


clf
fplot(sols.x)
hold on
fplot(sols.y)
grid on
xlim([T0-2 T1])
ylim([-250, 250])
xlabel('t')