clc
clear


% https://uk.mathworks.com/help/symbolic/solve-a-system-of-differential-equations.html
syms x(t) y(t)
L = [1 2; -1 1];
w = [x; y];
odeSystem = diff(w) == L*w;
ini = w(0) == [2; -1];

sols = dsolve(odeSystem,ini);
sols.x
sols.y

clf
fplot(sols.x)
hold on
fplot(sols.y)
grid on
legend('sols.x','sols.y','Location','best')