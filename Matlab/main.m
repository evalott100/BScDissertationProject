clc
clear


% https://uk.mathworks.com/help/symbolic/solve-a-system-of-differential-equations.html
syms x(t) y(t)
L = [1 2; -1 1];
w = [x; y];
odeSystem = diff(w) == L*w;
sols = dsolve(odeSystem);
[sols.x, sols.y]


% Our solution with matrix exponential
syms C3 C4
expform = simplify(rewrite(expm(L*t), 'sincos'))*[C3; C4];
[expform(1), expform(2)]
