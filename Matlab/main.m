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
syms C1 C2
expform = [C1, C2]*simplify(rewrite(expm(L*t), 'sincos'))






function ret = errorCalc ()
end