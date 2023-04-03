function y_value = analyt_Q2()
syms y(t) a
eqn = diff(y,t) == -10e6*(y-t^2) + t;
cond = y(0) == 1;
S = simplify(dsolve(eqn,cond));

y_value = @(time_vector) subs(S,{t},{time_vector}); % true solution