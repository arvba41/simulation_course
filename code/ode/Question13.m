%% Question 13
clc
clear all

syms y(t) a
eqn = diff(y,t) == -a*(y-cos(t));
cond = y(0) == 0;
S = simplify(dsolve(eqn,cond)); 

% considering a = 10
y_value = @(time_vector) subs(S,{a,t},{10,time_vector}); % true solution

% function definition: y'(t) = -10*y(t), y(0) = 1

ts = 0; % start time
tf = 6; % finish time

a = 10; % cosntantv


%% Solving of equations
h = 1e-2;
tt = ts:h:tf;
y0 = 0; % initital value

fy = @(t,y) -a*(y-cos(t));

% ode45
% tic; [tm,ym,yh1] = rk45(@(t,y) -a*(y-cos(t)),y0,[ts tf],2);toc;
% ode45
tic; [t1,y1] = ode45(@(t,y) fy(t,y),[ts tf],y0);toc;
% true solution
y_exact = y_value(tt);

figure(1)
clf

plot(t1,y1,tt,y_exact)
grid on
ylabel('$y(t)$')
xlabel('$t$')
legend('ode45-mine','ode45','y')
figtex(gca,2);
