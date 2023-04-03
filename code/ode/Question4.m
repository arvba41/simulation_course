clc;
clear all

h = 10e-2; % step size

%% Question 1
% y' = −10e6*y, y(0) = 1, t \in [0, 10e−6]
% Not enough samples to calculate the numerical solution 

%% Question 2
% y' = −10e6*(y-t^2) + t, y(0) = 1, t \in [0, 1]

syms y(t) a
eqn = diff(y,t) == -10e6*(y-t^2) + t;
cond = y(0) == 1;
S = simplify(dsolve(eqn,cond));

y_value = @(time_vector) subs(S,{t},{time_vector}); % true solution

ts = 0; tf = 1; % start and finist times, respectively 
y1_0 = 1; % initital condition

t1 = ts:h:tf; % time vector 

%%% Euler forward method
y1_f = euler_ft(@(y,t) -10e6*(y-t^2) + t,1,t1,h); 
%%% Euler backward method
y1_b = euler_bt(@(y,t) -10e6*(y-t^2) + t,1,t1,h); 
%%% Trapizoidal method
y1_fb = euler_fbt(@(y,t) -10e6*(y-t^2) + t,1,t1,h);

%%% Exact solution
y1_exact = y_value(t1);

%%%% Plottings
figure(1)
clf
tiledlayout(1,3)

nexttile;
plot(t1,[y1_f; y1_b; y1_fb],t1,y1_exact,'--r')
legend('analytical [$E_f$]','analytical [$E_b$]','analytical [$E_{fb}$]','excact')
grid on
figtex(gca,1)
ylabel('$$y$$')
xlabel('time [s]')
ylim([-10 10])

nexttile;
semilogy(t1,abs(y1_exact - [y1_f; y1_b; y1_fb]))
grid on
figtex(gca,[])
ylabel('$$|y_n - y(t)|$$')
xlabel('time [s]')

nexttile;
semilogy(t1,[y1_f; y1_b; y1_fb]./1)
grid on
figtex(gca,[])
ylabel('$$y_n/y_0$$')
xlabel('time [s]')

%% Question 3
% y' = -10e6*(y-sin(10^6*t)) + cos(10^6*t), y(0) = 1, t \in [0, 1]

syms y(t) a
eqn = diff(y,t) == -10e6*(y-sin(10^6*t)) + cos(10^6*t);
cond = y(0) == 1;
S = simplify(dsolve(eqn,cond));

y_value = @(time_vector) subs(S,{t},{time_vector}); % true solution

ts = 0; tf = 1; % start and finist times, respectively 
y2_0 = 1; % initital condition

t2 = ts:h:tf; % time vector 

%%% Euler forward method
y2_f = euler_ft(@(y,t) -10e6*(y-sin(10^6*t)) + cos(10^6*t) + t,y2_0,t2,h); 
%%% Euler backward method
y2_b = euler_bt(@(y,t) -10e6*(y-sin(10^6*t)) + cos(10^6*t) + t,y2_0,t2,h); 
%%% Trapizoidal method
y2_fb = euler_fbt(@(y,t) -10e6*(y-sin(10^6*t)) + cos(10^6*t) + t,y2_0,t2,h);

%%% Exact solution
y2_exact = y_value(t2);

%%
%%%% Plottings
figure(2)
clf
tiledlayout(1,3)

nexttile;
plot(t2,[y2_f; y2_b; y2_fb],t2,y2_exact,'--r')
legend('analytical [$E_f$]','analytical [$E_b$]','analytical [$E_{fb}$]','excact')
grid on
figtex(gca,2);
ylabel('$$y$$')
xlabel('time [s]')
ylim([-10 10])

nexttile;
semilogy(t2,abs(y2_exact - [y2_f; y2_b; y2_fb]))
grid on
figtex(gca,[]);
ylabel('$$|y_n - y(t)|$$')
xlabel('time [s]')

nexttile;
semilogy(t2,[y2_f; y2_b; y2_fb]./y2_0,t2,y2_exact./y2_0,'--r')
grid on
figtex(gca,[]);
ylabel('$$\frac{y_n}{y_0}$$')
xlabel('time [s]')

%%% The problem is itself a stiff problem (i.e, the model)
