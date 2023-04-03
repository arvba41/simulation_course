%% Question 2
clc
clear all
syms y(t) a
eqn = diff(y,t) == -a*(y-cos(t));
cond = y(0) == 0;
S = simplify(dsolve(eqn,cond));

% considering a = 10
y_value = @(time_vector) subs(S,{a,t},{10,time_vector}); % true solution

%% Question 3
% function definition: y'(t) = -a*(y-cos(t)), y(0) = 0

ts = 0; % start time
tf = 6; % finish time

a = 10; % cosntant
%% Solving of equations

tic
%%%%% Time step 1
h = 1; % time step
t1 = ts:h:tf; % time vector 
%%% Euler forward method
y1_f = euler_ft(@(y,t) -a*(y-cos(t)),0,t1,h); 
%%% Euler backward method
y1_b = euler_bt(@(y,t) -a*(y-cos(t)),0,t1,h); 
%%% Trapz method
y1_fb = euler_fbt(@(y,t) -a*(y-cos(t)),0,t1,h); 

y1_exact = y_value(t1); % true solution

%%%%% Time step 2
h = 0.1; % time step
t2 = ts:h:tf; % time vector 
%%% Euler forward method
y2_f = euler_ft(@(y,t) -a*(y-cos(t)),0,t2,h); 
%%% Euler backward method
y2_b = euler_bt(@(y,t) -a*(y-cos(t)),0,t2,h); 
%%% Trapz method
y2_fb = euler_fbt(@(y,t) -a*(y-cos(t)),0,t2,h); 

y2_exact = y_value(t2); % true solution

%%%%% Time step 3
h = 0.01; % time step
t3 = ts:h:tf; % time vector 
%%% Euler forward method
y3_f = euler_ft(@(y,t) -a*(y-cos(t)),0,t3,h); 
%%% Euler backward method
y3_b = euler_bt(@(y,t) -a*(y-cos(t)),0,t3,h); 
%%% Trapz method
y3_fb = euler_fbt(@(y,t) -a*(y-cos(t)),0,t3,h); 

y3_exact = y_value(t3); % true solution

%%%%% Time step 4
h = 0.001; % time step
t4 = ts:h:tf; % time vector 
%%% Euler forward method
y4_f = euler_ft(@(y,t) -a*(y-cos(t)),0,t4,h); 
%%% Euler backward method
y4_b = euler_bt(@(y,t) -a*(y-cos(t)),0,t4,h); 
%%% Trapz method
y4_fb = euler_fbt(@(y,t) -a*(y-cos(t)),0,t4,h); 

y4_exact = y_value(t4); % true solution
toc 
%% Plotting
figure(1)
clf
tiledlayout(2,4)

tic 
% Time step 1
nexttile;
plot(t1,[y1_f; y1_b; y1_fb],t1,y1_exact,'--r')
legend('analytical [$E_f$]','analytical [$E_b$]','analytical [$E_{fb}$]','excact')
figtex(gca,1);
grid on
title('step size $h_n = 1$')
ylabel('$$y(t)$$')

% Time step 2
nexttile;
plot(t2,[y2_f; y2_b; y2_fb],t2,y2_exact,'--r')
legend('analytical [$E_f$]','analytical [$E_b$]','analytical [$E_{fb}$]','excact')
figtex(gca,1);
grid on
title('step size $h_n = 0.1$')

% Time step 3
nexttile;
plot(t3,[y3_f; y3_b; y3_fb],t3,y3_exact,'--r')
legend('analytical [$E_f$]','analytical [$E_b$]','analytical [$E_{fb}$]','excact')
figtex(gca,1);
grid on
title('step size $h_n = 0.01$')

% Time step 4
nexttile;
plot(t4,[y4_f; y4_b; y4_fb],t4,y4_exact,'--r')
legend('analytical [$E_f$]','analytical [$E_b$]','analytical [$E_{fb}$]','excact')
figtex(gca,1);
grid on
title('step size $h_n = 0.001$')

% Time step 1
nexttile;
semilogy(t1,abs(y1_exact - [y1_f; y1_b; y1_fb]))
figtex(gca,[]);
grid on
ylabel('$$y_n - y(t_n)$$')
xlabel('time [s]')

% Time step 2
nexttile;
semilogy(t2,abs(y2_exact - [y2_f; y2_b; y2_fb]))
figtex(gca,[]);
grid on
xlabel('time [s]')

% Time step 3
nexttile;
semilogy(t3,abs(y3_exact - [y3_f; y3_b; y3_fb]))
figtex(gca,[]);
grid on
xlabel('time [s]')

% Time step 4
nexttile;
semilogy(t4,abs(y4_exact - [y4_f; y4_b; y4_fb]))
figtex(gca,[]);
grid on
xlabel('time [s]')
toc 

%% Second part of the problem 

a = 1000; % constant
%%%%% Time step 1
h = 0.01; % time step
ts = 0; % start time
tf = 2; % start time
t_vec = ts:h:tf; % time vector 

tic 
%%% Euler backward method
y_b = euler_bt(@(y,t) -a*(y-cos(t)),0,t_vec,h); 
%%% Trapz method
y_fb = euler_fbt(@(y,t) -a*(y-cos(t)),0,t_vec,h); 

y_exact = y_value(t_vec); % true solution

figure(2)
tiledlayout(1,2)

nexttile;
plot(t_vec,[y_b; y_fb],t_vec,y_exact,'--r')
legend('analytical [$E_b$]','analytical [$E_{fb}$]','excact')
grid on
figtex(gca,1);
ylabel('$$y(t)$$')
xlabel('time [s]')

nexttile;
semilogy(t_vec,abs(y_exact - [y_b; y_fb]))
grid on
figtex(gca,[]);
ylabel('$$y_n - y(t_n)$$')
xlabel('time [s]')

%%%% Comments 
% the solver seems to overcorrect at the beginning but the euler backward
% does not do the same 