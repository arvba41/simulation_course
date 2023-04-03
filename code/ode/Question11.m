%% implement different methods for solving the equations
clear all
clc

% function definition: y'(t) = -10*y(t), y(0) = 1

ts = 0; % start time
tf = 10; % finish time

%% Solving of equations

%%%%% Time step 1
h = 1; % time step
t1 = ts:h:tf; % time vector 
%%% Euler forward method
y1_f = euler_f(@(y) 3*y,1,t1,h); 
%%% Euler backward method
y1_b = euler_b(@(y) 3*y,1,t1,h); 
%%% Trapz method
y1_fb = euler_fb(@(y) 3*y,1,t1,h); 

y1_exact = exp(3*t1); % true solution


% %%%%% Time step 2
% h = 0.1; % time step
% t2 = ts:h:tf; % time vector 
% %%% Euler forward method
% y2_f = euler_f(@(y) 3*y,1,t2,h); 
% %%% Euler backward method
% y2_b = euler_b(@(y) 3*y,1,t2,h); 
% %%% Trapz method
% y2_fb = euler_fb(@(y) 3*y,1,t2,h); 
% 
% y2_exact = exp(3*t2); % true solution
% 
% 
% Time step 1
h = 0.01; % time step
t3 = ts:h:tf; % time vector 
%%% Euler forward method
y3_f = euler_f(@(y) 3*y,1,t3,h); 
%%% Euler backward method
y3_b = euler_b(@(y) 3*y,1,t3,h); 
%%% Trapz method
y3_fb = euler_fb(@(y) 3*y,1,t3,h); 

y3_exact = exp(3*t3); % true solution
% 
% % Time step 1
% h = 0.001; % time step
% t4 = ts:h:tf; % time vector 
% %%% Euler forward method
% y4_f = euler_f(@(y) 3*y,1,t4,h); 
% %%% Euler backward method
% y4_b = euler_b(@(y) 3*y,1,t4,h); 
% %%% Trapz method
% y4_fb = euler_fb(@(y) 3*y,1,t4,h); 
% 
% y4_exact = exp(3*t4); % true solution

%% difference operator
y1 = [y1_f; y1_b; y1_fb]; % output vector with time step 1
% y2 = [y2_f; y2_b; y2_fb]; % output vector with time step 2
y3 = [y3_f; y3_b; y3_fb]; % output vector with time step 3
% y4 = [y4_f; y4_b; y4_fb]; % output vector with time step 4

delta_y1 = zeros(size(y1)); % difference operator time step 1
% delta_y2 = zeros(size(y2)); % difference operator time step 2
delta_y3 = zeros(size(y3)); % difference operator time step 3
% delta_y4 = zeros(size(y4)); % difference operator time step 4

delta_y1(:,2:end) = diff(y1,1,2); 
% delta_y2(:,2:end) = diff(y2,1,2); 
delta_y3(:,2:end) = diff(y3,1,2); 
% delta_y4(:,2:end) = diff(y4,1,2); 

%% Plotting
figure(1)
clf
tiledlayout(3,2)

%---- Solutions 
% Time step 1
nexttile;
plot(t1,[y1_f; y1_b; y1_fb],t1,y1_exact,'--r')
legend('analytical [$E_f$]','analytical [$E_b$]','analytical [$E_{fb}$]','excact')
figtex(gca,1);
grid on
title('step size $h_n = 1$')
ylabel('$$y(t)$$')

% % Time step 2
% nexttile;
% plot(t2,[y2_f; y2_b; y2_fb],t2,y2_exact,'--r')
% legend('analytical [$E_f$]','analytical [$E_b$]','analytical [$E_{fb}$]','excact')
% figtex(gca,1);
% grid on
% title('step size $h_n = 0.1$')
% 
% Time step 3
nexttile;
plot(t3,[y3_f; y3_b; y3_fb],t3,y3_exact,'--r')
legend('analytical [$E_f$]','analytical [$E_b$]','analytical [$E_{fb}$]','excact')
figtex(gca,1);
grid on
title('step size $h_n = 0.01$')
% 
% % Time step 4
% nexttile;
% plot(t4,[y4_f; y4_b; y4_fb],t4,y4_exact,'--r')
% legend('analytical [$E_f$]','analytical [$E_b$]','analytical [$E_{fb}$]','excact')
% figtex(gca,1);
% grid on
% title('step size $h_n = 0.001$')
% %---- Error

% Time step 1
nexttile;
semilogy(t1,abs(y1_exact - [y1_f; y1_b; y1_fb]))
figtex(gca,[]);
grid on
ylabel('$$|y_n - y(t_n)|$$')
xlabel('time [s]')

% % Time step 2
% nexttile;
% semilogy(t2,abs(y2_exact - [y2_f; y2_b; y2_fb]))
% figtex(gca,[]);
% grid on
% 
% Time step 3
nexttile;
semilogy(t3,abs(y3_exact - [y3_f; y3_b; y3_fb]))
figtex(gca,[]);
grid on
% 
% % Time step 4
% nexttile;
% semilogy(t4,abs(y4_exact - [y4_f; y4_b; y4_fb]))
% figtex(gca,[]);
% grid on

%---- difference operator
% Time step 1
nexttile;
plot(t1,delta_y1/0.25,t1,3*y1_exact,'--r')
figtex(gca,[]);
grid on
ylabel('$$\frac{\Delta y_n}{\Delta t} $$')
xlabel('time [s]')

% Time step 2
% nexttile;
% plot(t2,delta_y2/0.1,t2,3*y2_exact,'--r')
% figtex(gca,[]);
% grid on
% xlabel('time [s]')
% 
% Time step 3
nexttile;
plot(t3,delta_y3/0.01,t3,3*y3_exact,'--r')
figtex(gca,[]);
grid on
xlabel('time [s]')
% 
% % Time step 4
% nexttile;
% plot(t4,delta_y4/0.001,t4,3*y4_exact,'--r')
% figtex(gca,[]);
% grid on
% xlabel('time [s]')