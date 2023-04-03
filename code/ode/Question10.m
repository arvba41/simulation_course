clear all
clc

syms y(t) u(t) a 
gamma = 1.3;
y_0 = 1;

eqn = diff(y,t) == gamma*y/(cos(t) + 1.1)*(-sin(t));
cond = y(-pi) == y_0;
S = simplify(dsolve(eqn,cond));

y_value = @(time_vector) subs(S,{t},{time_vector}); % true solution
% y' = gamma*y/u(t)*u', y(-pi) = y_0
% u(t) = cos(t) + 1.1, u' = -sin(t)

%%


u = @(t) cos(t) - 1.1;
dot_u = @(t) -sin(t);
f = @(t,y) gamma*y./u(t).*dot_u(t);


%% Solving of equations

%%%%% Time step 1
h = 0.1;
ts = -pi;
tf = pi;

t1 = ts:h:tf; 
%%% Euler forward method
y1_f = euler_ft(@(y,t) f(t,y),y_0,t1,h); 
%%% Euler backward method
y1_b = euler_bt(@(y,t) f(t,y),y_0,t1,h); 
%%% Trapz method
y1_fb = euler_fbt(@(y,t) f(t,y),y_0,t1,h); 

y1_exact = y_value(t1); % true solution
f1 = f(t1,y1_exact)/y_0;

%%%%% Time step 1
h = 0.01;
ts = -pi;
tf = pi;

t2 = ts:h:tf; 
%%% Euler forward method
y2_f = euler_ft(@(y,t) f(t,y),y_0,t2,h); 
%%% Euler backward method
y2_b = euler_bt(@(y,t) f(t,y),y_0,t2,h); 
%%% Trapz method
y2_fb = euler_fbt(@(y,t) f(t,y),y_0,t2,h); 

y2_exact = y_value(t2); % true solution
f2 = f(t2,y2_exact)/y_0;

%% difference operator
y1 = [y1_f; y1_b; y1_fb]; % output vector with time step 1
y2 = [y2_f; y2_b; y2_fb]; % output vector with time step 1

delta_y1 = zeros(size(y1)); % difference operator time step 1
delta_y2 = zeros(size(y2)); % difference operator time step 1

delta_y1(:,2:end) = diff(y1,1,2); 
delta_y2(:,2:end) = diff(y2,1,2); 

%% Plots
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
title('step size $h_n = 0.1$')
ylabel('$$y(t)$$')
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'})

% Time step 2
nexttile;
plot(t2,[y2_f; y2_b; y2_fb],t2,y2_exact,'--r')
legend('analytical [$E_f$]','analytical [$E_b$]','analytical [$E_{fb}$]','excact')
figtex(gca,1);
grid on
title('step size $h_n = 0.01$')
ylabel('$$y(t)$$')
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'})

%---- Error
% Time step 1
nexttile;
semilogy(t1,abs(y1_exact - [y1_f; y1_b; y1_fb]))
figtex(gca,[]);
grid on
ylabel('$$|y_n - y(t_n)|$$')
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'})

% Time step 2
nexttile;
semilogy(t2,abs(y2_exact - [y2_f; y2_b; y2_fb]))
figtex(gca,[]);
grid on
ylabel('$$|y_n - y(t_n)|$$')
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'})

%---- difference operator
% Time step 1
nexttile;
plot(t1,delta_y1/y_0,t1,f1,'--r')
figtex(gca,[]);
grid on
ylabel('$$\frac{\Delta y_n}{\Delta t} $$')
xlabel('$\omega$ t [s]')
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'})

% Time step 1
nexttile;
plot(t2,delta_y2/y_0,t2,f2,'--r')
figtex(gca,[]);
grid on
ylabel('$$\frac{\Delta y_n}{\Delta t} $$')
xlabel('$\omega$ t [s]')
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'})