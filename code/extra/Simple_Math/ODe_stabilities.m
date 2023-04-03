clear all; clc; 

%% considering a simple unstable system 
T = 10; % stop time 

f = @(t,y) y-sin(t)-cos(t); % function 
ts = linspace(0,T,40); % Ns = 40 % number of steps
tss = linspace(0,T,1000); % Ns = 1000 % number of steps

ys = ef(f,1,ts,mean(diff(ts))); % solving ode
yss = ef(f,1,tss,mean(diff(tss))); % solving ode

tv = linspace(0,T,500); % exact solution time steps

tval = linspace(0,T,20); % exact solution time steps
yval = linspace(-1,6,20); % y value for the dir field

% changing the initial conditions but keeing the same step size
ys_in1 = ef(f,1,ts,mean(diff(ts))); % solving ode
ys_in2 = ef(f,0,ts,mean(diff(ts))); % solving ode

figure(1); clf;
subplot(3,2,1);
plot(tv,cos(tv),'k',ts,ys,'bo-',tss,yss,'r.-'); hold on; 
dirplotter(f,tval,yval); hold off; legend('$y(t)$','$y_n$','$y_n$ @low $h_n$')
ax(1) = figtex(gca,1); ylabel('$y(t)$'), xlabel('$t$'); grid on;
title('solution and field lines');

subplot(3,2,2);
plot(tv,cos(tv),'k',ts,ys_in1,'bo-',ts,ys_in2,'g.-'); hold on; 
dirplotter(f,tval,yval); hold off; legend('$y(t)$','$y_n\ y(0) = 1$ ','$y_n\ y(0) = 0$')
ax(1) = figtex(gca,1); ylabel('$y(t)$'), xlabel('$t$'); grid on;
title('solution and field lines');

subplot(3,2,3);
semilogy(ts,abs(cos(ts) - ys),'bo-',tss,abs(cos(tss) - yss),'r.-'); 
legend('$|y(t) - y_n|$','$|y(t) - y_n$ @low $h_n|$'); title('Global error')
ax(2) = figtex(gca,1); ylabel('$e(t)$'); xlabel('t'); grid on;

subplot(3,2,4);
semilogy(ts,abs(cos(ts) - ys_in1),'bo-',ts,abs(cos(ts) - ys_in2),'g.-'); 
legend('$|y(t) - y_n\ y(0) = 1|$','$|y(t) - y_n\ y(0) = 1$|'); title('Global error')
ax(2) = figtex(gca,1); ylabel('$e(t)$'); xlabel('t'); grid on;

subplot(3,2,6);
semilogy(ts, abs(ys_in1 - ys_in2)); grid on
title('error between the changes')
%% considering a simple stable system 
T = 10; % stop time 

f = @(t,y) -y-sin(t)-cos(t); % function 
ts = linspace(0,T,40); % Ns = 40 % number of steps
tss = linspace(0,T,1000); % Ns = 1000 % number of steps

ys = ef(f,1,ts,mean(diff(ts))); % solving ode
yss = ef(f,1,tss,mean(diff(tss))); % solving ode

tv = linspace(0,T,500); % exact solution time steps

tval = linspace(0,T,20); % exact solution time steps
yval = linspace(-1,6,20); % y value for the dir field

% changing the initial conditions but keeing the same step size
ys_in1 = ef(f,1,ts,mean(diff(ts))); % solving ode
ys_in2 = ef(f,0,ts,mean(diff(ts))); % solving ode

figure(2); clf;
subplot(3,2,1);
plot(tv,-sin(tv),'k',ts,ys,'bo-',tss,yss,'r.-'); hold on; 
dirplotter(f,tval,yval); hold off; legend('$y(t)$','$y_n$','$y_n$ @low $h_n$')
ax(1) = figtex(gca,1); ylabel('$y(t)$'), xlabel('$t$'); grid on;
title('solution and field lines');

subplot(3,2,2);
plot(tv,-sin(tv),'k',ts,ys_in1,'bo-',ts,ys_in2,'g.-'); hold on; 
dirplotter(f,tval,yval); hold off; legend('$y(t)$','$y_n\ y(0) = 1$ ','$y_n\ y(0) = 0$')
ax(1) = figtex(gca,1); ylabel('$y(t)$'), xlabel('$t$'); grid on;
title('solution and field lines');

subplot(3,2,3);
semilogy(ts,abs(-sin(ts) - ys),'bo-',tss,abs(-sin(tss) - yss),'r.-'); 
legend('$|y(t) - y_n|$','$|y(t) - y_n$ @low $h_n|$'); title('Global error')
ax(2) = figtex(gca,1); ylabel('$e(t)$'); xlabel('t'); grid on;

subplot(3,2,4);
semilogy(ts,abs(-sin(ts) - ys_in1),'bo-',ts,abs(-sin(ts) - ys_in2),'g.-'); 
legend('$|y(t) - y_n\ y(0) = 1|$','$|y(t) - y_n\ y(0) = 1$|'); title('Global error')
ax(2) = figtex(gca,1); ylabel('$e(t)$'); xlabel('t'); grid on;

subplot(3,2,6);
semilogy(ts, abs(ys_in1 - ys_in2)); grid on
title('error between the changes')


%%%%%%--------------------------------------------------------------------
%%% unstabe ode
% T = 5;
% f = @(t,y) y-sin(t)-cos(t);
% [ts,ys] = Euler(f,[0,T],1,20);    % use N=20 steps
% 
% dirfield(f,0:.2:T,-1.1:.25:6); hold on
% tv = linspace(0,T,500);
% plot(tv,cos(tv),'k',ts,ys,'bo-'); hold off