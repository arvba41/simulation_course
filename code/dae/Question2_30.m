clear all; clc; 
% the DAE system for the pendulum
% parameters
m = 2.6; % mass of the pendulam
g = 9.8; 
l = 1; % length of the pendulam

% The Model
% y(1) = x
% y(2) = u
% y(3) = y'
% y(4) = y''
% y(5) = y
% y(6) = v
% y(7) = u'
% y(8) = v'
% y(9) = lambda
pendulam_idx1 = @(t,y) [ ...
    y(2);...
    y(7);...
    y(6) - y(3);...
    y(8) - y(4);...
    y(1)*y(9) - m*y(7);...
    y(5)*y(9) - m*g - m*y(8);...
    y(1)^2 + y(5)^2 - l^2;...
    - y(5)*y(3);...
    y(2)^2 + y(6)^2 + y(1)*y(7) + y(5)*y(4)]; 
% simulation 
tsim_l = [0 300]; % long simulation time

M_idx1 = zeros(9,9); % mass matrix
M_idx1(1,1) = 1; M_idx1(2,2) = 1; 
options_idx1 = odeset('Mass',M_idx1); 

y0 = -l; y0_idx0 = [0;1;0;0;y0;0;0;0;0]; % initial conditions

[t_idx1_l,y_idx1_l] = ode15s(@(t,y) pendulam_idx1(t,y),tsim_l,...
    y0_idx0,options_idx1);

%% Plots
figure(1)
clf

subplot(121)
plot(y_idx1_l(:,1),y_idx1_l(:,5)) %,y_idx0_s(:,1),y_idx0_s(:,2),'--')
grid on; legend('Index 1 -- Ode15s') %,'Index 0 -- Ode45')
ax(1) = figtex(gca,1);
xlabel('$x$'); ylabel('$y$'); title('Pendulum $t$ = 300s');
ylim([-l*1.1 0])

subplot(122)
semilogy(t_idx1_l(2:end),diff(t_idx1_l)) %,t_idx0_s(2:end),diff(t_idx0_s),'--')
grid on; ax(3) = figtex(gca);
xlabel('$t$ [s]'); ylabel('$h$');

figsize(1,0.35,16,24);
saveas(gcf,'Figures/Ugf2_30','epsc');

%% Modelica plots
addpath('Modelica results\')
data1_5s = csvread('exportedVariables.csv',1,0); 
data2_5min = csvread('exportedVariables5min.csv',1,0); 

figure(2); clf

subplot(221)
plot(data1_5s(:,2),data1_5s(:,3)) %,y_idx0_s(:,1),y_idx0_s(:,2),'--')
grid on; % legend('Index 1 -- Ode15s') %,'Index 0 -- Ode45')
ax(1) = figtex(gca);
xlabel('$x$'); ylabel('$y$'); title('Pendulum $t$ = 5s');
% ylim([-l*1.1 0])

subplot(222)
plot(data2_5min(:,2),data2_5min(:,3)) %,y_idx0_l(:,1),y_idx0_l(:,2),'--')
grid on; % legend('Index 1 -- Ode15s') %,'Index 0 -- Ode45')
ax(2) = figtex(gca);
xlabel('$x$'); ylabel('$y$'); title('Pendulum $t$ = 5min');
% ylim([-l*2 0])

subplot(223)
semilogy(data1_5s(2:end,1),diff(data1_5s(:,1))) %,t_idx0_s(2:end),diff(t_idx0_s),'--')
grid on; ax(3) = figtex(gca);
xlabel('$t$ [s]'); ylabel('$h$');

subplot(224)
semilogy(data2_5min(2:end,1),diff(data2_5min(:,1))) %,t_idx0_l(2:end),diff(t_idx0_l),'--')
grid on; ax(4) = figtex(gca);
xlabel('$t$ [s]'); ylabel('$h$');

% linkaxes(ax(1:2),'xy')
% linkaxes(ax(3:4),'y')

figsize(1,0.4);
saveas(gcf,'Figures/Ugf2_29','epsc');