clear all; clc; close all;
%% Setting up the system

p1 = 0.04; p2 = 10e4; p3 = 3e7; % Constants 
y0 = [1; 2e-4; 3e-1]; % initial conditions

% function definition
robertsdae = @(t,y) [-p1*y(1) + p2*y(2).*y(3); 
                p1*y(1) - p2*y(2).*y(3) - p3*y(2).^2;
                p3*y(2).^2];

% simulation
tsim = [0 100]; % simulation time
% M = [1 0 0; 0 1 0; 0 0 1]; % Mass matrix
% options = odeset('Mass',M);
% play with the tols and try with ode15i

[t,y] = ode45(@(t,y) robertsdae(t,y),tsim,y0);

figure(1)
plot(t,sum(y,2)); grid on;
ylabel('$x_1 + x_2 + x_3$'); xlabel('t [s]');
ax = figtex(gca); figsize(1,0.25);
title('Solver: \texttt{ode45}')

saveas(gca,'Figures/Ugf18c','epsc');

%% the circle drawer
circledrawer = @(t,y) [-y(2);...
                        y(1)];

y0 = [1; 0]; % initial conditions                    
% simulation
tsim = [0 100]; % simulation time

[t,y] = ode45(@(t,y) circledrawer(t,y),tsim,y0);

figure(2)
plot(t,sum(y.^2,2)); grid on;
ylabel('$y_1^2 + y_2^2$'); xlabel('t [s]');
ax = figtex(gca); figsize(1,0.25);
title('Solver: \texttt{ode45}')

saveas(gca,'Figures/Ugf18d','epsc');

