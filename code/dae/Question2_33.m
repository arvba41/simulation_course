clc; clear all;

p1 = 0.04; p2 = 10e4; p3 = 3e7; % Constants 
y0 = [1; 0; 0]; % initial conditions

% function definition
robertsdae = @(t,y) [-p1*y(1) + p2*y(2).*y(3); 
                p1*y(1) - p2*y(2).*y(3) - p3*y(2).^2;
                y(1) + y(2) + y(3) - 1 ];

% simulation
tsim = [0 40e10]; % simulation time
M = [1 0 0; 0 1 0; 0 0 0]; % Mass matrix
options = odeset('Mass',M);
% play with the tols and try with ode15i

[t,y] = ode15s(@(t,y) robertsdae(t,y),tsim,y0,options);

% plot the simulation
figure(1)

subplot(2,1,1)
y(:,2) = 1e4*y(:,2); % scaling y_2
semilogx(t,y);
title('Robertson DAE problem with a Conservation Law, solved by ODE15S');
legend('$y_1$','$y_2$','$y_3$');
ax(1) = figtex(gca,1);
xlabel('t [s]'); ylabel('$y$');
grid on

subplot(2,1,2)
loglog(t(2:end),diff(t),'*-');
ax(2) = figtex(gca);
xlabel('t [s]'); ylabel('$h$');
grid on

figsize(0.9,0.3*0.9);
saveas(gcf,'Figures/Ugf233','epsc');
