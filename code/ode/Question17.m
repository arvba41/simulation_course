clear all
clc

% variables
sigma = 10;
r = 28;
b = 8/3; % buffer variable

y0 = [0; 1; 0]; % initial conditions

%% function definition
fy = @(t,y) [sigma*(y(2) - y(1));...
            r*y(1) - y(2) - y(1)*y(3);...
            y(1)*y(2) - b*y(3)];

%% part A
ts = 0; tf = 100;
options = odeset('AbsTol',1e-6,'RelTol',1e-6);

fprintf('Lorenz equations with tolerance %e. ',1e-6);
tic; [t6,y6] = ode45(@(t,y) fy(t,y),[ts tf],y0,options); toc;

%%% plots
figure(1)
clf
tiledlayout(3,1,'TileSpacing','tight','Padding','compact')

nexttile([2 1]);
plot(y6(:,1),y6(:,3))
grid on
ylabel('$y_3$'); xlabel('$y_1$');
ax1=figtex(gca);

nexttile;
plot(t6,y6(:,2))
grid on
ylabel('$y_2$'); xlabel('$t$');
ax2=figtex(gca);

%% Part B
figure(2)
plot3(y6(:,1),y6(:,2),y6(:,3))
grid on
xlabel('$y_1$'); ylabel('$y_2$'); zlabel('$y_3$');
ax3=figtex(gca);

%% part C
options = odeset('AbsTol',1e-5,'RelTol',1e-5);
fprintf('Lorenz equations with tolerance %e. ',1e-5);
tic; [t5,y5] = ode45(@(t,y) fy(t,y),[ts tf],y0,options); toc;

options = odeset('AbsTol',1e-7,'RelTol',1e-7);
fprintf('Lorenz equations with tolerance %e. ',1e-7);
tic; [t7,y7] = ode45(@(t,y) fy(t,y),[ts tf],y0,options); toc;

options = odeset('AbsTol',1e-8,'RelTol',1e-8);
fprintf('Lorenz equations with tolerance %e. ',1e-8);
tic; [t8,y8] = ode45(@(t,y) fy(t,y),[ts tf],y0,options); toc;

%%% plots
figure(3)
clf
tiledlayout(5,5,'TileSpacing','tight','Padding','compact')

nexttile([2 2]);
plot(y5(:,1),y5(:,3))
grid on
ylabel('$y_3$'); xlabel('$y_1$');
title('1e-5');
ax4=figtex(gca);

nexttile([2 2]);
plot(y6(:,1),y6(:,3))
grid on
ylabel('$y_3$'); xlabel('$y_1$');
title('1e-6');
ax5=figtex(gca);

nexttile([2 2]);
plot(y7(:,1),y7(:,3))
grid on
ylabel('$y_3$'); xlabel('$y_1$');
title('1e-7');
ax6=figtex(gca);

nexttile([2 2]);
plot(y8(:,1),y8(:,3))
grid on
ylabel('$y_3$'); xlabel('$y_1$');
title('1e-8');
ax7=figtex(gca);

nexttile;
plot(t5(2:end),diff(t5))
grid on
title('1e-5');
ax8a=figtex(gca);

nexttile;
plot(t6(2:end),diff(t6))
grid on
title('1e-6');
ax9a=figtex(gca);

nexttile;
plot(t7(2:end),diff(t7))
grid on
title('1e-7');
ax10a=figtex(gca);

nexttile;
plot(t8(2:end),diff(t8))
grid on
ylabel('$h$'); title('1e-8');
ax11a=figtex(gca);

nexttile;
plot(t5,y5(:,2))
grid on
ylabel('$y_2$'); xlabel('$t$');
title('1e-5');
ax8=figtex(gca);

nexttile;
plot(t6,y6(:,2))
grid on
xlabel('$t$');
title('1e-6');
ax9=figtex(gca);

nexttile;
plot(t7,y7(:,2))
grid on
xlabel('$t$'); title('1e-7');
ax10=figtex(gca);

nexttile;
plot(t8,y8(:,2))
grid on
xlabel('$t$'); title('1e-8');
ax11=figtex(gca);

nexttile;
plot3(y5(:,1),y5(:,2),y5(:,3))
hold on
plot3(y6(:,1),y6(:,2),y6(:,3))
plot3(y7(:,1),y7(:,2),y7(:,3))
plot3(y8(:,1),y8(:,2),y8(:,3))
grid on
xlabel('$y_1$'); ylabel('$y_2$'); zlabel('$y_3$');
ax12=figtex(gca);

linkaxes([ax4 ax5 ax6 ax7],'xy')
linkaxes([ax8 ax9 ax10 ax11],'xy')
linkaxes([ax8a ax9a ax10a ax11a],'xy')

Y100 = [y5(100,:); y6(100,:); y7(100,:); y8(100,:)]
