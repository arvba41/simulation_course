clear all
clc

% variables
% epsilon = [4e-2, 4e-3];% // \epsilon
ts = 0; tf = 15;
q = 8e-4;
f = 2/3;
e = 4e-2; % buffer variable

% problem
% [ \epsilon x'   = [ x(1-x) + f \frac{q-x}{q+x} z
%       z'     ]                x - z              ]
% 

fy = @(t,y,ii) [1/e*(y(1)*(1 - y(1)) + f*((q-y(1))/(q+y(1)))*y(2));...
            y(1) - y(2)];


%% solutions
y0 = [0.4231; 0.4231];

%     e = 4e-2;
fprintf ('epsilon = %0.2f non-stiff ode23 solver\n',e)
    tic; [t_nonstiff_1,y_nonstiff_1] = ode23(@(t,y) fy(t,y),[ts tf],y0); toc;% non-stiff solver
fprintf ('epsilon = %0.2f stiff ode15s solver\n',e)
    tic; [t_stiff_1,y_stiff_1] = ode15s(@(t,y) fy(t,y),[ts tf],y0); toc;% stiff solver
%     e = 4e-3;
e = 4e-3; % buffer variable
fy = @(t,y,ii) [1/e*(y(1)*(1 - y(1)) + f*((q-y(1))/(q+y(1)))*y(2));...
            y(1) - y(2)];
        
fprintf ('epsilon = %0.3f non-stiff ode23 solver\n',e)
    tic; [t_nonstiff_2,y_nonstiff_2] = ode23(@(t,y) fy(t,y),[ts tf],y0); toc;% non-stiff solver
fprintf ('epsilon = %0.3f stiff ode15s solver\n',e)
    tic; [t_stiff_2,y_stiff_2] = ode15s(@(t,y) fy(t,y),[ts tf],y0); toc;% stiff solver

%% plottings
figure(1)
tiledlayout(4,4,'TileSpacing','tight','Padding','compact')

nexttile([2 2]);
plot(y_nonstiff_1(:,1),y_nonstiff_1(:,2),y_stiff_1(:,1),y_stiff_1(:,2),'--')
grid on
ylabel('$z$')
xlabel('$x$')
legend('ode23','ode15s')
ax21 = figtex(gca,1);

nexttile([2 2]);
plot(y_nonstiff_2(:,1),y_nonstiff_2(:,2),y_stiff_2(:,1),y_stiff_2(:,2),'--')
grid on
xlabel('$x$')
figtex(gca);
ax22 = figtex(gca);

linkaxes([ax21 ax22],'xy')

nexttile([1 2]);
semilogy(t_nonstiff_1(2:end),diff(t_nonstiff_1),t_stiff_1(2:end),diff(t_stiff_1),'--')
grid on
ylabel('$h$')
xlabel('$t$')
title('$\epsilon = 4\cdot10^{-2}$')
ax3 = figtex(gca);

nexttile([1 2]);
semilogy(t_nonstiff_2(2:end),diff(t_nonstiff_2),t_stiff_2(2:end),diff(t_stiff_2),'--')
grid on
xlabel('$t$')
title('$\epsilon = 4\cdot10^{-3}$')
ax4 = figtex(gca);

linkaxes([ax3 ax4],'x')

nexttile;
plot(t_nonstiff_1,y_nonstiff_1(:,1),t_stiff_1,y_stiff_1(:,1),'--')
grid on
ylabel('$x(t)$')
xlabel('$t$')
title('$\epsilon = 4\cdot10^{-2}$')
ax5 = figtex(gca);

nexttile;
plot(t_nonstiff_2,y_nonstiff_2(:,1),t_stiff_2,y_stiff_2(:,1),'--')
grid on
xlabel('$t$')
title('$\epsilon = 4\cdot10^{-3}$')
ax6 = figtex(gca);

nexttile;
plot(t_nonstiff_1,y_nonstiff_1(:,2),t_stiff_1,y_stiff_1(:,2),'--')
grid on
ylabel('$z(t)$')
xlabel('$t$')
title('$\epsilon = 4\cdot10^{-2}$')
ax7 = figtex(gca);

nexttile;
plot(t_nonstiff_2,y_nonstiff_2(:,2),t_stiff_2,y_stiff_2(:,2),'--')
grid on
xlabel('$t$')
title('$\epsilon = 4\cdot10^{-3}$')
ax8 = figtex(gca);
linkaxes([ax5 ax6],'y')
linkaxes([ax7 ax8],'y')
linkaxes([ax3 ax4 ax5 ax6 ax7 ax8],'x')
