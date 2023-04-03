%% 2.14(b)
clc; clear all;

% constants
m = 2.6; % mass of the pendulam
g = 9.8; 
l = 1; % length of the pendulam

% function definition
% in pendulam_idx1 y(1) = x, y(2) = y, y(3) = u, y(4) = v, y(5) = lambda
% lambda from algebraic constraint
pendulam_idx1 = @(t,y) [y(3); y(4); y(5)*y(1); y(5)*y(2) - m*g;
                        y(3)^2 + y(4)^2 + y(5)*(y(1)^2 + y(2)^2)/m - g*y(2)];           

M_idx1 = [1 0 0 0 0;
          0 1 0 0 0;
          0 0 m 0 0;
          0 0 0 m 0;
          0 0 0 0 0];
options_idx1 = odeset('Mass',M_idx1); % Mass matrix for idx 1

M_idx0 = [1 0 0 0 0;
          0 1 0 0 0;
          0 0 m 0 0;
          0 0 0 m 0;
          0 0 0 0 1];
options_idx0 = odeset('Mass',M_idx0); % Mass matrix for idx 0


% simulation
tsim_s = [0 5]; % short simulation time
tsim_l = [0 300]; % long simulation time
x0 = -0.6;% initital cordinate
    if (abs(x0) > l)
        error("Choose initital value less than 'l'");
    end
y0_idx1 = [x0;sign(x0)*sqrt(l^2 - x0^2);0;0;0]; 
    % initital condition based on constraint y(1)^2 + y(2)^2 = l^2
y0_idx0 = [x0;sign(x0)*sqrt(l^2 - x0^2);0;0;0]; 

[t_idx1_s,y_idx1_s] = ode15s(@(t,y) pendulam_idx1(t,y),tsim_s,y0_idx1,options_idx1);
[t_idx1_l,y_idx1_l] = ode15s(@(t,y) pendulam_idx1(t,y),tsim_l,y0_idx1,options_idx1);
% [t_idx0_s,y_idx0_s] = ode45(@(t,y) pendulam_idx0(t,y),tsim_s,y0_idx0,options_idx0);
% [t_idx0_l,y_idx0_l] = ode45(@(t,y) pendulam_idx0(t,y),tsim_l,y0_idx0,options_idx0);

%% Plots
figure(1)
clf

subplot(221)
plot(y_idx1_s(:,1),y_idx1_s(:,2)) %,y_idx0_s(:,1),y_idx0_s(:,2),'--')
grid on; legend('Index 1 -- Ode15s') %,'Index 0 -- Ode45')
ax(1) = figtex(gca,1);
xlabel('$x$'); ylabel('$y$'); title('Pendulum $t$ = 5s');
ylim([-l*1.1 0])

subplot(222)
plot(y_idx1_l(:,1),y_idx1_l(:,2)) %,y_idx0_l(:,1),y_idx0_l(:,2),'--')
grid on; legend('Index 1 -- Ode15s') %,'Index 0 -- Ode45')
ax(2) = figtex(gca,1);
xlabel('$x$'); ylabel('$y$'); title('Pendulum $t$ = 5min');
ylim([-l*2 0])

subplot(223)
semilogy(t_idx1_s(2:end),diff(t_idx1_s)) %,t_idx0_s(2:end),diff(t_idx0_s),'--')
grid on; ax(3) = figtex(gca);
xlabel('$t$ [s]'); ylabel('$h$');

subplot(224)
semilogy(t_idx1_l(2:end),diff(t_idx1_l)) %,t_idx0_l(2:end),diff(t_idx0_l),'--')
grid on; ax(4) = figtex(gca);
xlabel('$t$ [s]'); ylabel('$h$');

linkaxes(ax(1:2),'xy')
linkaxes(ax(3:4),'y')

figsize(1,0.25);
saveas(gcf,'Figures/Ugf2_14a','epsc');

%% 2.14(c)
% Baumgarte stabilization

% extra constants 
alpha = 1; beta = 1; gamma = 1;

% in pendulam_idx1 y(1) = x, y(2) = y, y(3) = u, y(4) = v, y(5) = lambda
% lambda from algebraic constraint
pendulam_BS_idx1 = @(t,y) [y(3); y(4); y(5)*y(1); y(5)*y(2) - m*g; ...
                        (y(3)^2 + y(4)^2 + y(5)*(y(2)^2 + y(1)^2)/m - g*y(2)) + ...
                        alpha*(y(3)*y(1) + y(4)*y(2)) + beta*(y(1)^2 + y(2)^2 - l^2)];    
% in pendulam_idx0 y(1) = u, y(2) = v, y(3) = x, y(4) = y, y(5) = lambda
pendulam_BS_idx0 = @(t,y) [y(3); y(4); y(5)*y(1); y(5)*y(2) - m*g;
                        -(y(1)^2 + y(2)^2 - l^2 + alpha*(y(3)*y(1) + y(4)*y(2)) + ...
                        beta*(y(3)^2 + y(4)^2 + y(5)*(y(2)^2 + y(1)^2)/m - g*y(2)) + 3*g*y(4)*m*gamma)/l^2];
                    
[t_idx1_BS_s,y_idx1_BS_s] = ode15s(@(t,y) pendulam_BS_idx1(t,y),tsim_s,y0_idx1,options_idx1);
[t_idx1_BS_l,y_idx1_BS_l] = ode15s(@(t,y) pendulam_BS_idx1(t,y),tsim_l,y0_idx1,options_idx1);
[t_idx0_BS_s,y_idx0_BS_s] = ode45(@(t,y) pendulam_BS_idx0(t,y),tsim_s,y0_idx0,options_idx0);
[t_idx0_BS_l,y_idx0_BS_l] = ode45(@(t,y) pendulam_BS_idx0(t,y),tsim_l,y0_idx0,options_idx0);

%% Plots
figure(2)
clf

subplot(221)
plot(y_idx1_BS_s(:,1),y_idx1_BS_s(:,2)) %,y_idx0_s(:,1),y_idx0_s(:,2),'--')
grid on; legend('Index 1 -- Ode15s') %,'Index 0 -- Ode45')
ax(1) = figtex(gca,1);
xlabel('$x$'); ylabel('$y$'); title('Pendulum $t$ = 5s');
ylim([-l*1.1 0])

subplot(222)
plot(y_idx1_BS_l(:,1),y_idx1_BS_l(:,2)) %,y_idx0_l(:,1),y_idx0_l(:,2),'--')
grid on; legend('Index 1 -- Ode15s') %,'Index 0 -- Ode45')
ax(2) = figtex(gca,1);
xlabel('$x$'); ylabel('$y$'); title('Pendulum $t$ = 5min');
ylim([-l*1.1 0])

subplot(223)
semilogy(t_idx1_BS_s(2:end),diff(t_idx1_BS_s)) %,t_idx0_s(2:end),diff(t_idx0_s),'--')
grid on; ax(3) = figtex(gca);
xlabel('$t$ [s]'); ylabel('$h$');

subplot(224)
semilogy(t_idx1_BS_l(2:end),diff(t_idx1_BS_l)) %,t_idx0_l(2:end),diff(t_idx0_l),'--')
grid on; ax(4) = figtex(gca);
xlabel('$t$ [s]'); ylabel('$h$');

linkaxes(ax(1:2),'xy')
linkaxes(ax(3:4),'y')

figsize(1,0.25);
saveas(gcf,'Figures/Ugf2_14c','epsc');

%% 2.14(d)

epsilon = [0 1e-3 1e-6]; 

% function definition
% in pendulam_idx1 y(1) = x, y(2) = y, y(3) = u, y(4) = v, y(5) = lambda
% lambda from algebraic constraint
pendulam_idx1_zz = @(t,y) [y(3); y(4); y(5)*y(1); y(5)*y(2) - m*g;
                        y(3)^2 + y(4)^2 + y(5)*(y(1)^2 + y(2)^2)/m - g*y(2)]; 
        % Direct index reduction technique
        
M_idx1_zz = [1 0 0 0 0;
          0 1 0 0 0;
          0 0 m 0 0;
          0 0 0 m 0;
          0 0 0 0 0]; % Mass matrix

M_idx1_zz(end) = epsilon(1);
options_idx1_zz = odeset('Mass',M_idx1_zz); 
[t_idx1_l_zz_1,y_idx1_l_zz_1] = ode15s(@(t,y) pendulam_idx1_zz(t,y),tsim_s,y0_idx1,options_idx1_zz);

M_idx1_zz(end) = epsilon(2);
options_idx1_zz = odeset('Mass',M_idx1_zz); 
[t_idx1_l_zz_2,y_idx1_l_zz_2] = ode45(@(t,y) pendulam_idx1_zz(t,y),tsim_s,y0_idx1,options_idx1_zz);

M_idx1_zz(end) = epsilon(3);
options_idx1_zz = odeset('Mass',M_idx1_zz); 
[t_idx1_l_zz_3,y_idx1_l_zz_3] = ode45(@(t,y) pendulam_idx1_zz(t,y),tsim_s,y0_idx1,options_idx1_zz);

%% Plots
figure(3)
clf

subplot(231)
plot(y_idx1_l_zz_1(:,1),y_idx1_l_zz_1(:,2)) %,y_idx0_s(:,1),y_idx0_s(:,2),'--')
grid on; 
ax(1) = figtex(gca);
xlabel('$x$'); ylabel('$y$'); title('$\epsilon$ = 0');
ylim([-l*1.1 0])

subplot(232)
plot(y_idx1_l_zz_2(:,1),y_idx1_l_zz_2(:,2)) %,y_idx0_l(:,1),y_idx0_l(:,2),'--')
grid on;
ax(2) = figtex(gca);
xlabel('$x$'); ylabel('$y$'); title('$\epsilon$ = 1m');

subplot(233)
plot(y_idx1_l_zz_3(:,1),y_idx1_l_zz_3(:,2)) %,y_idx0_l(:,1),y_idx0_l(:,2),'--')
grid on;
ax(3) = figtex(gca);
xlabel('$x$'); ylabel('$y$'); title('$\epsilon$ = 1$\mu$');

subplot(234)
semilogy(t_idx1_l_zz_1(2:end),diff(t_idx1_l_zz_1)) %,t_idx0_l(2:end),diff(t_idx0_l),'--')
grid on; ax(4) = figtex(gca);
xlabel('$t$ [s]'); ylabel('$h$');

subplot(235)
semilogy(t_idx1_l_zz_2(2:end),diff(t_idx1_l_zz_2)) %,t_idx0_l(2:end),diff(t_idx0_l),'--')
grid on; ax(5) = figtex(gca);
xlabel('$t$ [s]'); ylabel('$h$');

subplot(236)
semilogy(t_idx1_l_zz_3(2:end),diff(t_idx1_l_zz_3)) %,t_idx0_l(2:end),diff(t_idx0_l),'--')
grid on; ax(6) = figtex(gca);
xlabel('$t$ [s]'); ylabel('$h$');

% linkaxes(ax(1:3),'xy')
linkaxes(ax(4:6),'x')
figsize(1,0.3);
saveas(gcf,'Figures/Ugf2_14d','epsc');

%% Sensitivity on alpha and beta
% alpha_range = -1:1; beta_range = -1:1;
% % stab_func = @(t,s,a,b) s^2 + a*s + b; % stability function (pole zero plot function)
% 
% for ii = 1:length(alpha_range)
%     for jj = 1:length(beta_range)
%         alpha = alpha_range(ii);
%         beta = beta_range(jj);
%         % simulation 
%         [t_idx1_BS_l_range(:,ii,jj),y_idx1_BS_l_range(:,:,ii,jj)] = ode15s(@(t,y) pendulam_BS_idx1(t,y),tsim_l,y0_idx1,options_idx1);
%     end
% end
% 
% %% Plots
% figure(3)
% clf
% 
% subplot(231)
% plot(squeeze(y_idx1_BS_l_range(:,1,find(alpha_range==-1),find(beta_range==1))),squeeze(y_idx1_BS_l_range(:,2,find(alpha_range==-1),find(beta_range==1))))
% grid on; ax(1) = figtex(gca);
% xlabel('$x$'); ylabel('$y$'); title('Pendulum $t$ = 5min, $\alpha$ = -1, $\beta$ = 1');
% ylim([-l*1.1 0]); xlim([-1.1 1.1]);
% 
% subplot(232)
% plot(squeeze(y_idx1_BS_l_range(:,1,find(alpha_range==0),find(beta_range==0))),squeeze(y_idx1_BS_l_range(:,2,find(alpha_range==0),find(beta_range==0))))
% grid on; %,'Index 0 -- Ode45')
% ax(2) = figtex(gca);
% xlabel('$x$'); ylabel('$y$'); title('Pendulum $t$ = 5min, $\alpha$ = 0, $\beta$ = -1');
% 
% subplot(233)
% plot(squeeze(y_idx1_BS_l_range(:,1,find(alpha_range==1),find(beta_range==-1))),squeeze(y_idx1_BS_l_range(:,2,find(alpha_range==1),find(beta_range==-1))))
% grid on; %,'Index 0 -- Ode45')
% ax(3) = figtex(gca);
% xlabel('$x$'); ylabel('$y$'); title('Pendulum $t$ = 5min, $\alpha$ = 1, $\beta$ = -1');
% 
% subplot(234)
% semilogy(squeeze(t_idx1_BS_l_range(2:end,find(alpha_range==-1),find(beta_range==1))),diff(squeeze(t_idx1_BS_l_range(:,find(alpha_range==-1),find(beta_range==1)))))%,t_idx0_s(2:end),diff(t_idx0_s),'--')
% grid on; ax(4) = figtex(gca);
% xlabel('$t$ [s]'); ylabel('$h$');
% 
% subplot(235)
% semilogy(squeeze(t_idx1_BS_l_range(2:end,find(alpha_range==0),find(beta_range==0))),diff(squeeze(t_idx1_BS_l_range(:,find(alpha_range==0),find(beta_range==0)))))%,t_idx0_s(2:end),diff(t_idx0_s),'--')
% grid on; ax(5) = figtex(gca);
% xlabel('$t$ [s]'); ylabel('$h$');
% 
% subplot(236)
% semilogy(squeeze(t_idx1_BS_l_range(2:end,find(alpha_range==1),find(beta_range==-1))),diff(squeeze(t_idx1_BS_l_range(:,find(alpha_range==1),find(beta_range==-1)))))%,t_idx0_s(2:end),diff(t_idx0_s),'--')
% grid on; ax(6) = figtex(gca);
% xlabel('$t$ [s]'); ylabel('$h$');
% 
% linkaxes(ax(1:3),'xy')
% linkaxes(ax(3:6),'y')
