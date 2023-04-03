clear all; clc; close all;
%% the function definitiion

% the parameters
p = [21.8893;... p1(0)
      2.14e9;... p2(0)
      32.318;... p3(0)
      21.893;... p4(0)
      1.07e9;... p5(0)
      7.65e-18;... p6(0)
      4.03e-11;... p7(0)
      5.32e-18;... p8(0)
      ];

% the parameters
uo1 = 1.5776;
u7 = 0.5*(p(7) + sqrt(p(7)^2 + 4*p(7)*uo1)); u8 = u7;
u0 = [uo1; 8.32; 0; 0; 0; 0.0131; u7; u8; 0; 0];


% the batch reactor problem
batch_reactor = @(t,u) [...
        -p(3)*u(2)*u(8)                                           ; ...u1'
        -p(1)*u(2)*u(6) + p(2)*u(10)     - p(3)*u(2)*u(8)         ;... u2'
         p(3)*u(2)*u(8) + p(4)*u(4)*u(6) - p(5)*u(9)              ;... u3'
        -p(4)*u(4)*u(6) + p(5)*u(9)                               ;... u4'
         p(1)*u(2)*u(6) - p(2)*u(10)                              ;... u5'
        -p(1)*u(2)*u(6) - p(4)*u(4)*u(6) + p(2)*u(10) + p(5)*u(9) ;... u6'
        -0.0131 + u(6) + u(8) + u(9) + u(10)                      ;... u7'    
        u(8)  - p(7)*u(1)/(p(7) + u(7))                           ;... u8
        u(9)  - p(8)*u(3)/(p(8) + u(7))                           ;... u9
        u(10) - p(6)*u(5)/(p(6) + u(7))                           ;... u10
    ];

% the mass matrix
Nstep = diag([ones([1,7]) 0 0 0]);

%% the simualtion

% the time step
ts = [0 10]; % step-time of 0.001 s

options = odeset('Mass',Nstep); % the mass matrix

% BDF implicit solver
[tsim,usim] = ode15s(@(t,y) batch_reactor(t,y),ts,u0,options); 

%% the plots

figure(1)
clf

subplot(211)
plot(tsim,[usim(:,1), usim(:,3), usim(:,4), usim(:,5)])
ylabel('$u$'); xlabel('$t$ [hours]'); grid on;
legend('$u_1$','$u_3$','$u_4$','$u_5$')
ax1(1) = figtex(gca,1); figsize(1,0.3);
title('Bath-Reactor proboelm state-variables: $u_6$,$u_8$')

subplot(212)
plot(tsim,[usim(:,6), usim(:,8)])
ylabel('$u$'); xlabel('$t$ [hours]'); grid on;
legend('$u_6$','$u_8$')
ax1(2) = figtex(gca,1); figsize(1,0.3);
title('Bath-Reactor proboelm state-variables: $u_6$,$u_8$')

saveas(gcf,'Figures/Ugf2_9c1','epsc');
linkaxes(ax1,'x')

%% the sensitivity analysis
% using symbolic toolbox to find out all the derivatives

% clear all% temp

p_syms = sym('p',size(p)); % parameters 
u_syms = sym('u',size(u0)); % inputs 

% xp_syms = sym('xp_%d_%d',[10 8]); % states

f_syms = [...
        -p_syms(3)*u_syms(2)*u_syms(8)                                                                   ; ...u1'
        -p_syms(1)*u_syms(2)*u_syms(6) + p_syms(2)*u_syms(10)     - p_syms(3)*u_syms(2)*u_syms(8)        ;... u2'
         p_syms(3)*u_syms(2)*u_syms(8) + p_syms(4)*u_syms(4)*u_syms(6) - p_syms(5)*u_syms(9)             ;... u3'
        -p_syms(4)*u_syms(4)*u_syms(6) + p_syms(5)*u_syms(9)                                             ;... u4'
         p_syms(1)*u_syms(2)*u_syms(6) - p_syms(2)*u_syms(10)                                            ;... u5'
        -p_syms(1)*u_syms(2)*u_syms(6) - p_syms(4)*u_syms(4)*u_syms(6) + p_syms(2)*u_syms(10) + p_syms(5)*u_syms(9) ;... u6'
        -0.0131 + u_syms(6) + u_syms(8) + u_syms(9) + u_syms(10)                                         ;... u7'    
        u_syms(8)  - p_syms(7)*u_syms(1)/(p_syms(7) + u_syms(7))                                         ;... u8
        u_syms(9)  - p_syms(8)*u_syms(3)/(p_syms(8) + u_syms(7))                                         ;... u9
        u_syms(10) - p_syms(6)*u_syms(5)/(p_syms(6) + u_syms(7))                                         ;... u10
    ];% f_syms(t) = u_syms'(t)

for ii = 1:length(p_syms)
    fp_syms(:,ii) = diff(f_syms,p_syms(ii)); % fp
end

for ii = 1:length(u_syms)
    fx_syms(:,ii) = diff(f_syms,u_syms(ii)); % fx
end

nx = length(u_syms);
np = length(p_syms);
%% simulation
% f = matlabFunction(f_syms,'Vars',[{u_syms(1)} {u_syms(2)} {u_syms(3)} {u_syms(4)} {u_syms(5)} {u_syms(6)} {u_syms(7)} {u_syms(8)} {u_syms(9)} {u_syms(10)} {params_syms}]);
% f = matlabFunction(f_syms,'Vars',[{u_syms} {params_syms}]);

% % time step
% Nstep = 500;
% tspan = linspace(0, 10, Nstep);
% xp = zeros(np, Nstep, 2 * nx);
% 
% % mass matrix
% M = diag([ones([1,7]) 0 0 0 ones([1,7]) 0 0 0]);
% options = odeset('Mass',M); % the mass matrix
% 
% for ii = 1:np
%     P = sym('P', [nx, 1]);
%     
%     % Add sensitivity functions
%     f_sens_sym = [f_syms; fx_syms(:, :) * P + fp_syms(:, ii)];
%     F = matlabFunction(f_sens_sym, 'vars', [{u_syms} {p_syms} {P}]);
%     f_func = @(t, y) F(y(1:nx), p, y(nx + 1:end));
%     
%     % initial conditions
%     xp0 = [u0;zeros(10,1)];
%     
%     % simulation 
%     [~, xp(ii, :, :)] = ode15s(@(t, y) f_func(t,y), tspan, xp0, options);
%     
% end
%% simulation 
clear tsens Psens Ysens 
Xp0 = zeros(length(u0),length(p));
tvec = 0:0.001:10; % time vector 
% f = waitbar(0,'Started Simulation');

% for jj = 1:length(p)
for ii = 1:length(p)
        % \dot x = f_x\,x_p + f_p
%         xp0 = zeros(size(fx_syms(:,jj)))';
        
%         fx = matlabFunction(fx_syms(ii,:),'Vars',[{u_syms} {params_syms}]);
%         fxn = @(t) fx(interp1(tsim,usim,t,'spline'),p);
%         
% %         fp = matlabFunction(fp_syms(ii,jj),'Vars',[{u_syms} {params_syms}]);
%         fpn = @(t) fp(interp1(tsim,usim,t,'spline'),p);
%         
%         xp = @(t,xp) fxn(t)*xp +  fpn(t);
%         Ysens = euler_f_2dvec(@(t,y) xp(t,y),tvec,Xp0);
%         [tsens,Ysens] = ode45(@(t,y) xp(ts,y),tsim,Xp0);
%         Psens(jj,ii).res = Ysens;
%         Psens(jj,ii).res = tsens;
%         waitbar((ii+jj*length(Tem_req))/(length(Tem_req) * length(wem_req)),f,'Simulating %3.2f',(ii+jj*length(Tem_req))/(length(Tem_req) * length(wem_req))*100);
%     end
% end
fp = matlabFunction(fp_syms(:,ii),'Vars',[{u_syms} {p_syms}]);
fx = matlabFunction(fx_syms,'Vars',[{u_syms} {p_syms}]);

batch_reactor_sens = @(t,x) [...
        -p(3)*x(2)*x(8)                                           ; ...u1'
        -p(1)*x(2)*x(6) + p(2)*x(10)     - p(3)*x(2)*x(8)         ;... u2'
         p(3)*x(2)*x(8) + p(4)*x(4)*x(6) - p(5)*x(9)              ;... u3'
        -p(4)*x(4)*x(6) + p(5)*x(9)                               ;... u4'
         p(1)*x(2)*x(6) - p(2)*x(10)                              ;... u5'
        -p(1)*x(2)*x(6) - p(4)*x(4)*x(6) + p(2)*x(10) + p(5)*x(9) ;... u6'
        -0.0131 + x(6) + x(8) + x(9) + x(10)                      ;... u7'    
        x(8)  - p(7)*x(1)/(p(7) + x(7))                           ;... u8
        x(9)  - p(8)*x(3)/(p(8) + x(7))                           ;... u9
        x(10) - p(6)*x(5)/(p(6) + x(7))                           ;... u10
        ... batch reactorproblem 
        fx(x(1:10),p)*x(11:20) + fp(x(1:10),p);... xp1'
    ];

    %% sensitivity simulation 

    xp0 = [u0;fp(u0,p)];

    M = diag([ones([1,7]) 0 0 0 ones([1,7]) 0 0 0]);
    options = odeset('Mass',M); % the mass matrix

    % BDF implicit solver
[xpsim(ii).t,xpsim(ii).y] = ode15s(@(t,y) batch_reactor_sens(t,y),ts,xp0,options); 
end
%% plotting sensitivity analysis 

for jj=1:length(p)
figure(1+jj); clf

subplot(211)
for ii = 11:15
    plot(xpsim(jj).t,xpsim(jj).y(:,ii)*p(jj)); hold on;
end
ytext = char("$p_" + num2str(jj) + "du/dp_" + num2str(jj) + "$");
ylabel(ytext); xlabel('$t$ [hours]'); grid on; box off;
legend('$u_1$','$u_2$','$u_3$','$u_4$','$u_5$')
% ax1(1) = figtex(gca,1); figsize(1,0.3);
title('Bath-Reactor proboelm state-variables: $u_1$,$u_2$,$u_3$,$u_4$,and $u_5$');

subplot(212)
for ii = 15:20
    plot(xpsim(jj).t,xpsim(jj).y(:,ii)*p(jj)); hold on;
end
ytext = char("$p_{" + num2str(jj) + "} du/dp_{" + num2str(jj) + "}$");
ylabel(ytext); xlabel('$t$ [hours]'); grid on; box off;
legend('$u_6$','$u_7$','$u_8$','$u_9$','$u_{10}$')
% ax1(1) = figtex(gca,1); figsize(1,0.3);
title('Bath-Reactor proboelm state-variables: $u_1$,$u_2$,$u_3$,$u_4$,and $u_5$');
end


% %%
% % Xp' = Fx*Xp + Fp
% Xp0 = zeros(10,8); % initital condition
% 
% Fsens = @(t,Xp) fx(interp1(tsim,usim,t,'spline'),p').*squeeze(Xp)' + fp(interp1(tsim,usim,t,'spline')',p');
% 
% [tsens,Psens] = ode15s(@(t,y) Fsens(t,y),ts,Xp0,options);


% ask erik about this actually