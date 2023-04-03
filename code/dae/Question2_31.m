clc; clear all;
% tol. parameters
RelTol = 1e-16;
AbsTol = 1e-16;
%% Analytical solution 
syms y(t) theta
y0 = 2; theta_0 = 1;
eqn = diff(y,t) + theta*y == 0;
cond = y(0) == y0;
S = simplify(dsolve(eqn,cond));

y_value_theta_0 = @(time_vector) subs(S,{theta,t},{theta_0,time_vector}); % true solution
y_value = @(time_vector) subs(S,{t},{time_vector}); % true solution with theta variable

syms G(theta) t 
T = 5;
dG_fun = diff(int(y_value(t)^2,t,[0,T]),theta); % diff(G,theta) G = int(y_value(t)^2,t,[0,T])
dG_value = subs(diff(int(y_value(t)^2,t,[0,T]),theta),{theta},{theta_0});

%% Forward way
clear theta
y0 = [2, 0, 0]; theta = 1;
% y = @(t,y) -theta*y; % X' = f(t,X,P);
% 
% fx = -theta; fp = -y; Xp = -theta*y0*exp(-theta*t);s
% Xp' = fx*Xp + fp
% 
% g = y^2 = yo*exp(-theta*t);
% dG/dp = int(Xp*gx+gp,t);
% gx = 2*y; gp = -theta*y0*exp(-theta*t)

f = @(t,y) [-theta*y(1); ... y
    -theta*y(2) - y(1);... Xp
    2*y(1)*y(2)... dG/dp
    ]; 

options = odeset('RelTol',RelTol,'AbsTol',AbsTol); % setting finer tolerences
[t,yval] = ode45(@(t,y) f(t,y),[0 5],y0);

figure(1)
clf

subplot(131)
plot(t,yval(:,1))
grid on; a(1) = figtex(gca);
xlabel('$t$'); ylabel('$y$');

subplot(132)
plot(t,yval(:,2))
grid on; a(2) = figtex(gca);
xlabel('$t$'); ylabel('$x_p$');

subplot(133)
plot(t,yval(:,3))
grid on; a(3) = figtex(gca);
xlabel('$t$'); ylabel('$dG/d\theta$');

figsize(1,0.25);
saveas(gcf,'Figures/Ugf21b','epsc');

%% Adjoint way
% discuss with Jian and gang
T = 5; % final time 
theta = 1; y0 = 2;

f1 = @(t,y) -theta*y; ... y

options = odeset('RelTol',RelTol,'AbsTol',AbsTol); % TOL setting
[tsim,yval_2] = ode45(@(t,y) f1(t,y),[0 5],y0,options); % solving for y

% F = y' + theta*y;
% Fx' = 1; Fx = theta;
% g = y^2 --> gx = 2y;
% (-lambda*Fx')' + lambda*Fx + gx = 0; % backwards from T
% defining z(tau) = x(T - tau), t = T - tau;
% x'(t) = -y(T-t) = -y'(tau); x(t) = y(T-tau) = y(tau)
% here, -y(tau) = -lambda(T-tau); y(tau) = theta*lambda(T-tau);
% gx = 2y(T-tau); % using interplation 

gx = @(t) interp1(tsim,2*yval_2,t,'spline'); % 2*y(tau) //ERROR should be Tau not t
f2 = @(tau,lambda) -theta*lambda + gx(T - tau); % DE for lambda

options = odeset('RelTol',RelTol,'AbsTol',AbsTol); % Setting TOL
lambda_fin = 0; % initital guess based on lambda*Fx' = 0; since Fx' is 1, lambda = 0.
[tsim_lambda,lambda_val] = ode45(@(t,y) f2(t,y),[0 5],lambda_fin,options); % solving for lambda

tsim_lambda_t = flipud(T - tsim_lambda);
lambda_t = flipud(lambda_val)'; % flipping the lambda and making it w.r.t t

% lambda(0)*Fx(0)*Xp(0)
% Xp(0) = 0 // I think not super sure
% lambda(0) = lambda_t(0)
% Fx(0) = theta
x0 = 1*theta*lambda_t(1);% Fx'(0)*Fx(0)*lambda(0)

% gp = -theta*y0*exp(-theta*t); Fp = y;
% function dG/dP = int(gp - lambda(t)*Fp)dt

options = odeset('RelTol',RelTol,'AbsTol',AbsTol);
f3 = @(t,x) -theta*interp1(tsim,yval_2,t,'linear') - interp1(tsim_lambda_t,lambda_t,t,'spline')*interp1(tsim,yval_2,t,'linear'); % dG/dp
[tsim_dGdp,dGdp_val] = ode45(@(t,y) f3(t,y),[0 5],x0,options); % just integrating


%% plots
figure(2)
clf

subplot(131)
plot(tsim,yval_2)
grid on; a(1) = figtex(gca);
xlabel('$t$'); ylabel('$y$');

subplot(132)
plot(tsim_lambda,lambda_t)
grid on; a(1) = figtex(gca);
xlabel('$t$'); ylabel('$\lambda^T(t)$');

subplot(133)
plot(tsim_dGdp,dGdp_val)
grid on; a(1) = figtex(gca);
ylabel('$\frac{dG}{d\theta}$'); xlabel('$t$');

figsize(1,0.25);
saveas(gcf,'Figures/Ugf31c','epsc');
%% Display the final sensitivities
fprintf("Analytical method dG/dP = %4.4f\n",double(dG_value))
fprintf("Forward method dG/dP = %4.4f\n",yval(end,end))
fprintf("Forward method dG/dP = %4.4f\n",dGdp_val(end))
