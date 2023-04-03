clear all
clc

params.J = 1.09; % induction motor inertia [kgm^2]
params.np = 20; % number of pole pairs
params.Nr = 1000; % output rotor speed [rpm] (nominal)
params.Ld = 48.351e-6; params.Lq = 48.351e-6; 
params.Rs = 0.0523;
params.psi = 0.12607*3.0298;

f = params.np*params.Nr/60; % electrical frequency

K = 1;

v_pk = 800;
T_l = 0;

% w_ref = params.Nr*pi/30*params.np;
%% function definitions

% input funcrion
    va = @(t) v_pk*sin(2*pi*f*t); vb = @(t) v_pk*sin(2*pi*f*t + 2*pi/3);
        vc = @(t) v_pk*sin(2*pi*f*t + 4*pi/3);
    % Clarke tranformation
    valpha = @(t) 2/3*K*(va(t) - 1/2*vb(t)  - 1/2*vc(t));
    vbeta = @(t) 2/3*K*(sqrt(3)/2*vb(t)  - sqrt(3)/2*vc(t));
    % Park tranformation
    vd = @(t) cos(2*pi*f*t).*valpha(t) - sin(2*pi*f*t).*vbeta(t);
    vq = @(t) cos(2*pi*f*t).*vbeta(t) + sin(2*pi*f*t).*valpha(t);
    % speed 
    w_em = 1000*pi/30*params.np;

% function fy = [id; iq; Te]
fy = @(t,y) [1/params.Ld*(vd(t)-params.Rs*y(1)+params.Lq*y(2)*w_em);...
    1/params.Lq*(vq(t)-params.Rs*y(2)-(params.Ld*y(1) + params.psi)*w_em);...
    3/(2*K^2)*w_em*(params.psi*y(2) + (params.Ld - params.Lq)*y(1)*y(2))];

y0 = [0;0;0];
%% Solution
ts = 0;
h = 0.0001;
tf = 0.1;

t = ts:h:tf;
%%%% Euler forward method
tic; yf = ef(fy,y0,t,h); toc ;
%%%% Euler backward method
tic; yb = eb(fy,y0,t,h); toc ;
%%%% Trpizoidal method
tic; ytrapz = trapz(fy,y0,t,h); toc ;

%%%% ode23s
% tic; [time_ode23s,yode23s] = ode23s(fy,t,y0); toc ;
%%% ode45
tic; [time_ode45,yode45] = ode45(fy,t,y0); toc ;
%%%% ode113
% tic; [time_ode15s,yode15s] = ode15s(fy,t,y0); toc ;

is_f = abs(yf(1,:) + 1j*yf(2,:));
is_b = abs(yb(1,:) + 1j*yb(2,:));
is_trapz = abs(ytrapz(1,:) + 1j*ytrapz(2,:));
is_ode45 = abs(yode45(:,1) + 1j*yode45(:,2));

iabc_f = dq2abc(yf(1:2,:),t,f);
iabc_b = dq2abc(yb(1:2,:),t,f);
iabc_trapz = dq2abc(ytrapz(1:2,:),t,f);
% iabc_odc45 = dq2abc(yode45(:,1:2),t,f);

Te_f = yf(3,:);
Te_b = yb(3,:);
Te_trapz = ytrapz(3,:);
Te_ode45 = yode45(:,3);

% inputs
% vabc = v_pk*cos(2*pi*f*t + [0;4*pi/3;2*pi/3]);
vabc = [va(t);vb(t);vc(t)];
% vdq0 = abc2dq(vabc,t,f,K);
vab0 = [valpha(t);vbeta(t)];
vdq0 = [vd(t);vq(t)];


%% Plots
figure(1)
clf

tiledlayout(2,2)

nexttile;
plot(t,is_f,t,is_b,t,is_trapz,time_ode45,is_ode45)
grid on
ylabel('$i_{dq}(t)$ [A]')
xlabel('$t$')
legend('ef','eb','trapz','ode45')
figtex(gca,1);

nexttile;
% plot(t,wout,time_ode45,wout_yode45)
plot(t,0*t + w_em*30/pi/params.np)
% plot(t,iabc_f)
grid on
ylabel('$\omega_{r}(t)$ [rpm]')
xlabel('$t$')
figtex(gca);

nexttile;
plot(t,Te_f,t,Te_b,t,Te_trapz,time_ode45,Te_ode45)
grid on
ylabel('$T_e(t)$ [Nm]')
xlabel('$t$')
figtex(gca);

nexttile;
plot(t,vdq0)
grid on
ylabel('$v_{dq}(t)$')
xlabel('$t$')
figtex(gca);