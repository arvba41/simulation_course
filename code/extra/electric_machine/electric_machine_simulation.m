clear all
clc

%% parameters for the EM
params.J = 1.09; % induction motor inertia [kgm^2]
params.np = 20; % number of pole pairs
params.Nr = 1000; % output rotor speed [rpm] (nominal)
params.Ld = 48.351e-6; params.Lq = 48.351e-6; 
params.Rs = 0.0523;
params.psi = 0.12607*3.0298;

freq = params.np*params.Nr/60; % electrical frequency

K = 1;

v_pk = 800;
T_l = 0;

% w_ref = params.Nr*pi/30*params.np;
%% function definitions

% input funcrion
    va = @(t) v_pk*sin(2*pi*freq*t); 
    vb = @(t) v_pk*sin(2*pi*freq*t + 2*pi/3);
    vc = @(t) v_pk*sin(2*pi*freq*t + 4*pi/3);    
% Clarke tranformation
    valpha = @(t) 2/3*K*(va(t) - 1/2*vb(t)  - 1/2*vc(t));
    vbeta = @(t) 2/3*K*(sqrt(3)/2*vb(t)  - sqrt(3)/2*vc(t));
% Park tranformation
    vd = @(t) cos(2*pi*freq*t).*valpha(t) - sin(2*pi*freq*t).*vbeta(t);
    vq = @(t) cos(2*pi*freq*t).*vbeta(t) + sin(2*pi*freq*t).*valpha(t);
% function fy = [id; iq; wr]
fy = @(t,y) [1/params.Ld*(vd(t)-params.Rs*y(1)+params.Lq*y(2)*y(3));...
    1/params.Lq*(vq(t)-params.Rs*y(2)-(params.Ld*y(1) + params.psi)*y(3));...
    1/params.J*(3*params.np/(2*K^2)*((params.Ld - params.Lq)*y(1)*y(2)+...
        y(2)*params.psi)-T_l)];

y0 = [0;0;0]; % initial condition

%% Solution
ts = 0;
h = 0.0001;
tf = 1;

t = ts:h:tf;
%%%% Euler forward method
tic; yf = ef(fy,y0,t,h); toc ;
%%%% Euler backward method
% tic; yb = eb(fy,y0,t,h); toc ;
%%%% Trpizoidal method
% tic; ytrapz = trapz(fy,y0,t,h); toc ;

%%%% ode23s
% tic; [time_ode23s,yode23s] = ode23s(fy,t,y0); toc ;
%%% ode45
tic; [time_ode45,yode45] = ode45(fy,t,y0); toc ;
%%%% ode113
% tic; [time_ode15s,yode15s] = ode15s(fy,t,y0); toc ;

is_f = abs(yf(1,:) + 1j*yf(2,:));
% is_b = abs(yb(1,:) + 1j*yb(2,:));
% is_trapz = abs(ytrapz(1,:) + 1j*ytrapz(2,:));
is_ode45 = abs(yode45(:,1) + 1j*yode45(:,2));

iabc_f = dq2abc(yf(1:2,:),t,freq);
% iabc_b = dq2abc(yb(1:2,:),t,f);
% iabc_trapz = dq2abc(ytrapz(1:2,:),t,f);
iabc_odc45 = dq2abc(yode45(:,1:2)',t,freq);

% calculating the torques
Te_f = 3*params.np/(2*K^2)*((params.Ld - params.Lq)*yf(1,:).*yf(2,:) + yf(2,:)*params.psi);
% Te_b = 3*params.np/(2*K^2)*((params.Ld - params.Lq)*yb(1,:).*yb(2,:) + yb(2,:)*params.psi);
% Te_trapz = 3*params.np/(2*K^2)*((params.Ld - params.Lq)*ytrapz(1,:).*ytrapz(2,:) + ytrapz(2,:)*params.psi);
Te_ode45 = 3*params.np/(2*K^2)*((params.Ld - params.Lq)*yode45(:,1).*yode45(:,2) + yode45(:,2)*params.psi);

% inputs
% vabc = v_pk*cos(2*pi*f*t + [0;4*pi/3;2*pi/3]);
vabc = [va(t);vb(t);vc(t)];
% vdq0 = abc2dq(vabc,t,f,K);
vab0 = [valpha(t);vbeta(t)];
vdq0 = [vd(t);vq(t)];

vs = vd(t) + 1j*vq(t);

% speed conversion
wout = yf(3,:)*30/pi./params.np; 
wout_yode45 = yode45(:,3)*30/pi/params.np;
%% Plots
figure(1)
clf

tiledlayout(2,2)

nexttile;
plot(t,is_f,time_ode45,is_ode45)
% plot(t,iabc_f)
grid on
ylabel('$i_{dq}(t)$ [A]')
xlabel('$t$')
legend('ef','ode45')
ax1 = figtex(gca,1);

nexttile;
plot(t,wout,time_ode45,wout_yode45)
% plot(t,0*t + w_ref*30/pi)
% plot(t,iabc_f)
grid on
ylabel('$\omega_{r(e)}(t)$ [rpm]')
xlabel('$t$')
ax2 = figtex(gca);

nexttile;
plot(t,Te_f,time_ode45,Te_ode45)
grid on
ylabel('$T_e(t)$ [Nm]')
xlabel('$t$')
ax3 = figtex(gca);

nexttile;
plot(t,abs(vs))
grid on
ylabel('$v_{dq}(t)$')
xlabel('$t$')
ax4 = figtex(gca);

linkaxes([ax1 ax2 ax3 ax4],'x')