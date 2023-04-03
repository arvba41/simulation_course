clear all; clc

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

y01 = [0;0;0]; % initial condition

%% Solution
ts = 0;
h = 1e-6/5; h1 = 1e-8/5;
tf = 5;

t = ts:h:tf;

%%%% Euler forward method
tic; [yf,Ny] = ef_med_Nd(fy,y0,t,h); toc ;
tic; [yf1,Ny1] = ef_med_Nd(fy,y0,t,h1); toc ;

%%% ode45
tic; [time_ode45,yode45] = ode45(fy,t,y0); toc ;

tic; [time_ode451,yode451] = ode45(fy,t,y0); toc ;

% simulation data interption.
[is_f,Te_f,wout_ef] = em_readable_params(yf,params,K); 
[is_ode45,Te_ode45,wout_yode45] = em_readable_params_odeMatlab(yode45,params,K); 

[is_f1,Te_f1,wout_ef1] = em_readable_params(yf1,params,K); 
[is_ode451,Te_ode451,wout_yode451] = em_readable_params_odeMatlab(yode451,params,K); 

% input vectors
% vabc = [va(t);vb(t);vc(t)];
% vab0 = [valpha(t);vbeta(t)];
% vdq0 = [vd(t);vq(t)];
vs = vd(t) + 1j*vq(t);

%% Plots
figure(1)
clf

tiledlayout(2,3)

nexttile;
plot(t,is_f,t,is_f1)
% plot(t,iabc_f)
grid on
ylabel('$|i_{dq}(t)|$ [A]')
xlabel('$t$'); 
lgd1 = char("$h=" + string(h) + "$"); lgd2 = char("$h=" + string(h1) + "$"); 
legend(lgd1,lgd2); ax1 = figtex(gca,1);

nexttile;
plot(t,wout_ef,t,wout_ef1)
% plot(t,0*t + w_ref*30/pi)
% plot(t,iabc_f)
grid on
ylabel('$\omega_{r(e)}(t)$ [rpm]')
xlabel('$t$')
ax2 = figtex(gca);

% nexttile;
% plot(t,Te_f)
% grid on
% ylabel('$T_e(t)$ [Nm]')
% xlabel('$t$')
% ax3 = figtex(gca);

nexttile;
plot(t,abs(vs))
grid on
ylabel('$|v_{dq}(t)|$ [V]')
xlabel('$t$')
ax3 = figtex(gca); ylim([790 810])

nexttile;
semilogy(t(2:end),Ny(1,:),t(2:end),Ny1(1,:));
grid on; ax4 = figtex(gca);
ylabel('$\mathcal{N}_h(i_d)_hi_d(t)$'); xlabel('$t$')

nexttile;
semilogy(t(2:end),Ny(2,:),t(2:end),Ny1(2,:));
grid on; ax5 = figtex(gca);
ylabel('$\mathcal{N}_h(i_q)_hi_q(t)$'); xlabel('$t$')

nexttile;
semilogy(t(2:end),Ny(3,:),t(2:end),Ny1(3,:));
grid on; ax6 = figtex(gca);
ylabel('$\mathcal{N}_h(w_m)_hw_m(t)$'); xlabel('$t$')

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ],'x')
linkaxes([ax4 ax5 ax6],'y')

%% The direction field
% y_vec = [linspace(0,6000,50); linspace(0,12000,50); linspace(0, 1500, 50)];
% t_vec = linspace(0,10,50); 
% 
% fid = @(t,y) 1/params.Ld*(vd(t)-params.Rs*y(1)+params.Lq*y(2)*y(3));
% 
% fiq = @(t,y) 1/params.Lq*(vq(t)-params.Rs*y(2) ... 
%     - (params.Ld*y(1) + params.psi)*y(3));
% 
% fw = @(t,y) 1/params.J*(3*params.np/(2*K^2) ... 
%     *((params.Ld - params.Lq)*y(1)*y(2) + y(2)*params.psi)-T_l);

%% Solution to the ODE variable step size
ts = 0;
tf = 1;

% h1 = 0.001; h2 = 0.0001; h3 = 0.00001;
% h_vec = [0.001, 0.0001, 0.00001];

h_vec = 10.^[-1:-1:-5];
figure(3); clf; tiledlayout(2,3)
for ii = 1:length(h_vec)
    t_h = ts:h_vec(ii):tf; % t2 = ts:h2:tf; 

    vs_h = vd(t_h) + 1j*vq(t_h); % vs_h2 = vd(t2) + 1j*vq(t2);

    %%%% Euler forward method
    tic; yf_h = ef(fy,y0,t_h,h_vec(ii)); toc ;
    % tic; yf_h2 = ef(fy,y0,t2,h2); toc ;

    % simulation data interption.
    [~,Te_fh,wout_fh] = em_readable_params(yf_h,params,K); 
    % [~,Te_fh2,~] = em_readable_params(yf_h2,params,K); 

    % plottings
    em_plots(t_h,yf_h,Te_fh,wout_fh,vs_h)
end
% em_plots(t2,yf_h2,Te_fh2,vs_h2)

%% Solution to the ODE variable initial condition
ts = 0;
tf = 1;

% h1 = 0.001; h2 = 0.0001; h3 = 0.00001;
h_val = 0.0001;

% h_vec = 10.^[-1:-1:-5];
y0_vec = [0,0,0; 100,0,0; 0,100,0; 0,0,500; 500,500,0; 500,500,1500];

[len_y0_vec,~] = size(y0_vec);
figure(4); clf; tiledlayout(2,3)
for ii = 1:len_y0_vec
    t_h = ts:h_val:tf; % t2 = ts:h2:tf; 

    vs_h = vd(t_h) + 1j*vq(t_h); % vs_h2 = vd(t2) + 1j*vq(t2);

    %%%% Euler forward method
    tic; yf_h = ef(fy,y0_vec(ii,:)',t_h,h_val); toc ;
    % tic; yf_h2 = ef(fy,y0,t2,h2); toc ;

    % simulation data interption.
    [~,Te_fh,wout_fh] = em_readable_params(yf_h,params,K); 
    % [~,Te_fh2,~] = em_readable_params(yf_h2,params,K); 

    % plottings
    em_plots(t_h,yf_h,Te_fh,wout_fh,vs_h)
end
% em_plots(t2,yf_h2,Te_fh2,vs_h2)