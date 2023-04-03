clear all; clc

%% parameters for the EM
params.J = 1.09; % induction motor inertia [kgm^2]
params.np = 10; % number of pole pairs
params.Nr = 1000; % output rotor speed [rpm] (nominal)
params.Ld = 48.351e-6; params.Lq = 48.351e-6; 
params.Rs = 0.0523;
params.psi = 0.12607*3.0298*2;

freq = params.np*params.Nr/60; % electrical frequency

K = 1;

v_pk = 800;
T_l = 0;

ts = 0; tf = 0.2; % simulation start and stop time

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

% y01 = [0;0;0]; % initial condition

%% fixed step method

% varying the step-size
hstart = 1; hend =7;
h_vec = 10.^(-1*(hstart:hend));

% varying initial conditions
y0start = [-10;-10]; y0end = [10;10];
y0_vec = [y0start(1):5:y0end(1);...
    y0start(2):5:y0end(2);...
    0 0 0 0 0];

ysimef = cell(length(h_vec),length(y0_vec));
ferroref = zeros(length(h_vec),length(y0_vec));

jjlen = size(y0_vec);
% simulation and plots
for ii=1:length(h_vec)
    for jj=1:jjlen(2)
        %% simulation
        tsim = ts:h_vec(ii):tf;
        tStart = tic; 
        yf = ef(@(t,y) fy(t,y),y0_vec(:,jj),tsim,h_vec(ii)); 
        tEnd = toc(tStart);
        %% data manupulation
        fprintf('time taken for h = %g, inter no. %d is %.4f s.\n',...
            h_vec(ii),jj,tEnd);
        [is,Te,wout] = em_readable_params(yf,params,K);
        %% data storage
        ysimef{ii,jj}.y = yf; ysimef{ii,jj}.is = is; ysimef{ii,jj}.Te= Te;
        ysimef{ii,jj}.wout = wout; ysimef{ii,jj}.t = tsim;
%         ysimef{ii,jj}.ferror = abs(yf(end) - trusol(tsim(end))); 
%         ferroref(ii,jj) = ysimef{ii,jj}.ferror;
        %% clean up
        clear yf tsim is Te wout
    end
end
clear yf tsim is Te wout
clear ii jj jjlen

%% plots
clear ax 

figure(3); clf
count = 1;
tiledlayout(3,3)
for ii = [3 5 7]
    nexttile;
    plot(ysimef{ii,3}.t,ysimef{ii,3}.y(1,:),'LineWidth',2); grid on;  
    title('$i_{d}(t)$ [A]'); xlabel('$t$ [s]');
    ax(count) = figtex(gca);
    
    nexttile;
    plot(ysimef{ii,3}.t,ysimef{ii,3}.y(2,:),'LineWidth',2); grid on;  
    title('$i_{q}(t)$ [A]'); xlabel('$t$ [s]'); 
    ax(count + 1) = figtex(gca);

    nexttile;
    plot(ysimef{ii,3}.t,ysimef{ii,3}.wout,'LineWidth',2); grid on;
    title('$\omega_{r(m)}(t)$ [rpm]'); xlabel('$t$ [s]'); 
    lgd = "$h=\ $"+string(h_vec(ii));
    legend(char(lgd)); 
    ax(count + 2) = figtex(gca,1);
    
    count = count + 3;
    
end
clear count ii

linkaxes(ax,'x');
linkaxes([ax(1) ax(4) ax(7)],'y');
linkaxes([ax(2) ax(5) ax(8)],'y');
linkaxes([ax(3) ax(6) ax(9)],'y'); 
ax(1).YLim = [0 4e3]; ax(2).YLim = [0 4e3];

%% plots 2.0


figure(4); clf
count = 1;
tiledlayout(4,3)
for ii = [3 4 5 6]
    nexttile;
    plot(ysimef{ii,3}.t,ysimef{ii,3}.y(1,:) ...
        - interpn(ysimef{end,3}.t,ysimef{end,3}.y(1,:),ysimef{ii,3}.t),...
        'LineWidth',2); grid on;  
    title('$i_{d}(t)$ [A]'); xlabel('$t$ [s]');
    ax(count) = figtex(gca);
    
    nexttile;
    plot(ysimef{ii,3}.t,ysimef{ii,3}.y(2,:) ...
        - interpn(ysimef{end,3}.t,ysimef{end,3}.y(2,:),ysimef{ii,3}.t)...
        ,'LineWidth',2); grid on;  
    title('$i_{q}(t)$ [A]'); xlabel('$t$ [s]'); 
    ax(count + 1) = figtex(gca);

    nexttile;
    plot(ysimef{ii,3}.t,ysimef{ii,3}.wout ...
        - interpn(ysimef{end,3}.t,ysimef{end,3}.wout,ysimef{ii,3}.t)...
        ,'LineWidth',2); grid on;
    title('$\omega_{r(m)}(t)$ [rpm]'); xlabel('$t$ [s]'); 
    lgd = "$h=\ $"+string(h_vec(ii));
    legend(char(lgd)); 
    ax(count + 2) = figtex(gca,1);
    
    count = count + 3;
    
end
clear count ii

%%

subplot(331); ylim([0 1000]); xlim([0 5]);
% title(char(join(lgd)))


subplot(332); ylim([0 100]); xlim([0 5]);
subplot(333); ylim([0 1500]); xlim([0 5]); 
ax(9).XLim = [ts tf];
% saveas(gcf,'../Figures/ef_EMrunCase','epsc'); 
%
jjlen = size(y0_vec);
is_end = zeros(length(h_vec),length(y0_vec));
Te_end = zeros(length(h_vec),length(y0_vec));
wout_end = zeros(length(h_vec),length(y0_vec));
% simulation and plots
for ii=1:length(h_vec)
    for jj=1:jjlen(2)
        is_end(ii,jj) = ysimef{ii,jj}.is(end);
        Te_end(ii,jj) = ysimef{ii,jj}.Te(end);
        wout_end(ii,jj) = ysimef{ii,jj}.wout(end);
    end
end

%% 



% saveas(gcf,'../Figures/ef_MotorCase','epsc')
%% variable step method

% varying relative tolarance 
RTOLstart = -3; RTOLend = 5;
RTOLvec = 10.^(RTOLstart:RTOLend);

% varying absalute tolarance 
ATOLstart = -3; ATOLend = 5;
ATOLvec = 10.^(ATOLstart:ATOLend);


% varying initial conditions
y0start = [-10;-10]; y0end = [10;10];
% y0_vec = [y0start(1):5:y0end(1);...
%     y0start(2):5:y0end(2);...
%     0 0 0 0 0];
y0_vec = [0; 0; 0];

ysimode45 = cell(length(ATOLvec),length(RTOLvec));
ferroref = zeros(length(ATOLvec),length(RTOLvec));
for ii=1:length(ATOLvec)
    for jj=1:length(RTOLvec)
        %% simulation
        tStart = tic;
        ode45option = odeset('AbsTol',ATOLvec(ii),'RelTol',RTOLvec(jj));
        [tfode45,yfode45] = ode45(@(t,y) fy(t,y),[ts tf],y0_vec,ode45option);
        tEnd = toc(tStart); 
        %% data manupulation
        fprintf('time taken for ATol = %g and RTol = %g , is %.4f s.\n',...
            ATOLvec(ii),RTOLvec(jj),tEnd);
        [is,Te,wout] = em_readable_params_odeMatlab(yfode45,params,K);
        %% data storage
        ysimode45{ii,jj}.y = yfode45; ysimode45{ii,jj}.is = is; ysimode45{ii,jj}.Te= Te;
        ysimode45{ii,jj}.wout = wout; ysimode45{ii,jj}.t = tfode45;
%         ysimef{ii,jj}.ferror = abs(yf(end) - trusol(tsim(end))); 
%         ferroref(ii,jj) = ysimef{ii,jj}.ferror;
        %% clean up
        clear yf tsim is Te wout
    end
end

%% plots
clear ax 

figure(4); clf
count = 1;
for ii = [1 5 9]
    subplot(3,4,count);
    plot(ysimode45{ii,ii}.t,real(ysimode45{ii,ii}.is),'LineWidth',2); grid on;  
    title('$i_{d}(t)$ [A]'); xlabel('$t$ [s]');
    ax(count) = figtex(gca);
    
    subplot(3,4,count + 1);
    plot(ysimode45{ii,ii}.t,imag(ysimode45{ii,ii}.is),'LineWidth',2); grid on;  
    title('$i_{q}(t)$ [A]'); xlabel('$t$ [s]'); 
    lgd = "ATol = "+string(ATOLvec(ii));
    legend(char(lgd)); 
    ax(count + 1) = figtex(gca,1);

    subplot(3,4,count + 2);
    plot(ysimode45{ii,ii}.t,ysimode45{ii,ii}.wout,'LineWidth',2); grid on;
    title('$\omega_{r(m)}(t)$ [rpm]'); xlabel('$t$ [s]'); 
    ax(count + 2) = figtex(gca);
    
    subplot(3,4,count + 3);
    plot(ysimode45{ii,ii}.t(2:end),diff(ysimode45{ii,ii}.t),'LineWidth',2); grid on;
    title('$h$ [-]'); xlabel('$t$ [s]'); 
    ax(count + 3) = figtex(gca);

    count = count + 4;
    
end
clear count ii


linkaxes(ax,'x');
linkaxes([ax(1) ax(5) ax(9)],'y');
linkaxes([ax(2) ax(6) ax(10)],'y');
linkaxes([ax(3) ax(7) ax(11)],'y');
linkaxes([ax(4) ax(8) ax(12)],'y');

% subplot(331); ylim([0 1000]); xlim([0 5]);
% title(char(join(lgd)))


% subplot(332); ylim([0 100]); xlim([0 5]);
% subplot(333); ylim([0 1500]); xlim([0 5]); 
ax(9).XLim = [ts tf];
% ax(8).YLim = [1e-5 1];
ax(7).YLim = [0 1500];
ax(1).YLim = [0 1000];
ax(2).YLim = [0 1000];
% saveas(gcf,'../Figures/ef_EMrunCase','epsc'); 
%
jjlen = size(y0_vec);
is_end = zeros(length(h_vec),length(y0_vec));
Te_end = zeros(length(h_vec),length(y0_vec));
wout_end = zeros(length(h_vec),length(y0_vec));
% simulation and plots
for ii=1:length(h_vec)
    for jj=1:jjlen(2)
        is_end(ii,jj) = ysimode45{ii,jj}.is(end);
        Te_end(ii,jj) = ysimode45{ii,jj}.Te(end);
        wout_end(ii,jj) = ysimode45{ii,jj}.wout(end);
    end
end

% saveas(gcf,'../Figures/ef_MotorCase_ode45','epsc')

%% selecting plots

figure(5); clf;
ii = 1; jj = 1; 

    subplot(5,4,1);
    plot(ysimode45{ii,jj}.t,real(ysimode45{ii,jj}.is),'LineWidth',2); grid on;  
    title('$i_{d}(t)$ [A]'); xlabel('$t$ [s]');
    ax(1) = figtex(gca); hold on
    
    subplot(5,4,1 + 1);
    plot(ysimode45{ii,jj}.t,imag(ysimode45{ii,jj}.is),'LineWidth',2); grid on;  
    title('$i_{q}(t)$ [A]'); xlabel('$t$ [s]'); 
    ax(1 + 1) = figtex(gca); hold on

    subplot(5,4,1 + 2);
    plot(ysimode45{ii,jj}.t,ysimode45{ii,jj}.wout,'LineWidth',2); grid on;
    title('$\omega_{r(m)}(t)$ [rpm]'); xlabel('$t$ [s]'); 
    ax(1 + 2) = figtex(gca); hold on
    
    subplot(5,4,1 + 3);
    semilogy(ysimode45{ii,jj}.t(2:end),diff(ysimode45{ii,jj}.t),'LineWidth',2); grid on;
    title('$h$ [-]'); xlabel('$t$ [s]'); 
    lgd = "ATOL = "+string(ATOLvec(ii)) + ", RTOL = " + string(RTOLvec(jj));
    legend(char(lgd),'Location','northwest'); 
    ax(1 + 3) = figtex(gca,1); hold on

ii = 1; jj = 4; 

    subplot(5,4,5);
    plot(ysimode45{ii,jj}.t,real(ysimode45{ii,jj}.is),'LineWidth',2); grid on;  
    title('$i_{d}(t)$ [A]'); xlabel('$t$ [s]');
    ax(5) = figtex(gca); hold on
    
    subplot(5,4,5 + 1);
    plot(ysimode45{ii,jj}.t,imag(ysimode45{ii,jj}.is),'LineWidth',2); grid on;  
    title('$i_{q}(t)$ [A]'); xlabel('$t$ [s]'); 
    ax(5 + 1) = figtex(gca); hold on

    subplot(5,4,5 + 2);
    plot(ysimode45{ii,jj}.t,ysimode45{ii,jj}.wout,'LineWidth',2); grid on;
    title('$\omega_{r(m)}(t)$ [rpm]'); xlabel('$t$ [s]'); 
    ax(5 + 2) = figtex(gca); hold on
    
    subplot(5,4,5 + 3);
    semilogy(ysimode45{ii,jj}.t(2:end),diff(ysimode45{ii,jj}.t),'LineWidth',2); grid on;
    title('$h$ [-]'); xlabel('$t$ [s]');
    lgd = "ATOL = "+string(ATOLvec(ii)) + ", RTOL = " + string(RTOLvec(jj));
    legend(char(lgd),'Location','northwest'); 
    ax(5 + 3) = figtex(gca,1); hold on
    
ii = 8; jj = 1; 

    subplot(5,4,9);
    plot(ysimode45{ii,jj}.t,real(ysimode45{ii,jj}.is),'LineWidth',2); grid on;  
    title('$i_{d}(t)$ [A]'); xlabel('$t$ [s]');
    ax(9) = figtex(gca); 
    
    subplot(5,4,9 + 1);
    plot(ysimode45{ii,jj}.t,imag(ysimode45{ii,jj}.is),'LineWidth',2); grid on;  
    title('$i_{q}(t)$ [A]'); xlabel('$t$ [s]');  
    ax(9 + 1) = figtex(gca); 

    subplot(5,4,9 + 2);
    plot(ysimode45{ii,jj}.t,ysimode45{ii,jj}.wout,'LineWidth',2); grid on;
    title('$\omega_{r(m)}(t)$ [rpm]'); xlabel('$t$ [s]'); 
    ax(9 + 2) = figtex(gca); 
    
    subplot(5,4,9 + 3);
    semilogy(ysimode45{ii,jj}.t(2:end),diff(ysimode45{ii,jj}.t),'LineWidth',2); grid on;
    title('$h$ [-]'); xlabel('$t$ [s]'); 
    lgd = "ATOL = "+string(ATOLvec(ii)) + ", RTOL = " + string(RTOLvec(jj));
    legend(char(lgd),'Location','northwest'); 
    ax(9 + 3) = figtex(gca,1); 
    
ii = 9; jj = 3; 

    subplot(5,4,13);
    plot(ysimode45{ii,jj}.t,real(ysimode45{ii,jj}.is),'LineWidth',2); grid on;  
    title('$i_{d}(t)$ [A]'); xlabel('$t$ [s]');
    ax(13) = figtex(gca); 
    
    subplot(5,4,13 + 1);
    plot(ysimode45{ii,jj}.t,imag(ysimode45{ii,jj}.is),'LineWidth',2); grid on;  
    title('$i_{q}(t)$ [A]'); xlabel('$t$ [s]');  
    ax(13 + 1) = figtex(gca); 

    subplot(5,4,13 + 2);
    plot(ysimode45{ii,jj}.t,ysimode45{ii,jj}.wout,'LineWidth',2); grid on;
    title('$\omega_{r(m)}(t)$ [rpm]'); xlabel('$t$ [s]'); 
    ax(13 + 2) = figtex(gca); 
    
    subplot(5,4,13 + 3);
    semilogy(ysimode45{ii,jj}.t(2:end),diff(ysimode45{ii,jj}.t),'LineWidth',2); grid on;
    title('$h$ [-]'); xlabel('$t$ [s]'); 
    lgd = "ATOL = "+string(ATOLvec(ii)) + ", RTOL = " + string(RTOLvec(jj));
    legend(char(lgd),'Location','northwest'); 
    ax(13 + 3) = figtex(gca,1); 

ii = 9; jj = 4; 

    subplot(5,4,17);
    plot(ysimode45{ii,jj}.t,real(ysimode45{ii,jj}.is),'LineWidth',2); grid on;  
    title('$i_{d}(t)$ [A]'); xlabel('$t$ [s]');
    ax(17) = figtex(gca); 
    
    subplot(5,4,17 + 1);
    plot(ysimode45{ii,jj}.t,imag(ysimode45{ii,jj}.is),'LineWidth',2); grid on;  
    title('$i_{q}(t)$ [A]'); xlabel('$t$ [s]');  
    ax(17 + 1) = figtex(gca); 

    subplot(5,4,17 + 2);
    plot(ysimode45{ii,jj}.t,ysimode45{ii,jj}.wout,'LineWidth',2); grid on;
    title('$\omega_{r(m)}(t)$ [rpm]'); xlabel('$t$ [s]'); 
    ax(17 + 2) = figtex(gca); 
    
    subplot(5,4,17 + 3);
    semilogy(ysimode45{ii,jj}.t(2:end),diff(ysimode45{ii,jj}.t),'LineWidth',2); grid on;
    title('$h$ [-]'); xlabel('$t$ [s]'); 
    lgd = "ATOL = "+string(ATOLvec(ii)) + ", RTOL = " + string(RTOLvec(jj));
    legend(char(lgd),'Location','northwest'); 
    ax(17 + 3) = figtex(gca,1);  
    
linkaxes(ax,'x'); ax(20).XLim = [ts tf];
linkaxes([ax(1) ax(5) ax(9) ax(13) ax(17)],'y');
linkaxes([ax(2) ax(6) ax(10) ax(14) ax(18)],'y');
linkaxes([ax(3) ax(7) ax(11) ax(15) ax(19)],'y');
linkaxes([ax(4) ax(8) ax(12) ax(16) ax(20)],'y');    

% ax(8).YLim = [1e-5 1];
ax(1).YLim = [0 1000];
ax(2).YLim = [0 1000];
ax(3).YLim = [0 2000];
ax(12).YLim = [1e-4 3e-3];
% ax(2).YLim = [0 1000];

% saveas(gcf,'../Figures/ef_MotorCase_ode45_Atol_Rtol','epsc')
% saveas(gcf,'../Figures/ef_MotorCase_ode45_Atol_Rtol','svg')

% %% Solution
% ts = 0;
% h = 1e-6/5; h1 = 1e-8/5;
% tf = 5;
% 
% t = ts:h:tf;
% 
% %%%% Euler forward method
% tic; [yf,Ny] = ef_med_Nd(fy,y0,t,h); toc ;
% tic; [yf1,Ny1] = ef_med_Nd(fy,y0,t,h1); toc ;
% 
% %%% ode45
% tic; [time_ode45,yode45] = ode45(fy,t,y0); toc ;
% 
% tic; [time_ode451,yode451] = ode45(fy,t,y0); toc ;
% 
% % simulation data interption.
% [is_f,Te_f,wout_ef] = em_readable_params(yf,params,K); 
% [is_ode45,Te_ode45,wout_yode45] = em_readable_params_odeMatlab(yode45,params,K); 
% 
% [is_f1,Te_f1,wout_ef1] = em_readable_params(yf1,params,K); 
% [is_ode451,Te_ode451,wout_yode451] = em_readable_params_odeMatlab(yode451,params,K); 
% 
% % input vectors
% % vabc = [va(t);vb(t);vc(t)];
% % vab0 = [valpha(t);vbeta(t)];
% % vdq0 = [vd(t);vq(t)];
% vs = vd(t) + 1j*vq(t);
% 
% %% Plots
% figure(1)
% clf
% 
% tiledlayout(2,3)
% 
% nexttile;
% plot(t,is_f,t,is_f1)
% % plot(t,iabc_f)
% grid on
% ylabel('$|i_{dq}(t)|$ [A]')
% xlabel('$t$'); 
% lgd1 = char("$h=" + string(h) + "$"); lgd2 = char("$h=" + string(h1) + "$"); 
% legend(lgd1,lgd2); ax1 = figtex(gca,1);
% 
% nexttile;
% plot(t,wout_ef,t,wout_ef1)
% % plot(t,0*t + w_ref*30/pi)
% % plot(t,iabc_f)
% grid on
% ylabel('$\omega_{r(e)}(t)$ [rpm]')
% xlabel('$t$')
% ax2 = figtex(gca);
% 
% % nexttile;
% % plot(t,Te_f)
% % grid on
% % ylabel('$T_e(t)$ [Nm]')
% % xlabel('$t$')
% % ax3 = figtex(gca);
% 
% nexttile;
% plot(t,abs(vs))
% grid on
% ylabel('$|v_{dq}(t)|$ [V]')
% xlabel('$t$')
% ax3 = figtex(gca); ylim([790 810])
% 
% nexttile;
% semilogy(t(2:end),Ny(1,:),t(2:end),Ny1(1,:));
% grid on; ax4 = figtex(gca);
% ylabel('$\mathcal{N}_h(i_d)_hi_d(t)$'); xlabel('$t$')
% 
% nexttile;
% semilogy(t(2:end),Ny(2,:),t(2:end),Ny1(2,:));
% grid on; ax5 = figtex(gca);
% ylabel('$\mathcal{N}_h(i_q)_hi_q(t)$'); xlabel('$t$')
% 
% nexttile;
% semilogy(t(2:end),Ny(3,:),t(2:end),Ny1(3,:));
% grid on; ax6 = figtex(gca);
% ylabel('$\mathcal{N}_h(w_m)_hw_m(t)$'); xlabel('$t$')
% 
% linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ],'x')
% linkaxes([ax4 ax5 ax6],'y')
% 
% %% The direction field
% % y_vec = [linspace(0,6000,50); linspace(0,12000,50); linspace(0, 1500, 50)];
% % t_vec = linspace(0,10,50); 
% % 
% % fid = @(t,y) 1/params.Ld*(vd(t)-params.Rs*y(1)+params.Lq*y(2)*y(3));
% % 
% % fiq = @(t,y) 1/params.Lq*(vq(t)-params.Rs*y(2) ... 
% %     - (params.Ld*y(1) + params.psi)*y(3));
% % 
% % fw = @(t,y) 1/params.J*(3*params.np/(2*K^2) ... 
% %     *((params.Ld - params.Lq)*y(1)*y(2) + y(2)*params.psi)-T_l);
% 
% %% Solution to the ODE variable step size
% ts = 0;
% tf = 1;
% 
% % h1 = 0.001; h2 = 0.0001; h3 = 0.00001;
% % h_vec = [0.001, 0.0001, 0.00001];
% 
% h_vec = 10.^[-1:-1:-5];
% figure(3); clf; tiledlayout(2,3)
% for ii = 1:length(h_vec)
%     t_h = ts:h_vec(ii):tf; % t2 = ts:h2:tf; 
% 
%     vs_h = vd(t_h) + 1j*vq(t_h); % vs_h2 = vd(t2) + 1j*vq(t2);
% 
%     %%%% Euler forward method
%     tic; yf_h = ef(fy,y0,t_h,h_vec(ii)); toc ;
%     % tic; yf_h2 = ef(fy,y0,t2,h2); toc ;
% 
%     % simulation data interption.
%     [~,Te_fh,wout_fh] = em_readable_params(yf_h,params,K); 
%     % [~,Te_fh2,~] = em_readable_params(yf_h2,params,K); 
% 
%     % plottings
%     em_plots(t_h,yf_h,Te_fh,wout_fh,vs_h)
% end
% % em_plots(t2,yf_h2,Te_fh2,vs_h2)
% 
% %% Solution to the ODE variable initial condition
% ts = 0;
% tf = 1;
% 
% % h1 = 0.001; h2 = 0.0001; h3 = 0.00001;
% h_val = 0.0001;
% 
% % h_vec = 10.^[-1:-1:-5];
% y0_vec = [0,0,0; 100,0,0; 0,100,0; 0,0,500; 500,500,0; 500,500,1500];
% 
% [len_y0_vec,~] = size(y0_vec);
% figure(4); clf; tiledlayout(2,3)
% for ii = 1:len_y0_vec
%     t_h = ts:h_val:tf; % t2 = ts:h2:tf; 
% 
%     vs_h = vd(t_h) + 1j*vq(t_h); % vs_h2 = vd(t2) + 1j*vq(t2);
% 
%     %%%% Euler forward method
%     tic; yf_h = ef(fy,y0_vec(ii,:)',t_h,h_val); toc ;
%     % tic; yf_h2 = ef(fy,y0,t2,h2); toc ;
% 
%     % simulation data interption.
%     [~,Te_fh,wout_fh] = em_readable_params(yf_h,params,K); 
%     % [~,Te_fh2,~] = em_readable_params(yf_h2,params,K); 
% 
%     % plottings
%     em_plots(t_h,yf_h,Te_fh,wout_fh,vs_h)
% end
% % em_plots(t2,yf_h2,Te_fh2,vs_h2)

%% 
figure(6); clf
ii = 5; jj = 3; 

    subplot(1,4,1);
    plot(ysimode45{ii,jj}.t,real(ysimode45{ii,jj}.is),'LineWidth',2); grid on;  
    title('$i_{d}(t)$ [A]'); xlabel('$t$ [s]');
    ax(1) = figtex(gca); hold on
    
    subplot(1,4,1 + 1);
    plot(ysimode45{ii,jj}.t,imag(ysimode45{ii,jj}.is),'LineWidth',2); grid on;  
    title('$i_{q}(t)$ [A]'); xlabel('$t$ [s]'); 
    ax(1 + 1) = figtex(gca); hold on

    subplot(1,4,1 + 2);
    plot(ysimode45{ii,jj}.t,ysimode45{ii,jj}.wout,'LineWidth',2); grid on;
    title('$\omega_{r(m)}(t)$ [rpm]'); xlabel('$t$ [s]'); 
    ax(1 + 2) = figtex(gca); hold on
    
    subplot(1,4,1 + 3);
    semilogy(ysimode45{ii,jj}.t(2:end),diff(ysimode45{ii,jj}.t),'LineWidth',2); grid on;
    title('$h$ [-]'); xlabel('$t$ [s]'); 
    lgd = "ATOL = "+string(ATOLvec(ii)) + ", RTOL = " + string(RTOLvec(jj));
    legend(char(lgd),'Location','northwest'); 
    ax(1 + 3) = figtex(gca,1); hold on
    
    ax(1).YLim = [0 1000];
    ax(2).YLim = [0 1000];
    ax(3).YLim = [0 2000];
%     ax(4).YLim = [1e-4 3e-3];
    linkaxes(ax,'x'); ax(1).XLim = [ts tf];