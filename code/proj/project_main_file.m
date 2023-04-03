clear all
clc

% function definition: y'(t) = -10*y(t), y(0) = 1

ts = 0; % start time
tf = 1; % finish time

lambda = 10; 
fy = @(t,y) -lambda*y; % ode
trusol = @(t) exp(-lambda*t);
%% fixed step-size method 

% varying the step-size
hstart = 1; hend = 7;
h_vec = 10.^(-1*(hstart:hend));

% varying initial conditions
y0start = 1; y0end = 1;
y0_vec = y0start:0.1:y0end;

ysimef = cell(length(h_vec),length(y0_vec));
ferroref = zeros(length(h_vec),length(y0_vec));
% simulation and plots
for ii=1:length(h_vec)
    for jj=1:length(y0_vec)
        %% simulation
        tsim = ts:h_vec(ii):tf;
        yf = ef(@(t,y) fy(t,y),y0_vec(jj),tsim,h_vec(ii));
        
        %% difference method (ef)
        dn = ef_diff(@(t,y) fy(t,y),tsim,yf,h_vec(ii));
        
        %% data storage
        ysimef{ii,jj}.y = yf;
        ysimef{ii,jj}.t = tsim;
        ysimef{ii,jj}.dn = dn;
        ysimef{ii,jj}.ferror = abs(yf(end) - trusol(tsim(end))); 
        ferroref(ii,jj) = ysimef{ii,jj}.ferror;
        %% clean up
        clear yf tsim
    end
end
clear ii jj 

% plots
figure(1); clf;
% for ii = 1:length(y0_vec)
%     subplot(221);
%     plot(ysim{2,ii}.t,ysim{2,ii}.y,'LineWidth',2); grid on; hold on; 
%     
%     subplot(223);
%     semilogy(ysim{2,ii}.t,abs(ysim{2,ii}.y-trusol(ysim{2,ii}.t)),'LineWidth',2); grid on; hold on;
% end
% subplot(221); lgd = "$y_n(0)=\ $"+string(y0_vec(:));
% legend(char(lgd(1)),char(lgd(2)),char(lgd(3)),char(lgd(4)),char(lgd(5)),'Location','best'); 
% figtex(gca,1); xlabel('$t$ [s]'); ylabel('$y_n(t)$');
% subplot(223); 
% legend(char(lgd(1)),char(lgd(2)),char(lgd(3)),char(lgd(4)),char(lgd(5)),'Location','best'); 
% figtex(gca,1); xlabel('$t$ [s]'); ylabel('$|y_n - y(t)|$')

for ii = 1:length(h_vec)
    subplot(121);
    plot(ysimef{ii,1}.t,ysimef{ii,1}.y,'LineWidth',2); grid on; hold on; 
    
    subplot(122);
    semilogy(ysimef{ii,1}.t,abs(ysimef{ii,1}.y-trusol(ysimef{ii,1}.t)),'LineWidth',2); grid on; hold on;
end
subplot(121); lgd = "$h=\ $"+string(h_vec(:));
legend(char(lgd(1)),char(lgd(2)),char(lgd(3)),char(lgd(4)),char(lgd(5)),...
    char(lgd(6)),char(lgd(7)),'Location','best'); 
figtex(gca,1); xlabel('$t$ [s]'); ylabel('$y_n(t)$');
subplot(122); 
% legend(char(lgd(1)),char(lgd(2)),char(lgd(3)),char(lgd(4)),char(lgd(5)),...
%     char(lgd(6)),char(lgd(7)),'Location','best'); 
figtex(gca); xlabel('$t$ [s]'); ylabel('$|y_n - y(t)|$')

saveas(gcf,'../Figures/ef_testCase','epsc')

% figure(3); clf;


%% variable step method

y0ode45 = 1; % Initial value

% varying relative tolarance 
RTOLstart = 2; RTOLend = 8;
RTOLvec = 10.^(-1*(RTOLstart:RTOLend));

% varying absalute tolarance 
ATOLstart = 2; ATOLend = 8;
ATOLvec = 10.^(-1*(ATOLstart:ATOLend));

ysimode45 = cell(length(RTOLvec),length(ATOLvec));
ferrorode45 = zeros(length(RTOLvec),length(ATOLvec));
hode45 = zeros(length(RTOLvec),length(ATOLvec));
for ii=1:length(RTOLvec)
    for jj=1:length(ATOLvec)
        %% simulation
        ode45option = odeset('RelTol',RTOLvec(ii),'AbsTol',ATOLvec(jj));
        [tfode45,yfode45] = ode45(@(t,y) fy(t,y),[ts tf],y0ode45,ode45option);
        %% data storage
        ysimode45{ii,jj}.y = yfode45;
        ysimode45{ii,jj}.t = tfode45;
        ysimode45{ii,jj}.ferror = abs(yfode45(end) - trusol(tfode45(end))); 
        ferrorode45(ii,jj) = ysimode45{ii,jj}.ferror;
        hode45(ii,jj) = tfode45(end) - tfode45(end - 1); 
        %% clean up
        clear yf tsim ode45option
    end
end

% plots
figure(2); clf;
for ii = 1:length(RTOLvec)
    subplot(321);
    plot(ysimode45{ii,4}.t,ysimode45{ii,4}.y,'LineWidth',2); grid on; hold on; 
    
    subplot(323);
    semilogy(ysimode45{ii,4}.t(2:end),diff(ysimode45{ii,4}.t),'LineWidth',2); grid on; hold on;
    
    subplot(325);
    semilogy(ysimode45{ii,4}.t,abs(ysimode45{ii,4}.y-trusol(ysimode45{ii,4}.t)),'LineWidth',2); grid on; hold on;
end
subplot(321); lgd = "RTOL$=\ $"+string(RTOLvec(:));
legend(char(lgd(1)),char(lgd(2)),char(lgd(3)),char(lgd(4)),char(lgd(5)),...
    char(lgd(6)),char(lgd(7)),'Location','best'); 
figtex(gca,1); xlabel('$t$ [s]'); ylabel('$y_n(t)$'); 
title('Change in relative tolerance with ATOL = 1e$-$5');
subplot(323); figtex(gca); xlabel('$t$ [s]'); ylabel('$h$')
subplot(325); 
% legend(char(lgd(1)),char(lgd(2)),char(lgd(3)),char(lgd(4)),char(lgd(5)),...
%     char(lgd(6)),char(lgd(7)),'Location','best'); 
figtex(gca); xlabel('$t$ [s]'); ylabel('$|y_n - y(t)|$')

for ii = 1:length(ATOLvec)
    subplot(322);
    plot(ysimode45{4,ii}.t,ysimode45{4,ii}.y,'LineWidth',2); grid on; hold on; 
    
    subplot(324);
    semilogy(ysimode45{4,ii}.t(2:end),diff(ysimode45{4,ii}.t),'LineWidth',2); grid on; hold on;
    
    subplot(326);
    semilogy(ysimode45{4,ii}.t,abs(ysimode45{4,ii}.y-trusol(ysimode45{4,ii}.t)),'LineWidth',2); grid on; hold on;
end
subplot(322); lgd = "ATOL$=\ $"+string(ATOLvec(:));
legend(char(lgd(1)),char(lgd(2)),char(lgd(3)),char(lgd(4)),char(lgd(5)),...
    char(lgd(6)),char(lgd(7)),'Location','best'); 
figtex(gca,1); xlabel('$t$ [s]'); ylabel('$y_n(t)$');
title('Change in absolute tolerance with RTOL = 1e$-$5');
subplot(324); figtex(gca); xlabel('$t$ [s]'); ylabel('$h$')
subplot(326); 
% legend(char(lgd(1)),char(lgd(2)),char(lgd(3)),char(lgd(4)),char(lgd(5)),...
%     char(lgd(6)),char(lgd(7)),'Location','best'); 
figtex(gca); xlabel('$t$ [s]'); ylabel('$|y_n - y(t)|$')

saveas(gcf,'../Figures/ode45_testCase','epsc')

%% Shit        
%         subplot(212); 
%         semilogy(tsim,abs(yf-trusol(tsim)),'LineWidth',2); grid on; hold on; figtex(gca);
        

% %% Solving of equations
% 
% %%%%% Time step 1
% h = 1e-6; % time step
% t = ts:h:tf; % time vector 
% %%% Euler forward method
% tic; [yf,Ny] = ef_med_Nd(@(t,y) -10*y,1,t,h); toc ;
% 
% figure(10); clf;
% subplot(411); plot(t,yf,t,exp(-10*t)); grid on;
% ylabel('$y$'); xlabel('$t$'); 
% lgd1 = char("$h=" + string(h) + "$"); lgd2 = char("$h=\infty$"); 
% legend(lgd1,lgd2); figtex(gca,1); 
% title('Simulation Results')
% 
% subplot(412); semilogy(t(2:end),Ny); grid on;
% figtex(gca); ylabel('$\mathcal{N}_n(y)_ny(t)$'); xlabel('$t$'); 
% title('Local Trucation Error')
% 
% subplot(413); plot(t,abs(yf-exp(-10*t))); grid on; hold on;
% scatter(t(end),abs(yf(end)-exp(-10*t(end)))); hold off;
% fprintf("Global error: %e \n",abs(yf(end)-exp(-10*t(end))));
% figtex(gca); ylabel('$y_n - y(t)$'); xlabel('$t$'); 
% 
% %% Solving of equations
% 
% %%%%% Time step 1
% h = 1e-5; % time step
% t = ts:h:tf; % time vector 
% %%% Euler forward method
% tic; [yf,Ny] = ef_med_Nd(@(t,y) -1*y,1,t,h); toc ;
% %%% Rk 45 method 
% tic; sol_ode45 = ode45(@(t,y) -1*y,[ts tf],1,h); toc ;
% 
% figure(11); clf;
% 
% subplot(221); plot(t,yf,t,exp(-1*t),'--'); grid on;
% ylabel('$y$'); xlabel('$t$'); 
% lgd1 = char("numerical"); lgd2 = char("analytical");  
% legend(lgd1,lgd2); figtex(gca,1); title('Euler forward')
% 
% subplot(223); semilogy(t,abs(yf-exp(-1*t))); grid on; hold on;
% % scatter(t(end),abs(yf(end)-exp(-1*t(end)))); hold off;
% % fprintf("Global error Euler forward: %f \n",abs(yf(end)-exp(-10*t(end))));
% % fprintf("Global error ode45: %f \n",abs(sol_ode45.y(end)-exp(-10*sol_ode45.x(end))));
% figtex(gca); ylabel('$y_n - y(t)$'); xlabel('$t$'); 
% 
% subplot(222); plot(sol_ode45.x,sol_ode45.y,sol_ode45.x,exp(-1*sol_ode45.x),'--'); grid on;
% ylabel('$y$'); xlabel('$t$'); 
% lgd1 = char("numerical"); lgd2 = char("analytical");  
% legend(lgd1,lgd2); figtex(gca,1); title('ode45')
% 
% subplot(224); semilogy(sol_ode45.x,abs(sol_ode45.y-exp(-1*sol_ode45.x))); grid on; hold on;
% % scatter(t(end),abs(yf(end)-exp(-1*t(end)))); hold off;
% % fprintf("Global error Euler forward: %f \n",abs(yf(end)-exp(-10*t(end))));
% % fprintf("Global error ode45: %f \n",abs(sol_ode45.y(end)-exp(-10*sol_ode45.x(end))));
% figtex(gca); ylabel('$y_n - y(t)$'); xlabel('$t$'); 
