clear all
clc

% function definition: y'(t) = -10*y(t), y(0) = 1

ts = 1; % start time
tf = 25; % finish time

lambda = 1; 
fy = @(t,y) -5*t*y.^2 + 5/t - 1/(t^2); % ode
trusol = @(t) 1./t;

%% constant step-size method 

% varying the step-size
hstart = 0.19; hend = 0.21;
h_vec = 10.^(-1*[hstart:0.01:hend]);

% varying initial conditions
y0start = 1; y0end = 1.2;
y0_vec = y0start;

ysim = cell(length(h_vec),length(y0_vec));
% simulation and plots
for ii=1:length(h_vec)
    for jj=1:length(y0_vec)
        %% simulation
        tsim = ts:h_vec(ii):tf;
        yf = ef(@(t,y) fy(t,y),y0_vec(jj),tsim,h_vec(ii));
        %% data storage
        ysim{ii,jj}.y = yf;
        ysim{ii,jj}.t = tsim;
        ysim{ii,jj}.ferror = yf(end) - trusol(tsim(end)); 
        %% clean up
        clear yf tsim
    end
end
clear ii jj 

%% plots%
figure(10); clf;
for ii = 1:length(y0_vec)
    subplot(221);
    plot(ysim{2,ii}.t,ysim{2,ii}.y,'LineWidth',2); grid on; hold on; 
    
    subplot(223);
    semilogy(ysim{2,ii}.t,abs(ysim{2,ii}.y-trusol(ysim{2,ii}.t)),'LineWidth',2); grid on; hold on;
end
subplot(221); lgd = "$y_n(0)=\ $"+string(y0_vec(:));
% legend(char(lgd(1)),char(lgd(2)),char(lgd(3)),char(lgd(4)),char(lgd(5))); 
figtex(gca); xlabel('$t$ [s]'); ylabel('$y_n(t)$');
subplot(223); 
% legend(char(lgd(1)),char(lgd(2)),char(lgd(3)),char(lgd(4)),char(lgd(5))); 
figtex(gca); xlabel('$t$ [s]'); ylabel('$|y_n - y(t)|$')

for ii = 1:length(h_vec)
    subplot(222);
    plot(ysim{ii,1}.t,ysim{ii,1}.y,'LineWidth',2); grid on; hold on; 
    
    subplot(224);
    semilogy(ysim{ii,1}.t,abs(ysim{ii,1}.y-trusol(ysim{ii,1}.t)),'LineWidth',2); grid on; hold on;
end
subplot(222); lgd = "$h=\ $"+string(h_vec(:));
% legend(char(lgd(1)),char(lgd(2)),char(lgd(3)),char(lgd(4)),char(lgd(5))); 
figtex(gca); xlabel('$t$ [s]'); ylabel('$y_n(t)$');
subplot(224); 
% legend(char(lgd(1)),char(lgd(2)),char(lgd(3)),char(lgd(4)),char(lgd(5))); 
figtex(gca); xlabel('$t$ [s]'); ylabel('$|y_n - y(t)|$')




        
%         subplot(212); 
%         semilogy(tsim,abs(yf-trusol(tsim)),'LineWidth',2); grid on; hold on; figtex(gca);
        

%% Solving of equations

%%%%% Time step 1
h = 1e-6; % time step
t = ts:h:tf; % time vector 
%%% Euler forward method
tic; [yf,Ny] = ef_med_Nd(@(t,y) -10*y,1,t,h); toc ;

figure(10); clf;
subplot(411); plot(t,yf,t,exp(-10*t)); grid on;
ylabel('$y$'); xlabel('$t$'); 
lgd1 = char("$h=" + string(h) + "$"); lgd2 = char("$h=\infty$"); 
legend(lgd1,lgd2); figtex(gca,1); 
title('Simulation Results')

subplot(412); semilogy(t(2:end),Ny); grid on;
figtex(gca); ylabel('$\mathcal{N}_n(y)_ny(t)$'); xlabel('$t$'); 
title('Local Trucation Error')

subplot(413); plot(t,abs(yf-exp(-10*t))); grid on; hold on;
scatter(t(end),abs(yf(end)-exp(-10*t(end)))); hold off;
fprintf("Global error: %e \n",abs(yf(end)-exp(-10*t(end))));
figtex(gca); ylabel('$y_n - y(t)$'); xlabel('$t$'); 

%% Solving of equations

%%%%% Time step 1
h = 1e-5; % time step
t = ts:h:tf; % time vector 
%%% Euler forward method
tic; [yf,Ny] = ef_med_Nd(@(t,y) -1*y,1,t,h); toc ;
%%% Rk 45 method 
tic; sol_ode45 = ode45(@(t,y) -1*y,[ts tf],1,h); toc ;

figure(11); clf;

subplot(221); plot(t,yf,t,exp(-1*t),'--'); grid on;
ylabel('$y$'); xlabel('$t$'); 
lgd1 = char("numerical"); lgd2 = char("analytical");  
legend(lgd1,lgd2); figtex(gca,1); title('Euler forward')

subplot(223); semilogy(t,abs(yf-exp(-1*t))); grid on; hold on;
% scatter(t(end),abs(yf(end)-exp(-1*t(end)))); hold off;
% fprintf("Global error Euler forward: %f \n",abs(yf(end)-exp(-10*t(end))));
% fprintf("Global error ode45: %f \n",abs(sol_ode45.y(end)-exp(-10*sol_ode45.x(end))));
figtex(gca); ylabel('$y_n - y(t)$'); xlabel('$t$'); 

subplot(222); plot(sol_ode45.x,sol_ode45.y,sol_ode45.x,exp(-1*sol_ode45.x),'--'); grid on;
ylabel('$y$'); xlabel('$t$'); 
lgd1 = char("numerical"); lgd2 = char("analytical");  
legend(lgd1,lgd2); figtex(gca,1); title('ode45')

subplot(224); semilogy(sol_ode45.x,abs(sol_ode45.y-exp(-1*sol_ode45.x))); grid on; hold on;
% scatter(t(end),abs(yf(end)-exp(-1*t(end)))); hold off;
% fprintf("Global error Euler forward: %f \n",abs(yf(end)-exp(-10*t(end))));
% fprintf("Global error ode45: %f \n",abs(sol_ode45.y(end)-exp(-10*sol_ode45.x(end))));
figtex(gca); ylabel('$y_n - y(t)$'); xlabel('$t$'); 
