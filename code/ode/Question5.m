clear all
clc

h = 0.02;
ts = 0;
tf = 6;

t1 = ts:h:tf; 

% y' = [0 -1; 1 0]*y; y(0) = [1; 0];
y_value = @(time_vector) exp([-time_vector; time_vector]); % true solution

%%% Euler forward method
y1_f = euler_f_vec(@(y) [0 -1;1 0]*y,[1; 0],t1,h); 
%%% Euler backward method
y1_b = euler_b_vec(@(y) [0 -1;1 0]*y,[1; 0],t1,h); 
%%% Trapz method
y1_fb = euler_fb_vec(@(y) [0 -1;1 0]*y,[1; 0],t1,h);
%%% Exact solution
y1_exact = [cos(t1); sin(t1)]; 


%%
%%%% Plottings
figure(1)
clf
tiledlayout(4,3)

nexttile;
plot(y1_f(1,:),y1_f(2,:))
ylabel('$y_1$');xlabel('$y_2$')
figtex(gca);

nexttile;
plot(y1_b(1,:),y1_b(2,:))
ylabel('$y_1$');xlabel('$y_2$')
figtex(gca);

nexttile;
plot(y1_fb(1,:),y1_fb(2,:))
ylabel('$y_1$');xlabel('$y_2$')
figtex(gca);

nexttile;
plot(t1,y1_f,t1,y1_exact,'--')
legend('analytical','analytical','excact','excact','location','northwest')
grid on
figtex(gca,1)
ylabel('$$y$$')
% xlabel('time [s]')
title('Euler Forward')

nexttile;
plot(t1,[y1_b],t1,y1_exact,'--')
grid on
figtex(gca,[])
% ylabel('$$y$$')
% xlabel('time [s]')
title('Euler Backward')

nexttile;
plot(t1,[y1_fb],t1, y1_exact,'--')
grid on
figtex(gca,[])
% ylabel('$$y$$')
% xlabel('time [s]')
title('Trapizoidal')

% Error calalculations
nexttile;
semilogy(t1,abs(y1_exact - y1_f))
grid on
figtex(gca,[])
ylabel('$$|y_n - y(t)|$$')

nexttile;
semilogy(t1,abs(y1_exact - y1_b))
grid on
figtex(gca,[])

nexttile;
semilogy(t1,abs(y1_exact - y1_fb))
grid on
figtex(gca,[])

% Gradients calalculations
nexttile;
semilogy(t1,y1_f./[1; 0])
grid on
figtex(gca,[]);
ylabel('$$y_n/y_0$$')
xlabel('time [s]')

nexttile;
semilogy(t1,y1_b./[1; 0])
grid on
figtex(gca,[]);
xlabel('time [s]')

nexttile;
semilogy(t1,y1_fb./[1; 0])
grid on
figtex(gca,[]);
xlabel('time [s]')

% nexttile;
% semilogy(t1,abs(y1_exact - [y1_f; y1_b; y1_fb]))
% grid on
% figtex(gca,[])
% ylabel('$$|y_n - y(t)|$$')
% xlabel('time [s]')
% 
% nexttile;
% semilogy(t1,[y1_f; y1_b; y1_fb]./[1; 0])
% grid on
% figtex(gca,[])
% ylabel('$$y_n/y_0$$')
% xlabel('time [s]')