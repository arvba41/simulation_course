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

v_pk = 800; % peak voltage

ts = 0; tf = 0.05; % simulation start and stop time

w_ref = params.Nr*pi/30*params.np*0.98; % electrical speed of the em 
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
fy = @(t,y) [1/params.Ld*(vd(t)-params.Rs*y(1)+params.Lq*y(2)*w_ref);...
    1/params.Lq*(vq(t)-params.Rs*y(2)-(params.Ld*y(1) + params.psi)*w_ref);...
    ];

%% fixed step method

% varying the step-size
hstart = 2; hend =9;
h_vec = 10.^(-1*(hstart:hend));

% initial conditions
y0 = [0;0]; % 0 current start

for ii=1:length(h_vec)
    % simulation
        tsim = ts:h_vec(ii):tf;
        tStart = tic; 
        yf = ef(@(t,y) fy(t,y),y0,tsim,h_vec(ii)); 
        tEnd = toc(tStart);
    % data manupulation
        fprintf('time taken for h = %g is %.4f s.\n',...
            h_vec(ii),tEnd);
    % data storage
        ysimef{ii}.y = yf; ysimef{ii}.t = tsim;
    % clean up
        clear yf tsim is Te wout
end
%% plotting

textwidth = 15; textheight = 21;
figsize = [textwidth,  textheight*round(length(h_vec)/2)/6];
f = figure(10); clf; clear ax

tiledlayout(round(length(h_vec)/2),2);

for ii = 1:length(h_vec)
    ax(ii) = nexttile;
    titletext = "Step length h=" + num2str(h_vec(ii));
    plot(ysimef{ii}.t,ysimef{ii}.y,'LineWidth',2); 
    box off; title(char(titletext)); clear titletext;
    ylabel('$i(t)$ [A]'); xlabel('$t$ [s]'); 
end

legend('$i_d$','$i_q$','Interpreter','latex');

linkaxes(ax,'xy'); 
h = findall(f,'Type','axes'); % An array if you have subplots
set(h, 'TickLabelInterpreter', 'latex')

set(f, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
set(f, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);
ax(1).YLim = [0 300];

print -dpdf ../../doc/proj/Figures/plot1.pdf

%% plotting 2

textwidth = 15; textheight = 25;
figsize = [textwidth,  textheight];
f = figure(11); clf; clear ax

for ii = 1:length(h_vec)-1
    
    % id plots
    ax(1) = subplot(121);
    semilogy(ysimef{ii}.t,abs(ysimef{ii}.y(1,:) ...
        - interpn(ysimef{end}.t,ysimef{end}.y(1,:),ysimef{ii}.t)),...
        'LineWidth',2,'DisplayName',char("h = " + num2str(h_vec(ii)))); 
    hold on; box off; title('$i_d$ for different step lengths $h$'); 
    ylabel('$\Delta i_d(t)$ [A]'); xlabel('$t$ [s]'); 
    
    % id plots
    ax(2) = subplot(122);
    titletext = "Step length h=" + num2str(h_vec(ii));
    semilogy(ysimef{ii}.t,abs(ysimef{ii}.y(1,:) ...
        - interpn(ysimef{end}.t,ysimef{end}.y(1,:),ysimef{ii}.t)),...
        'LineWidth',2,'DisplayName',char("h = " + num2str(h_vec(ii)))); 
    hold on; box off; title('$i_q$ for different step lengths $h$');
    ylabel('$\Delta i_q(t)$ [A]'); xlabel('$t$ [s]'); 
end

legend('Interpreter','latex');

linkaxes(ax,'xy'); 
h = findall(f,'Type','axes'); % An array if you have subplots
set(h, 'TickLabelInterpreter', 'latex')

set(f, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
set(f, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

print -dpdf ../../doc/proj/Figures/plot2.pdf